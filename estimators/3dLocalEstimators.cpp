/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file 3dLocalEstimators.cpp
 * @ingroup Tools
 * @author Jeremy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), Universite de Lyon, France
 * LAboratoire de MAthematiques - LAMA (CNRS, UMR 5807), Universite de Savoie, France
 *
 * @date 2012/06/20
 *
 * DGtal 3D curvature shape comparator
 *
 * This file is part of the DGtalTools library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "DGtal/base/Common.h"
#include "DGtal/base/Clock.h"

//shapes
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/kernel/sets/SetPredicate.h"

//Digitizer
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/geometry/curves/GridCurve.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"


//Estimators
#include "DGtal/geometry/curves/estimation/TrueLocalEstimatorOnPoints.h"

#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantGaussianCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/FunctorOnCells.h"
#include "DGtal/topology/CanonicSCellEmbedder.h"

#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/images/imagesSetsUtils/SimpleThresholdForegroundPredicate.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"


//Vol Export
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/images/ImageHelper.h"

using namespace DGtal;

template <typename Shape, typename KSpace, typename ConstIterator>
std::vector<double>
estimateTrueMeanCurvatureQuantity( Shape aShape, KSpace K, double h, ConstIterator it, ConstIterator it_end )
{
  typedef typename Shape::RealPoint RealPoint;
  CanonicSCellEmbedder<KSpace> midpoint(K); 

  std::vector<double> values;
  for ( ; it != it_end; ++it )
  {
    RealPoint A = midpoint( *it ) * h;
    A = aShape.nearestPoint( A, 0.01, 200, 0.1 * h );
    values.push_back( aShape.meanCurvature( A ));
  }
  return values;
}

template <typename Shape, typename KSpace, typename ConstIterator>
std::vector<double>
estimateTrueGaussianCurvatureQuantity( Shape aShape, KSpace K, double h, ConstIterator it, ConstIterator it_end )
{
  typedef typename Shape::RealPoint RealPoint;
  CanonicSCellEmbedder<KSpace> midpoint(K); 

  std::vector<double> values;
  for ( ; it != it_end; ++it )
  {
    RealPoint A = midpoint( *it ) * h;
    A = aShape.nearestPoint( A, 0.01, 200, 0.1 * h );
    values.push_back( aShape.gaussianCurvature( A ));
  }
  return values;
}

template <typename Estimator, typename ConstIterator>
std::vector<typename Estimator::Quantity>
estimateQuantity( Estimator & estimator, 
      ConstIterator it, ConstIterator it_end )
{
  std::vector<typename Estimator::Quantity> values;

  std::back_insert_iterator< std::vector< typename Estimator::Quantity > > valuesIterator( values );
  estimator.eval ( it, it_end, valuesIterator );

  return values;
}

template <typename Space, typename Shape>
bool
compareShapeEstimators( const std::string & name,
                        Shape & aShape,
                        double border_min[],
                        double border_max[],
                        double h,
                        double radius_kernel,
                        const bool export_vol,
                        const std::string & pathToSaveVolFile = "",
                        const std::string & nameVolFile = "" )
{
  typedef typename Space::RealPoint RealPoint;
  typedef typename Space::Integer Integer;
  typedef GaussDigitizer< Space, Shape > Digitizer;
  typedef KhalimskySpaceND< Space::dimension, Integer > KSpace;
  typedef LightImplicitDigitalSurface< KSpace, Digitizer > LightImplicitDigSurface;
  typedef HyperRectDomain< Space > Domain;

  typedef DigitalSurface< LightImplicitDigSurface > MyDigitalSurface;
  typedef typename MyDigitalSurface::ConstIterator ConstIterator;
  typedef typename KSpace::Surfel Surfel;

  typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
  typedef GraphVisitorRange< Visitor > VisitorRange;
  typedef typename VisitorRange::ConstIterator SurfelConstIterator;

  Digitizer dig;
  dig.attach( aShape );
  dig.init( RealPoint( border_min ), RealPoint( border_max ), h );
  Domain domain = dig.getDomain();

  typedef typename ImageSelector< Domain, unsigned int >::Type Image;
  Image image( domain );
  DGtal::imageFromRangeAndValue( domain.begin(), domain.end(), image );

  KSpace K;
  bool ok = K.init( domain.lowerBound(), domain.upperBound(), true );
  if ( ! ok )
  {
    std::cerr << "[compareShapeEstimators]" << " error in creating KSpace." << std::endl;
    return false;
  }

  try
  {
    // Extracts shape boundary
    SurfelAdjacency< KSpace::dimension > SAdj ( true );
    Surfel bel = Surfaces<KSpace>::findABel ( K, dig, 10000 );

    LightImplicitDigSurface LightImplDigSurf ( K, dig, SAdj, bel );
    MyDigitalSurface surf ( LightImplDigSurf );


    std::cout << "# range size = " << surf.size() << std::endl;
    std::cout << "# h = " << h << std::endl;

    // Estimations
    // True values
    std::cout << "# True values computation" << std::endl;

    VisitorRange range( new Visitor( surf, *surf.begin() ) );
    SurfelConstIterator abegin = range.begin();
    SurfelConstIterator aend = range.end();

    std::vector<double> trueMeanCurvatures =
      estimateTrueMeanCurvatureQuantity ( aShape, K, h, abegin, aend );

    VisitorRange range2( new Visitor( surf, *surf.begin() ) );
    abegin = range2.begin();
    aend = range2.end();

    std::vector<double> trueGaussianCurvatures =
      estimateTrueGaussianCurvatureQuantity ( aShape, K, h, abegin, aend );

    typedef ImageToConstantFunctor< Image, Digitizer > MyPointFunctor;
    MyPointFunctor pointFunctor( &image, &dig, 1, true );

    // Integral Invariant Mean Curvature
    std::cout << "# Mean Curvature estimation from the Integral Invariant viewpoint" << std::endl;
    typedef FunctorOnCells< MyPointFunctor, KSpace > MyCellFunctor;

    MyCellFunctor functor ( pointFunctor, K );
    IntegralInvariantMeanCurvatureEstimator< KSpace, MyCellFunctor > IIMeanCurvatureEstimator ( K, functor );

    double re_convolution_kernel = radius_kernel * std::pow( h, 1.0/3.0 ); // to obtains convergence results, re must follow the rule re=kh^(1/3)
    std::cout << "# computed kernel radius = " << re_convolution_kernel << std::endl;

    Clock c;
    c.startClock();
    IIMeanCurvatureEstimator.init ( h, re_convolution_kernel );

    VisitorRange range3( new Visitor( surf, *surf.begin() ) );
    abegin = range3.begin();
    aend = range3.end();

    std::vector<double> IIMeanCurvatures =
      estimateQuantity( IIMeanCurvatureEstimator, abegin, aend );
    double TIIMeanCurv = c.stopClock();

    // Integral Invariant Gaussian Curvature
    std::cout << "# Gaussian Curvature estimation from the Integral Invariant viewpoint" << std::endl;

    IntegralInvariantGaussianCurvatureEstimator< KSpace, MyCellFunctor > IIGaussianCurvatureEstimator ( K, functor );
    std::cout << "# computed kernel radius = " << re_convolution_kernel << std::endl;

    c.startClock();
    IIGaussianCurvatureEstimator.init( h, re_convolution_kernel );

    VisitorRange range4( new Visitor( surf, *surf.begin() ) );
    abegin = range4.begin();
    aend = range4.end();

    std::vector<double> IIGaussianCurvatures =
      estimateQuantity( IIGaussianCurvatureEstimator, abegin, aend );
    double TIIGaussCurv = c.stopClock();

    // Output
    std::cout << "# Shape = " << name << std::endl
              << "# Time-IIMeanCurvature = " << TIIMeanCurv << std::endl
              << "# Time-IIGaussianCurvature = " << TIIGaussCurv << std::endl
              << "# id x y z TrueMeanCurvature TrueGaussianCurvature"
              << " IIMeanCurvature IIGaussianCurvature"
              << std::endl;

    VisitorRange range5( new Visitor( surf, *surf.begin() ) );
    abegin = range5.begin();
    aend = range5.end();

    Z3i::DigitalSet set3d( domain );

    for ( int i = 0; abegin != aend; ++abegin, ++i )
    {
      Surfel p = *abegin;
      std::cout << i << std::setprecision( 15 )
                << " " << p.myCoordinates[ 0 ] << " " << p.myCoordinates[ 1 ] << " " << p.myCoordinates[ 2 ]
                << " " << trueMeanCurvatures[ i ]
                << " " << trueGaussianCurvatures[ i ]
                << " " << IIMeanCurvatures[ i ]
                << " " << IIGaussianCurvatures[ i ]
                << std::endl;
      set3d.insert ( K.sCoords ( p ));
    }

    if ( export_vol )
    {
      typedef typename ImageSelector < Domain, unsigned char>::Type Image;
      Image mImage = ImageFromSet< Image >::create( set3d, 128, true);

      for ( typename Image::Range::Iterator it = mImage.begin(), itend = mImage.end();
            it != itend;
            ++it )
      {
        (*it) += 128;
      }

      std::stringstream sstm;
      sstm << pathToSaveVolFile << nameVolFile << "_h_" << h << ".vol";
      VolWriter< Image >::exportVol(sstm.str(), mImage );
    }
  }
  catch ( InputException e )
  {
    std::cerr << "[estimatorCurvatureComparator3D]"
              << " error."
              << e.what() << std::endl;
    return false;
  }
  return true;
}

void usage( int /*argc*/, char** argv )
{
  std::cerr << "Usage: " << argv[ 0 ] << " <Polynomial> <Px> <Py> <Pz> <Qx> <Qy> <Qz> <step> <radius> (<path_to_save>) (<name_of_vol_file>)" << std::endl;
  std::cerr << "\t - displays the boundary of a shape defined implicitly by a 3-polynomial <Polynomial>." << std::endl;
  std::cerr << "\t - P and Q defines the bounding box." << std::endl;
  std::cerr << "\t - step is the grid step." << std::endl;
  std::cerr << "\t - radius is the kernel support k radius." << std::endl;
  std::cerr << "\t - path is optional. It's the path where you want to generate a .vol file of the polynomial shape (for external computation). If no path was set, we don't export as a .vol file." << std::endl;
  std::cerr << "\t - name is optional. It's the name of your .vol file you want to generate (for external computation). If no name was set, we don't export as a .vol file." << std::endl;
}

int main( int argc, char** argv )
{
  if ( argc < 10 )
  {
    usage( argc, argv );
    return 1;
  }
  double border_min[ 3 ];
  double border_max[ 3 ];
  for ( unsigned int i = 0; i < 3; ++i )
  {
    border_min[ i ] = atof( argv[ 2+i ] );
    border_max[ i ] = atof( argv[ 5+i ] );
  }
  double h = atof( argv[ 8 ] );
  double radius = atof( argv[ 9 ] );

  typedef Z3i::Space::RealPoint RealPoint;
  typedef Z3i::Space::RealPoint::Coordinate Ring;
  typedef MPolynomial< 3, Ring > Polynomial3;
  typedef MPolynomialReader<3, Ring> Polynomial3Reader;
  typedef ImplicitPolynomial3Shape<Z3i::Space> ImplicitShape;

  /// Construction of the polynomial shape
  Polynomial3 poly;
  Polynomial3Reader reader;
  std::string poly_str = argv[ 1 ];
  std::string::const_iterator iter = reader.read( poly, poly_str.begin(), poly_str.end() );
  if ( iter != poly_str.end() )
  {
    std::cerr << "ERROR: I read only <"
              << poly_str.substr( 0, iter - poly_str.begin() )
              << ">, and I built P=" << poly << std::endl;
    return 1;
  }


  bool export_vol = false;
  std::string pathToSaveVolFile = "";
  std::string nameVolFile = "";
  if( argc >= 12 )
  {
    export_vol = true;
    pathToSaveVolFile = argv[ 10 ];
    nameVolFile = argv[ 11 ];
  }

  ImplicitShape shape( poly );

  /// Computation of 3D curvature estimators (mean & Gaussian) on the implicit shape
  compareShapeEstimators< Z3i::Space, ImplicitShape >
      (
        poly_str,
        shape,
        border_min, border_max,
        h,
        radius,
        export_vol, pathToSaveVolFile, nameVolFile
      );
}
