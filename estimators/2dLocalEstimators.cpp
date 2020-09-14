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
 * @file 2dLocalEstimators.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr)
 * LIRIS (CNRS, UMR 5205),
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS,
 * France
 * @author Jeremy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), Universite de Lyon, France
 * LAboratoire de MAthematiques - LAMA (CNRS, UMR 5807), Universite de Savoie, France
 *
 * @date 2011/07/04
 *
 * DGtal tangeant & curvature estimators comparator.
 * @WARNING IntegralInvariant curvature results are set in the reverse order in file. You need to reverse the order in order to compare with others.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>

#include "CLI11.hpp"

#include "DGtal/base/Common.h"
#include "DGtal/base/Clock.h"

//shapes
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/DigitalSurface.h"

//Digitizer
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/geometry/curves/GridCurve.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/geometry/volumes/KanungoNoise.h"


//Estimators
#include "DGtal/geometry/curves/estimation/TrueLocalEstimatorOnPoints.h"
#include "DGtal/geometry/curves/estimation/TrueGlobalEstimatorOnPoints.h"
#include "DGtal/geometry/curves/estimation/ParametricShapeCurvatureFunctor.h"
#include "DGtal/geometry/curves/estimation/ParametricShapeTangentFunctor.h"
#include "DGtal/geometry/curves/estimation/ParametricShapeArcLengthFunctor.h"

#include "DGtal/geometry/curves/BinomialConvolver.h"
#include "DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h"
#include "DGtal/geometry/curves/estimation/SegmentComputerEstimators.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/StabbingCircleComputer.h"

#include "DGtal/images/ImageHelper.h"
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h"

#include "DGtal/kernel/BasicPointFunctors.h"

using namespace DGtal;



/**
 @page Doc2dLocalEstimators 2dLocalEstimators
 
 @brief Compares local estimators on implicit shapes using DGtal library.

 @b Usage: 	2dlocalEstimators --output <output> --shape <shapeName> [required parameters] --estimators <binaryWord> --properties <binaryWord>

Below are the different available families of estimators: 
	 - True estimators
	 - Maximal DSS based estimators
	 - Maximal DCA based estimators
	 - Binomial convolver based estimators
	 - Integral Invariants based estimators

The i-th family of estimators is enabled if the i-th character of the binary word is not 0. The default binary word is '10000'. This means that the first family of estimators, ie. true estimators, is enabled, whereas the next ones are disabled. 
Below are the different available properties: 
	 - Tangent
	 - Curvature




 @b Allowed @b options @b are : 
 @code
  -h,--help                             Print this help message and exit
  -l,--list                             List all available shapes
  -o,--output TEXT REQUIRED             Output
  -s,--shape TEXT REQUIRED              Shape name
  -R,--radius FLOAT                     Radius of the shape
  -K,--kernelradius FLOAT=2.35162e-314  Radius of the convolution kernel (Integral invariants estimators)
  --alpha FLOAT=0.333333                Alpha parameter for Integral Invariant computation
  -A,--axis1 FLOAT                      Half big axis of the shape (ellipse)
  -a,--axis2 FLOAT                      Half small axis of the shape (ellipse)
  -r,--smallradius FLOAT=5              Small radius of the shape
  -v,--varsmallradius FLOAT=5           Variable small radius of the shape
  -k UINT=3                             Number of branches or corners the shape (default 3)
  --phi FLOAT=0                         Phase of the shape (in radian)
  -w,--width FLOAT=10                   Width of the shape
  -p,--power FLOAT=2                    Power of the metric
  -x,--center_x FLOAT=0                 x-coordinate of the shape center
  -y,--center_y FLOAT=0                 y-coordinate of the shape center
  -g,--gridstep FLOAT=1                 Gridstep for the digitization
  -n,--noise FLOAT=0                    Level of noise to perturb the shape
  --properties TEXT=11                  the i-th property is disabled iff there is a 0 at position i
  -e,--estimators TEXT=1000             the i-th estimator is disabled iff there is a 0 at position i
  -E,--exportShape TEXT                 Exports the contour of the source shape as a sequence of discrete points (.sdp)
  --lambda BOOLEAN=0                    Use the shape to get a better approximation of the surface (optional)
 @endcode

 @b Example: 
 With this tool you can easely compare several estimator with the real value:
 @code
$  2dlocalEstimators --output curvature --shape flower --radius 15 -v 5  --gridstep 1  --estimators 11100 --properties 01
 @endcode

You can display the result by using gnuplot:

@code
$ gnuplot
gnuplot> plot [][-0.4:0.35] 'curvature_True_curvature.dat' w lines title "true curvature" , 'curvature_MDCA_curvature.dat' w lines  title "Maximal DCA curvature estimator", 'curvature_MDSSl_curvature.dat' w lines title "Maximal DSS based estimators"
@endcode


You should obtain such a graph:


 @image html res2dLocalEstimators.png "Resulting visualization."
 
 @see
 @ref 2dLocalEstimators.cpp

 */



/**
 * Global vectors to describe the available shapes and their
 * parameters.
 */
std::vector<std::string> shapes2D;
std::vector<std::string> shapesDesc;
std::vector<std::string> shapesParam1;
std::vector<std::string> shapesParam2;
std::vector<std::string> shapesParam3;
std::vector<std::string> shapesParam4;

template< typename RealPoint >
struct OptionsIntegralInvariant
{
  double alpha; // <! Alpha parameter for the convolution kernel. 1/3 by default
  double radius; // <! Radius of the convolution kernel.
  RealPoint center; // <! Center of the shape.
  bool lambda_optimized;
};


/**
 * Create the static list of shapes.
 *
 */
void createList()
{
  shapes2D.push_back("ball");
  shapesDesc.push_back("Ball for the Euclidean metric.");
  shapesParam1.push_back("--radius [-R]");
  shapesParam2.push_back("");
  shapesParam3.push_back("");
  shapesParam4.push_back("");

  shapes2D.push_back("square");
  shapesDesc.push_back("square (no signature).");
  shapesParam1.push_back("--width [-w]");
  shapesParam2.push_back("");
  shapesParam3.push_back("");
  shapesParam4.push_back("");

  shapes2D.push_back("lpball");
  shapesDesc.push_back("Ball for the l_power metric (no signature).");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--power [-p]");
  shapesParam3.push_back("");
  shapesParam4.push_back("");

  shapes2D.push_back("flower");
  shapesDesc.push_back("Flower with k petals with radius ranging from R+/-v.");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--varsmallradius [-v],");
  shapesParam3.push_back("--k [-k],");
  shapesParam4.push_back("--phi");

  shapes2D.push_back("ngon");
  shapesDesc.push_back("Regular k-gon.");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--k [-k],");
  shapesParam3.push_back("--phi");
  shapesParam4.push_back("");

  shapes2D.push_back("accflower");
  shapesDesc.push_back("Accelerated Flower with k petals.");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--varsmallradius [-v],");
  shapesParam3.push_back("--k [-k],");
  shapesParam4.push_back("--phi");

  shapes2D.push_back("ellipse");
  shapesDesc.push_back("Ellipse.");
  shapesParam1.push_back("--axis1 [-A],");
  shapesParam2.push_back("--axis2 [-a],");
  shapesParam3.push_back("--phi");
  shapesParam4.push_back("");


}

/**
 * Display the shape list with parameters.
 *
 */
void displayList()
{
  trace.emphase()<<"2D Shapes:"<<std::endl;
  for(unsigned int i=0; i<shapes2D.size(); ++i)
    trace.info()<<"\t"<<shapes2D[i]<<"\t"
               <<shapesDesc[i]<<std::endl
              <<"\t\tRequired parameter(s): "
             << shapesParam1[i]<<" "
             << shapesParam2[i]<<" "
             << shapesParam3[i]<<" "
             << shapesParam4[i]<<std::endl;

}


/**
 * Check if a given shape is available. If not, we exit with an error.
 * If it is, we return the corresponding index in the global vectors.
 *
 * @param shapeName name of the shape to search.
 *
 * @return index of the shape in the shape vectors.
 */
unsigned int checkAndReturnIndex(const std::string &shapeName)
{
  unsigned int pos=0;

  while ((pos < shapes2D.size()) && (shapes2D[pos] != shapeName))
    pos++;

  if (pos == shapes2D.size())
  {
    trace.error() << "The specified shape has not found.";
    trace.info() << std::endl;
    exit(1);
  }

  return pos;
}

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam(std::string param)
{
  trace.error() <<" Parameter: "<<param<<" is required.";
  trace.info()<<std::endl;
  exit(1);
}

/**
 * Estimation error message.
 *
 * @param currentSize number of values returned by the estimator
 * @param expectedSize expected number of values
 */
void estimationError(int currentSize, int expectedSize)
{
  if (currentSize != expectedSize)
  {
    trace.error() << " error in the estimation"
                  << " (got " << currentSize << " values"
                  << " instead of " << expectedSize << ")";
    trace.info() << std::endl;
    exit(1);
  }

}

/**
 * Estimation. Merely call the init and eval methods of the
 * given estimator.
 *
 * @param estimator any local estimator
 * @param h the grid step
 * @param itb begin iterator
 * @param ite end iterator
 * @param ito output iterator on estimated quantities
 */
template <typename Estimator, typename ConstIterator, typename OutputIterator>
void
estimation( Estimator & estimator, double h,
            const ConstIterator& itb, const ConstIterator& ite, const OutputIterator& ito )
{
  Clock c;
  c.startClock();
  estimator.init( h, itb, ite );
  estimator.eval( itb, ite, ito );
  double time = c.stopClock();
  std::cout << "# Time: " << time << std::endl;
}


/**
 *
 * @return Euclidean radius for the convolver of Integral Invariant estimators
 */
template< typename ConstIteratorOnPoints, typename RPoint >
unsigned int suggestedSizeIntegralInvariant( const double h,
                                             const RPoint& center,
                                             const ConstIteratorOnPoints& itb,
                                             const ConstIteratorOnPoints& ite )
{
  ConstIteratorOnPoints it = itb;
  RPoint p( *it );
  RPoint distance = p - center;
  auto minRadius = distance.norm();
  ++it;

  for ( ; it != ite; ++it )
  {
    p = *it;
    distance = p - center;
    if ( distance.norm() < minRadius )
    {
      minRadius = distance.norm();
    }
  }

  return minRadius * h;
}


/**
 * Estimation of tangents and curvature
 * from several different methods
 *
 * @param filename name of a file to save results ( will be postfix by name of estimators )
 * @param aShape shape
 * @param h grid step
 * @param optionsII options for Integral Invariants estimators
 * @param options estimators to use (BC, II, MDCA, ...)
 * @param properties properties of estimators (curvature and/or tangeant)
 * @param noiseLevel level to noised the shape. 0 <= noiseLevel < 1
 */
template <typename Space, typename Shape>
bool
computeLocalEstimations( const std::string & filename,
                         Shape * aShape,
                         const double & h,
                         struct OptionsIntegralInvariant< Z2i::RealPoint > optionsII,
                         const std::string & options,
                         const std::string & properties,
                         const std::string & outShape,
                         double noiseLevel = 0.0 )
{
  // Types
  typedef typename Space::Vector Vector;
  typedef typename Space::RealPoint RealPoint;
  typedef typename Space::Integer Integer;
  typedef HyperRectDomain<Space> Domain;
  typedef KhalimskySpaceND<Space::dimension,Integer> KSpace;
  typedef typename KSpace::SCell SCell;
  typedef GaussDigitizer<Space,Shape> Digitizer;
  typedef KanungoNoise< Digitizer, Z2i::Domain > KanungoPredicate;

  bool withNoise = ( noiseLevel <= 0.0 ) ? false : true;
  /*if( withNoise )
        noiseLevel = std::pow(noiseLevel, h);*/

  ASSERT (( noiseLevel < 1.0 ));

  bool tangent = ( properties.at( 0 ) != '0' ) ? true : false;
  bool curvature = ( properties.at( 1 ) != '0' ) ? true : false;

  // Digitizer
  Digitizer* dig = new Digitizer();
  dig->attach( *aShape ); // attaches the shape.
  Vector vlow(-1,-1); Vector vup(1,1);
  dig->init( aShape->getLowerBound()+vlow, aShape->getUpperBound()+vup, h );
  Domain domain = dig->getDomain();

  //Noise

  Clock c;

  // Create cellular space
  KSpace K;
  bool ok = K.init( dig->getLowerBound(), dig->getUpperBound(), true );
  if ( ! ok )
  {
    std::cerr << "[2dLocalEstimators]"
              << " error in creating KSpace." << std::endl;
    return false;
  }
  try {

    // Extracts shape boundary
    SurfelAdjacency< KSpace::dimension > SAdj( true );
    SCell bel;
    std::vector< SCell > points;

    KanungoPredicate  *noisifiedObject;
    if ( withNoise )
    {
      noisifiedObject = new KanungoPredicate( *dig, domain, noiseLevel );
      bel = Surfaces< KSpace >::findABel( K, *noisifiedObject, 10000 );
      Surfaces< KSpace >::track2DBoundary( points, K, SAdj, *noisifiedObject, bel );

      double minsize = dig->getUpperBound()[0] - dig->getLowerBound()[0];
      while( points.size() < 2 * minsize )
      {
        points.clear();
        bel = Surfaces< KSpace >::findABel( K, *noisifiedObject, 10000 );
        Surfaces< KSpace >::track2DBoundary( points, K, SAdj, *noisifiedObject, bel );
      }
    }
    else
    {
      bel = Surfaces< KSpace >::findABel( K, *dig, 10000 );
      Surfaces< KSpace >::track2DBoundary( points, K, SAdj, *dig, bel );
    }

    // Create GridCurve
    GridCurve< KSpace > gridcurve;
    gridcurve.initFromSCellsVector( points );
    if(outShape != "")
      {
        std::ofstream outS;
        outS.open(outShape.c_str());
        for(const auto &p : points)
          {

            Dimension track = *( K.sDirs( p ) );
            SCell pointel = K.sIndirectIncident( p, track );
            outS << K.sCoords( pointel )[0] << " " << K.sCoords( pointel )[1] << std::endl;
          }
        outS.close();
      }
    // Ranges
    typedef typename GridCurve< KSpace >::MidPointsRange PointsRange;
    PointsRange pointsRange = gridcurve.getMidPointsRange();

    // Estimations
    if (gridcurve.isClosed())
    {
      if (options.at(0) != '0')
      {
        if( tangent )
        {
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_tangeant.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# True tangents computation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< RealPoint > out_it( file, "\n" );

          typedef ParametricShapeTangentFunctor< Shape > TangentFunctor;
          typedef typename PointsRange::ConstCirculator C;
          TrueLocalEstimatorOnPoints< C, Shape, TangentFunctor >
              trueTangentEstimator;
          trueTangentEstimator.attach( aShape );
          estimation( trueTangentEstimator, h,
                      pointsRange.c(), pointsRange.c(),
                      out_it );

          file.close();

        }

        if( curvature )
        {
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_curvature.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# True curvature computation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< double > out_it( file, "\n" );

          typedef ParametricShapeCurvatureFunctor< Shape > CurvatureFunctor;
          typedef typename PointsRange::ConstCirculator C;
          TrueLocalEstimatorOnPoints< C, Shape, CurvatureFunctor >
              trueCurvatureEstimator;
          trueCurvatureEstimator.attach( aShape );
          estimation( trueCurvatureEstimator, h,
                      pointsRange.c(), pointsRange.c(),
                      out_it );

          file.close();
        }
      }

      // Maximal Segments
      if (options.at(1) != '0')
      {
        if( tangent )
        {
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MDSS_tangeant.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Most centered maximal DSS tangent estimation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< RealPoint > out_it( file, "\n" );

          typedef typename GridCurve< KSpace >::PointsRange PointsRange2;
          PointsRange2 pointsRange2 = gridcurve.getPointsRange();

          typedef typename PointsRange2::ConstCirculator C;
          typedef ArithmeticalDSSComputer< C, int, 4 > SegmentComputer;
          typedef TangentFromDSSEstimator<SegmentComputer> SCFunctor;
          SegmentComputer sc;
          SCFunctor f;


          MostCenteredMaximalSegmentEstimator<SegmentComputer,SCFunctor> MDSSTangentEstimator(sc, f);
          estimation( MDSSTangentEstimator, h,
                      pointsRange2.c(), pointsRange2.c(),
                      out_it );

          file.close();
        }
        if( curvature )
        {
          c.startClock();

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MDSSl_curvature.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Most centered maximal DSS (length) curvature estimation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< double > out_it( file, "\n" );

          typedef typename GridCurve< KSpace >::PointsRange PointsRange2;
          PointsRange2 pointsRange2 = gridcurve.getPointsRange();

          typedef typename PointsRange2::ConstCirculator C;
          typedef ArithmeticalDSSComputer< C, int, 4 > SegmentComputer;
          typedef CurvatureFromDSSLengthEstimator<SegmentComputer> SCFunctor;
          SegmentComputer sc;
          SCFunctor f;
          MostCenteredMaximalSegmentEstimator<SegmentComputer,SCFunctor> MDSSCurvatureEstimator(sc, f);

          estimation( MDSSCurvatureEstimator, h,
                      pointsRange2.c(), pointsRange2.c(),
                      out_it );

          file.close();


          memset(&full_filename[0], 0, sizeof(full_filename));
          sprintf( full_filename, "%s%s", filename.c_str(), "_MDSSlw_curvature.dat" );
          file.open( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Most centered maximal DSS (length & width) curvature estimation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< double > out_it2( file, "\n" );

          typedef CurvatureFromDSSEstimator<SegmentComputer> SCFunctor2;
          SegmentComputer sc2;
          SCFunctor2 f2;
          MostCenteredMaximalSegmentEstimator<SegmentComputer,SCFunctor2> MDSSCurvatureEstimator2(sc2, f2);
          estimation( MDSSCurvatureEstimator2, h,
                      pointsRange2.c(), pointsRange2.c(),
                      out_it2 );

          double time = c.stopClock();
          file << "# Time: " << time << std::endl;

          file.close();

        }
      }

      //Maximal circular arcs
      if (options.at(2) != '0')
      {
        if( tangent )
        {
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MDCA_tangent.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Most centered maximal DCA tangents estimation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< RealPoint > out_it( file, "\n" );

          typedef typename GridCurve<KSpace>::IncidentPointsRange Range;
          typedef typename Range::ConstCirculator C;
          Range r = gridcurve.getIncidentPointsRange();
          typedef StabbingCircleComputer<C> SegmentComputer;
          typedef TangentFromDCAEstimator<SegmentComputer> SCFunctor;
          SegmentComputer sc;
          SCFunctor f;
          MostCenteredMaximalSegmentEstimator<SegmentComputer,SCFunctor> MDCATangentEstimator(sc, f);
          estimation( MDCATangentEstimator, h,
                      r.c(), r.c(),
                      out_it );

          file.close();
        }

        if( curvature )
        {
          c.startClock();

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MDCA_curvature.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Most centered maximal DCA curvature estimation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          std::ostream_iterator< double > out_it( file, "\n" );

          typedef typename GridCurve<KSpace>::IncidentPointsRange Range;
          typedef typename Range::ConstCirculator C;
          Range r = gridcurve.getIncidentPointsRange();
          typedef StabbingCircleComputer<C> SegmentComputer;
          typedef CurvatureFromDCAEstimator<SegmentComputer, false> SCFunctor;
          SegmentComputer sc;
          SCFunctor f;
          MostCenteredMaximalSegmentEstimator<SegmentComputer,SCFunctor> MDCACurvatureEstimator(sc, f);
          estimation( MDCACurvatureEstimator, h,
                      r.c(), r.c(),
                      out_it );

          double time = c.stopClock();
          file << "# Time: " << time << std::endl;

          file.close();
        }
      }

      //Binomial convolver
      if (options.at(3) != '0')
      {
        if( tangent )
        {
          c.startClock();

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_BC_tangeant.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Tangents estimation from binomial convolution" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          typedef typename PointsRange::ConstIterator I;
          typedef BinomialConvolver<I, double> MyBinomialConvolver;
          file << "# mask size = " <<
                  MyBinomialConvolver::suggestedSize( h, pointsRange.begin(), pointsRange.end() ) << std::endl;

          typedef TangentFromBinomialConvolverFunctor< MyBinomialConvolver, RealPoint >
              TangentBCFct;
          BinomialConvolverEstimator< MyBinomialConvolver, TangentBCFct> BCTangentEstimator;

          std::ostream_iterator< RealPoint > out_it( file, "\n" );

          BCTangentEstimator.init( h, pointsRange.begin(), pointsRange.end(), true );
          BCTangentEstimator.eval( pointsRange.begin(), pointsRange.end(), out_it );

          double time = c.stopClock();
          file << "# Time: " << time << std::endl;

          file.close();
        }

        if( curvature )
        {
          c.startClock();

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_BC_curvature.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Curvature estimation from binomial convolution" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          typedef typename PointsRange::ConstIterator I;
          typedef BinomialConvolver<I, double> MyBinomialConvolver;
          file << "# mask size = " <<
                  MyBinomialConvolver::suggestedSize( h, pointsRange.begin(), pointsRange.end() ) << std::endl;

          std::ostream_iterator< double > out_it( file, "\n" );

          typedef CurvatureFromBinomialConvolverFunctor< MyBinomialConvolver, double >
              CurvatureBCFct;
          BinomialConvolverEstimator< MyBinomialConvolver, CurvatureBCFct> BCCurvatureEstimator;

          BCCurvatureEstimator.init( h, pointsRange.begin(), pointsRange.end(), true );
          BCCurvatureEstimator.eval( pointsRange.begin(), pointsRange.end(), out_it );

          double time = c.stopClock();
          file << "# Time: " << time << std::endl;

          file.close();
        }
      }

      /// <! @WARNING IntegralInvariant curvature results are set in the reverse order in file. You need to reverse the order in order to compare with others.
      //Integral Invariants
      if (options.at(4) != '0')
      {
        if( curvature )
        {
          c.startClock();

          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_II_curvature.dat" );
          std::ofstream file( full_filename );

          file << "# h = " << h << std::endl;
          file << "# Integral Invariant curvature estimation" << std::endl;
          file << "# range size = " << pointsRange.size() << std::endl;
          if ( withNoise )
          {
            file << "# noise level (init) = " << noiseLevel/h << std::endl;
            file << "# noise level (current) = " << noiseLevel << std::endl;
          }

          if( optionsII.radius <= 0.0 )
          {
            optionsII.radius = suggestedSizeIntegralInvariant( h,  optionsII.center, pointsRange.begin(), pointsRange.end() );
            file << "# Estimated radius: " << optionsII.radius << std::endl;
          }
          double re = optionsII.radius * std::pow( h, optionsII.alpha );
          file << "# full kernel (digital) size (with alpha = " << optionsII.alpha << ") = " <<
                  re / h << std::endl;

          std::ostream_iterator< double > out_it( file, "\n" );

          if ( withNoise )
          {
            typedef LightImplicitDigitalSurface< KSpace, KanungoPredicate > LightImplicitDigSurface;
            typedef DigitalSurface< LightImplicitDigSurface > DigSurface;

            LightImplicitDigSurface LightImplDigSurf( K, *noisifiedObject, SAdj, bel );
            DigSurface surf( LightImplDigSurf );

            typedef DepthFirstVisitor< DigSurface > Visitor;
            typedef GraphVisitorRange< Visitor > VisitorRange;
            typedef typename VisitorRange::ConstIterator VisitorConstIterator;

            VisitorRange range( new Visitor( surf, *surf.begin() ));
            VisitorConstIterator ibegin = range.begin();
            VisitorConstIterator iend = range.end();

            typedef functors::IICurvatureFunctor<Z2i::Space> MyIICurvatureFunctor;
            typedef IntegralInvariantVolumeEstimator< KSpace, KanungoPredicate, MyIICurvatureFunctor > MyIICurvatureEstimator;
            
            MyIICurvatureFunctor curvatureFunctor;
            curvatureFunctor.init( h, re );

            MyIICurvatureEstimator curvatureEstimator( curvatureFunctor );
            curvatureEstimator.attach( K, *noisifiedObject );
            curvatureEstimator.setParams( re/h );
            curvatureEstimator.init( h, ibegin, iend );

            curvatureEstimator.eval( ibegin, iend, out_it );
          }
          else
          {
            typedef LightImplicitDigitalSurface< KSpace, Digitizer > LightImplicitDigSurface;
            typedef DigitalSurface< LightImplicitDigSurface > DigSurface;

            LightImplicitDigSurface LightImplDigSurf( K, *dig, SAdj, bel );
            DigSurface surf( LightImplDigSurf );

            typedef DepthFirstVisitor< DigSurface > Visitor;
            typedef GraphVisitorRange< Visitor > VisitorRange;
            typedef typename VisitorRange::ConstIterator VisitorConstIterator;

            VisitorRange range( new Visitor( surf, *surf.begin() ));
            VisitorConstIterator ibegin = range.begin();
            VisitorConstIterator iend = range.end();

            typedef functors::IICurvatureFunctor<Z2i::Space> MyIICurvatureFunctor;
            typedef IntegralInvariantVolumeEstimator< KSpace, Digitizer, MyIICurvatureFunctor > MyIICurvatureEstimator;
            
            MyIICurvatureFunctor curvatureFunctor;
            curvatureFunctor.init( h, re );

            MyIICurvatureEstimator curvatureEstimator( curvatureFunctor );
            curvatureEstimator.attach( K, *dig );
            curvatureEstimator.setParams( re/h );
            curvatureEstimator.init( h, ibegin, iend );

            curvatureEstimator.eval( ibegin, iend, out_it );
          }

          double time = c.stopClock();
          file << "# Time: " << time << std::endl;

          file.close();
        }
      }

      //delete noisifiedObject;
      delete dig;
    }
    else
    {
      //delete noisifiedObject;
      delete dig;
      std::cerr << "[computeLocalEstimations]"
                << " error: open digital curve found." << std::endl;
      return false;
    }
  }
  catch ( InputException e )
  {
    std::cerr << "[computeLocalEstimations]"
              << " error in finding a bel." << std::endl;
    return false;
  }
  return true;
}


///////////////////////////////////////////////////////////////////////////////


int main( int argc, char** argv )
{
  // parse command line CLI ----------------------------------------------
  CLI::App app;
  std::string shapeName;
  std::string filename;
  double radius, kernelradius;
  double power {2.0};
  double smallradius {5};
  double varsmallradius {5};
  double cx {0.0}, cy {0.0};
  double h {1.0};
  unsigned int k {3};
  double phi {0.0};
  double width {10.0};
  double axis1, axis2;
  double alpha {1.0/3.0};
  double noiseLevel {0.0};
  std::string properties {"11"};
  std::string outShape {""};
  bool lambda {false};
  std::string options {"1000"};

  app.description("Compares local estimators on implicit shapes using DGtal library.\n Typical use example:\n \t 2dlocalEstimators --output <output> --shape <shapeName> [required parameters] --estimators <binaryWord> --properties <binaryWord>\n");
  auto listOpt = app.add_flag("--list,-l","List all available shapes");
  auto outputOpt = app.add_option("--output,-o", filename, "Output")->required();
  auto shapeNameOpt = app.add_option("--shape,-s", shapeName, "Shape name")->required();
  auto radiusOpt = app.add_option("--radius,-R", radius, "Radius of the shape" );
  auto kernelradiusOpt = app.add_option("--kernelradius,-K", radius, "Radius of the convolution kernel (Integral invariants estimators)", true);
  auto alphaOpt = app.add_option("--alpha", alpha, "Alpha parameter for Integral Invariant computation", true);
  auto axis1Opt = app.add_option("--axis1,-A", axis1, "Half big axis of the shape (ellipse)" );
  auto axis2Opt = app.add_option("--axis2,-a", axis2, "Half small axis of the shape (ellipse)" );
  auto smallradiusOpt = app.add_option("--smallradius,-r", smallradius, "Small radius of the shape", true);
  auto varsmallradiusOpt = app.add_option("--varsmallradius,-v", varsmallradius, "Variable small radius of the shape", true );
  auto kOpt = app.add_option("-k", k, "Number of branches or corners the shape (default 3)", true );
  auto phiOpt = app.add_option("--phi", phi, "Phase of the shape (in radian)", true );
  auto widthOpt = app.add_option("--width,-w", width, "Width of the shape", true );
  auto powerOpt = app.add_option("--power,-p", power, "Power of the metric", true );
  app.add_option("--center_x,-x", cx, "x-coordinate of the shape center", true );
  app.add_option("--center_y,-y", cy, "y-coordinate of the shape center", true );
  app.add_option("--gridstep,-g", h, "Gridstep for the digitization", true );
  app.add_option("--noise,-n", noiseLevel, "Level of noise to perturb the shape", true);
  app.add_option("--properties", properties, "the i-th property is disabled iff there is a 0 at position i", true);
  app.add_option("--estimators,-e", options, "the i-th estimator is disabled iff there is a 0 at position i", true);
  app.add_option("--exportShape,-E", outShape, "Exports the contour of the source shape as a sequence of discrete points (.sdp)", true);
  app.add_option("--lambda", lambda, "Use the shape to get a better approximation of the surface (optional)", true);

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI

  //List creation
  createList();

  if ( listOpt->count() > 0 )
  {
    displayList();
    return 0;
  }

  unsigned int nb = 4; //number of available methods
  if (options.size() < nb)
  {
    trace.error() << " At least " << nb
                  << " characters are required "
                  << " with option --estimators.";
    trace.info() << std::endl;
    exit(1);
  }

  nb = 2; //number of available properties
  if (properties.size() < nb)
  {
    trace.error() << " At least " << nb
                  << " characters are required "
                  << " with option --properties.";
    trace.info() << std::endl;
    exit(1);
  }

  //We check that the shape is known
  unsigned int id = checkAndReturnIndex(shapeName);

  // standard types
  typedef Z2i::Space Space;
  typedef Space::RealPoint RealPoint;

  RealPoint center( cx, cy );
  
  struct OptionsIntegralInvariant< RealPoint > optII;
  optII.radius = kernelradius;
  optII.alpha = alpha;
  optII.lambda_optimized = lambda;
  optII.center = center;

  if (id ==0)
  {
    if (radiusOpt->count()==0) missingParam("--radius");
    //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
    
    Ball2D<Space> * ball = new Ball2D<Space>( center, radius);
    computeLocalEstimations<Space>( filename, ball, h, optII, options, properties, outShape, noiseLevel );
    delete ball;
  }
  else if (id ==1)
  {
    //if (widthOpt->count()==0) missingParam("--width");
    
    ImplicitHyperCube<Space> object(Z2i::Point(0,0), width/2);
    trace.error()<< "Not available.";
    trace.info()<<std::endl;
  }
  else if (id ==2)
  {
    //if (powerOpt->count()==0) missingParam("--power");
    if (radiusOpt->count()==0) missingParam("--radius");
    
    ImplicitRoundedHyperCube<Space> ball( Z2i::Point(0,0), radius, power );
    trace.error()<< "Not available.";
    trace.info()<<std::endl;
  }
  else if (id ==3)
  {
    //if (varsmallradiusOpt->count()==0) missingParam("--varsmallradius");
    if (radiusOpt->count()==0) missingParam("--radius");
    //if (kOpt->count()==0) missingParam("--k");
    //if (phiOpt->count()==0) missingParam("--phi");
    //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
    
    Flower2D<Space> * flower = new Flower2D<Space>( center, radius, varsmallradius, k, phi );
    computeLocalEstimations<Space>( filename, flower, h, optII, options, properties, outShape, noiseLevel );
    delete flower;
  }
  else if (id ==4)
  {
    if (radiusOpt->count()==0) missingParam("--radius");
    //if (kOpt->count()==0) missingParam("--k");
    //if (phiOpt->count()==0) missingParam("--phi");
    //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
    
    NGon2D<Space> * object = new NGon2D<Space>( center, radius, k, phi );
    computeLocalEstimations<Space>( filename, object, h, optII, options, properties, outShape, noiseLevel );
    delete object;
  }
  else if (id ==5)
  {
    //if (varsmallradiusOpt->count()==0) missingParam("--varsmallradius");
    if (radiusOpt->count()==0) missingParam("--radius");
    //if (kOpt->count()==0) missingParam("--k");
    //if (phiOpt->count()==0) missingParam("--phi");
    //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
  
    AccFlower2D<Space> * accflower = new AccFlower2D<Space>( center, radius, varsmallradius, k, phi );
    computeLocalEstimations<Space>( filename, accflower, h, optII, options, properties, outShape, noiseLevel );
    delete accflower;
  }
  else if (id ==6)
  {
    if (axis1Opt->count()==0) missingParam("--axis1");
    if (axis2Opt->count()==0) missingParam("--axis2");
    //if (phiOpt->count()==0) missingParam("--phi");
    //if (!(vm.count("kernelradius"))) missingParam("--kernelradius");
    
    Ellipse2D<Space> * ellipse = new Ellipse2D<Space>( center, axis1, axis2, phi );
    computeLocalEstimations<Space>( filename, ellipse, h, optII, options, properties, outShape, noiseLevel );
    delete ellipse;
  }
}
