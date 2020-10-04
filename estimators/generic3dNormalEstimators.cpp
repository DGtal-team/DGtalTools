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
 * @file genericNormalEstimator.cpp
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2014/05/02
 *
 * Estimates the normal vector field of an implicitly defined shape
 * for several estimators. The implicit shape is perturbated by a
 * Kanungo noise.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <sstream>

#include "CLI11.hpp"

#include "DGtal/base/Common.h"
#include "DGtal/base/CountedPtr.h"
#include "DGtal/base/CountedConstPtrOrConstPtr.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/math/Statistic.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/geometry/surfaces/estimation/CNormalVectorEstimator.h"
#include "DGtal/geometry/surfaces/estimation/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/geometry/surfaces/estimation/VCMDigitalSurfaceLocalEstimator.h"
#include "DGtal/geometry/surfaces/estimation/TrueDigitalSurfaceLocalEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantCovarianceEstimator.h"
#include "DGtal/geometry/volumes/KanungoNoise.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/ShapeGeometricFunctors.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/SimpleThresholdForegroundPredicate.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#ifdef WITH_VISU3D_QGLVIEWER
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/Display3DFactory.h"
#endif

using namespace std;
using namespace DGtal;



/**
 @page generic3dNormalEstimators generic3dNormalEstimators
 
 @brief  computes a normal vector field over a digitized 3D implicit surface for several estimators.

 @b Usage:  ./estimators/generic3dNormalEstimators -p <polynomial> [options]

Computes a normal vector field over a digitized 3D implicit surface
for several estimators (II|VCM|Trivial|True), specified with -e. You
may add Kanungo noise with option -N. These estimators are compared
with ground truth. You may then: 1) visualize the normals or the angle
deviations with -V (if WITH_QGL_VIEWER is enabled), 2) outputs them as
a list of cells/estimations with -n, 3) outputs them as a ImaGene file
with -O, 4) outputs them as a NOFF file with -O, 5) computes
estimation statistics with option -S.

 
 @b Allowed @b options @b are : 

 @code
 
 -h,--help                             Print this help message and exit
 -p,--polynomial TEXT REQUIRED         the implicit polynomial whose zero-level defines the shape of interest.
 -N,--noise FLOAT=0                    the Kanungo noise level l=arg, with l^d the probability that a point at distance d is flipped inside/outside.
 -a,--minAABB FLOAT=-10                the min value of the AABB bounding box (domain)
 -A,--maxAABB FLOAT=10                 the max value of the AABB bounding box (domain)
 -g,--gridstep FLOAT=1                 the gridstep that defines the digitization (often called h).
 -e,--estimator TEXT:{True,VCM,II,Trivial}=True
                                       the chosen normal estimator: True | VCM | II | Trivial
 -R,--R-radius FLOAT=5                 the constant for parameter R in R(h)=R h^alpha (VCM).
 -r,--r-radius FLOAT=3                 the constant for parameter r in r(h)=r h^alpha (VCM,II,Trivial).
 -k,--kernel TEXT=hat                  the function chi_r, either hat or ball.
 --alpha FLOAT=0                       the parameter alpha in r(h)=r h^alpha (VCM).
 -t,--trivial-radius FLOAT=3           the parameter t defining the radius for the Trivial estimator. Also used for reorienting the VCM.
 -E,--embedding INT=0                  the surfel -> point embedding for VCM estimator: 0: Pointels, 1: InnerSpel, 2: OuterSpel.
 -o,--output TEXT=output               the output basename. All generated files will have the form <arg>-*, for instance <arg>-angle-deviation-<gridstep>.txt, <arg>-normals-<gridstep>.txt, <arg>-cells-<gridstep>.txt, <arg>-noff-<gridstep>.off.
 -S,--angle-deviation-stats            computes angle deviation error and outputs them in file <basename>-angle-deviation-<gridstep>.txt, as specified by -o <basename>.
 -x,--export TEXT=None                 exports surfel normals which can be viewed with ImaGene tool 'viewSetOfSurfels' in file <basename>-cells-<gridstep>.txt, as specified by -o <basename>. Parameter <arg> is None|Normals|AngleDeviation. The color depends on the angle deviation in degree: 0 metallic blue, 5 light cyan, 10 light green, 15 light yellow, 20 yellow, 25 orange, 30 red, 35, dark red, 40- grey
 -n,--normals                          outputs every surfel, its estimated normal, and the ground truth normal in file <basename>-normals-<gridstep>.txt, as specified by -o <basename>.
 -O,--noff                             exports the digital surface with normals as NOFF file <basename>-noff-<gridstep>.off, as specified by -o <basename>..
 -V,--view TEXT                        view the digital surface with normals.  Parameter <arg> is None|Normals|AngleDeviation. The color depends on the angle deviation in degree: 0 metallic blue, 5 light cyan, 10 light green, 15 light yellow, 20 yellow, 25 orange, 30 red, 35, dark red, 40- grey.
 @endcode


@b Example @b of @b implicit @b surface (specified by -p):
 - ellipse  : 90-3*x^2-2*y^2-z^2 
 - torus    : -1*(x^2+y^2+z^2+6*6-2*2)^2+4*6*6*(x^2+y^2) 
 - rcube    : 6561-x^4-y^4-z^4
 - goursat  : 8-0.03*x^4-0.03*y^4-0.03*z^4+2*x^2+2*y^2+2*z^2
 - distel   : 10000-(x^2+y^2+z^2+1000*(x^2+y^2)*(x^2+z^2)*(y^2+z^2))
 - leopold  : 100-(x^2*y^2*z^2+4*x^2+4*y^2+3*z^2)
 - diabolo  : x^2-(y^2+z^2)^2
 - heart    : -1*(x^2+2.25*y^2+z^2-1)^3+x^2*z^3+0.1125*y^2*z^3
 - crixxi   : -0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3


@b Implemented @b estimators (specified by -e):
 - True     : supposed to be the ground truth for computations. Of course, it is only approximations.
 - VCM      : normal estimator by digital Voronoi Covariance Matrix. Radii parameters are given by -R, -r.
 - II       : normal estimator by Integral Invariants. Radius parameter is given by -r.
 - Trivial  : the normal obtained by average trivial surfel normals in a ball neighborhood. Radius parameter is given by -r.

@note :
     - This is a normal *direction* evaluator more than a normal vector evaluator. Orientations of normals are deduced from ground truth. This is due to the fact that II and VCM only estimates normal directions.
     - This tool only analyses one surface component, and one that contains at least as many surfels as the width of the digital bounding box. This is required when analysing noisy data, where a lot of the small components are spurious. The drawback is that you cannot analyse the normals on a surface with several components.



 @b Example: 

 - @b Example @b of @b normal @b comparisons:
 You can estimate the normal vectors and compare the error with the true normals (option  -V AngleDeviation):
  @code
 # apply the Integral Invariant estimator  (options  -e II  -r 8 ):
 $ generic3dNormalEstimators -p "10000-(x^2+y^2+z^2+1000*(x^2+y^2)*(x^2+z^2)*(y^2+z^2))" -a -10 -A 10 -e II  -r 8  -V AngleDeviation -g 0.05 
 # apply the Voronoi Covariance Measure based estimator  (options  -e VCM  -r 8 ): 
 $ generic3dNormalEstimators -p "10000-(x^2+y^2+z^2+1000*(x^2+y^2)*(x^2+z^2)*(y^2+z^2))"  -a -10 -A 10 -e VCM   -R 20 -r 8 -V AngleDeviation -g 0.05 
 # to display source digital surface  :  visualize the normals and tape Key E to export surface (in  exportedMesh.off) 
 $ generic3dNormalEstimators -p "10000-(x^2+y^2+z^2+1000*(x^2+y^2)*(x^2+z^2)*(y^2+z^2))"  -a -10 -A 10 -e VCM   -R 20 -r 8 -V Normals -g 0.05 
 # visualize generated file:
 $ meshViewer -i  exportedMesh.off

  @endcode
  You should obtain such a result:
  |type                   |  visualization                                 | 
  | :----:                | :-------------:                                |
  | digital surface       | ![ ](resGeneric3dNormalEstimatorsDiscrete.png) |
  |  VCM angle-deviation  | ![ ](resGeneric3dNormalEstimatorsVCM.png)      |
  | II angle-deviation    | ![ ](resGeneric3dNormalEstimatorsII.png)       |




 - @b Example @b of @b normal @b visualization with noise add:
 This tool allows to add some noise on initial shape (option -N) and it is also possiblt to  visualize the source shape by using the resulting normals:
 @code 
 # apply the Integral Invariant estimator  (options  -e II  -r 6 ): 
  generic3dNormalEstimators -p "100-(x^2*y^2*z^2+4*x^2+4*y^2+3*z^2)"  -a -10 -A 10 -e II  -r 6 -V Normals -g 0.1 -N 0.3
 # apply the Voronoi Covariance Measure based estimator  (options  -e VCM  -r 8 ):  
   generic3dNormalEstimators -p "100-(x^2*y^2*z^2+4*x^2+4*y^2+3*z^2)"  -a -10 -A 10 -e VCM   -R 10 -r 5 -V Normals -g 0.1 -N 0.3
 # apply the true normals:
 ./estimators/generic3dNormalEstimators -p "100-(x^2*y^2*z^2+4*x^2+4*y^2+3*z^2)"  -a -10 -A 10 -e True  -V Normals -g 0.1 -N 0.3
 # to display source digital surface  :  tape Key E on previous displays to export surface (in  exportedMesh.off) 
 # visualize generated file:
 $ meshViewer -i  exportedMesh.off

 @endcode 

 You should obtain such a result:
  |type                   |  visualization                                 | 
  | :----:                | :-------------:                                |
  | digital surface       | ![ ](resGeneric3dNormalEstimatorsNoiseDiscrete.png) |
  |  VCM estimator  | ![ ](resGeneric3dNormalEstimatorsNoiseVCM.png)      |
  | II estimator    | ![ ](resGeneric3dNormalEstimatorsNoiseII.png)       |
  | True Normals    | ![ ](resGeneric3dNormalEstimatorsNoiseTrue.png)       |




 @see
 @ref generic3dNormalEstimators.cpp



 

 */

struct AllParams{
  double noise {0.0};
  double minAABB {-10.0};
  double maxAABB {10.0};
  double gridstep {1.0};
  double R_radius {5.0};
  double r_radius {3.0};
  double alpha {0.0};
  double trivial_radius {3.0};
  int embedding {0};
  std::string polynomials;
  std::string estimator {"True"};
  std::string kernel {"hat"};
  std::string output {"output"};
  std::string exportX {"None"};
  std::string view {"None"};
  bool angle_deviation_stats;
  bool normals;
  bool noff;
};



template <typename SCell, typename RealVector>
struct GradientMapAdapter {
  typedef std::map<SCell,RealVector> SCell2RealVectorMap;
  typedef SCell                                 Argument;
  typedef RealVector                               Value;
  GradientMapAdapter( ConstAlias<SCell2RealVectorMap> map )
    : myMap( map ) {}
  RealVector operator()( const Argument& arg ) const
  {
    typename SCell2RealVectorMap::const_iterator it = myMap->find( arg );
    if ( it != myMap->end() ) return it->second;
    else return RealVector();
  }
  CountedConstPtrOrConstPtr<SCell2RealVectorMap> myMap;
};

template <typename SCellEmbedder>
struct SCellEmbedderWithNormal : public SCellEmbedder
{
  using SCellEmbedder::space;
  using SCellEmbedder::operator();
  typedef typename SCellEmbedder::KSpace          KSpace;
  typedef typename SCellEmbedder::SCell            SCell;
  typedef typename SCellEmbedder::RealPoint    RealPoint;
  typedef SCell                                 Argument;
  typedef RealPoint                                Value;
  typedef typename KSpace::Space::RealVector  RealVector;
  typedef std::map<SCell,RealVector> SCell2RealVectorMap;
  typedef GradientMapAdapter<SCell,RealVector> GradientMap;

  SCellEmbedderWithNormal( ConstAlias<SCellEmbedder> embedder,
                           ConstAlias<SCell2RealVectorMap> map )
    : SCellEmbedder( embedder ), myMap( map )
  {}

  GradientMap gradientMap() const
  {
    return GradientMap( myMap );
  }

  CountedConstPtrOrConstPtr<SCell2RealVectorMap> myMap;
};

template <typename DigitalSurface,
          typename Estimator>
void exportNOFFSurface( const DigitalSurface& surface,
                        const Estimator& estimator,
                        std::ostream& output )
{
  typedef typename DigitalSurface::KSpace KSpace;
  typedef typename DigitalSurface::ConstIterator ConstIterator;
  typedef typename DigitalSurface::Surfel Surfel;
  typedef typename Estimator::Quantity Quantity;
  const KSpace& ks = surface.container().space();
  std::map<Surfel,Quantity> normals;
  for ( ConstIterator it = surface.begin(), itE = surface.end(); it != itE; ++it )
    {
      Quantity n_est = estimator.eval( it );
      normals[ *it ] = n_est;
    }
  CanonicSCellEmbedder<KSpace> surfelEmbedder( ks );
  typedef SCellEmbedderWithNormal< CanonicSCellEmbedder<KSpace> > Embedder;
  Embedder embedder( surfelEmbedder, normals );
  surface.exportAs3DNOFF( output, embedder );
}


/**
   Computes the normal estimations. Outputs statistics or export cell geometry.
 */
template <typename KSpace,
          typename ImplicitShape,
          typename Surface,
          typename TrueEstimator,
          typename Estimator>
void computeEstimation
( const AllParams &params,     //< command-line parameters
  const KSpace& K,                 //< cellular grid space
  const ImplicitShape& shape,      //< implicit shape "ground truth"
  const Surface& surface,          //< digital surface approximating shape
  TrueEstimator& true_estimator,   //< "ground truth" estimator
  Estimator& estimator )           //< an initialized estimator
{
  typedef typename Surface::Surfel Surfel;
  typedef typename Estimator::Quantity Quantity;
  typedef double Scalar;
  typedef DepthFirstVisitor< Surface > Visitor;
  typedef GraphVisitorRange< Visitor > VisitorRange;
  typedef typename VisitorRange::ConstIterator VisitorConstIterator;

  std::string fname = params.output;
  string nameEstimator = params.estimator;
  trace.beginBlock( "Computing " + nameEstimator + "estimations." );
  CountedPtr<VisitorRange> range( new VisitorRange( new Visitor( surface, *(surface.begin()) )) );
  std::vector<Quantity> n_estimations;
  estimator.eval( range->begin(), range->end(), std::back_inserter( n_estimations ) );
  trace.info() << "- nb estimations  = " << n_estimations.size() << std::endl;
  trace.endBlock();

  trace.beginBlock( "Computing ground truth." );
  range = CountedPtr<VisitorRange>( new VisitorRange( new Visitor( surface, *(surface.begin()) )) );
  std::vector<Quantity> n_true_estimations;
  true_estimator.eval( range->begin(), range->end(), std::back_inserter( n_true_estimations ) );
  trace.info() << "- nb estimations  = " << n_true_estimations.size() << std::endl;
  trace.endBlock();

  trace.beginBlock( "Correcting orientations." );
  ASSERT( n_estimations.size() == n_true_estimations.size() );
  for ( unsigned int i = 0; i < n_estimations.size(); ++i )
    if ( n_estimations[ i ].dot( n_true_estimations[ i ] ) < 0 )
      n_estimations[ i ] = -n_estimations[ i ];
  trace.endBlock();

  DGtal::GradientColorMap<double> grad( 0.0, 40.0 );
  // 0 metallic blue, 5 light cyan, 10 light green, 15 light
  // yellow, 20 yellow, 25 orange, 30 red, 35, dark red, 40- grey
  grad.addColor( DGtal::Color( 128, 128, 255 ) ); // 0
  grad.addColor( DGtal::Color( 128, 255, 255 ) ); // 5
  grad.addColor( DGtal::Color( 128, 255, 128 ) ); // 10
  grad.addColor( DGtal::Color( 255, 255, 128 ) ); // 15
  grad.addColor( DGtal::Color( 255, 255, 0   ) ); // 20
  grad.addColor( DGtal::Color( 255, 128, 0   ) ); // 25
  grad.addColor( DGtal::Color( 255,   0, 0   ) ); // 30
  grad.addColor( DGtal::Color( 128,   0, 0   ) ); // 35
  grad.addColor( DGtal::Color( 128, 128, 128 ) ); // 40

  if ( params.angle_deviation_stats  )
    {
      trace.beginBlock( "Computing angle deviation error stats." );
      std::ostringstream adev_sstr;
      adev_sstr << fname << "-" << nameEstimator << "-angle-deviation-"
                << estimator.h() << ".txt";
      DGtal::Statistic<Scalar> adev_stat;
      unsigned int i = 0;
      range = CountedPtr<VisitorRange>( new VisitorRange( new Visitor( surface, *(surface.begin()) )) );
      for ( VisitorConstIterator it = range->begin(), itE = range->end(); it != itE; ++it, ++i )
        {
          Quantity n_est = n_estimations[ i ];
          Quantity n_true_est = n_true_estimations[ i ];
          Scalar angle_error = acos( n_est.dot( n_true_est ) );
          adev_stat.addValue( angle_error );
        }
      adev_stat.terminate();
      std::ofstream adev_output( adev_sstr.str().c_str() );
      adev_output << "# Average error X of the absolute angle between two vector estimations." << std::endl;
      adev_output << "# h L1 L2 Loo E[X] Var[X] Min[X] Max[X] Nb[X]" << std::endl;
      adev_output << estimator.h()
                  << " " << adev_stat.mean() // L1
                  << " " << sqrt( adev_stat.unbiasedVariance()
                                  + adev_stat.mean()*adev_stat.mean() ) // L2
                  << " " << adev_stat.max() // Loo
                  << " " << adev_stat.mean() // E[X] (=L1)
                  << " " << adev_stat.unbiasedVariance() // Var[X]
                  << " " << adev_stat.min() // Min[X]
                  << " " << adev_stat.max() // Max[X]
                  << " " << adev_stat.samples() // Nb[X]
                  << std::endl;
      adev_output.close();
      trace.endBlock();
    }
  if ( params.exportX != "None" )
    {
      trace.beginBlock( "Exporting cell geometry." );
      std::ostringstream export_sstr;
      export_sstr << fname << "-" << nameEstimator << "-cells-"
                  << estimator.h() << ".txt";
      std::ofstream export_output( export_sstr.str().c_str() );
      export_output << "# ImaGene viewer (viewSetOfSurfels) file format for displaying cells." << std::endl;
      bool adev =  params.exportX == "AngleDeviation";
      unsigned int i = 0;
      range = CountedPtr<VisitorRange>( new VisitorRange( new Visitor( surface, *(surface.begin()) )) );
      for ( VisitorConstIterator it = range->begin(), itE = range->end(); it != itE; ++it, ++i )
        {
          Quantity n_est = n_estimations[ i ];
          Quantity n_true_est = n_true_estimations[ i ];
          Scalar angle_error = acos( n_est.dot( n_true_est ) )*180.0 / 3.14159625;
          Surfel s = *it;
          export_output
            << "CellN"
            << " " << min( 1023, max( 512+K.sKCoord( s, 0 ), 0 ) )
            << " " << min( 1023, max( 512+K.sKCoord( s, 1 ), 0 ) )
            << " " << min( 1023, max( 512+K.sKCoord( s, 2 ), 0 ) )
            << " " << K.sSign( s );
          Color c = grad( 0 );
          if ( adev ) c = grad( max( 0.0, min( angle_error, 40.0 ) ) );
          export_output << " " << ((double) c.red() / 255.0 )
                        << " " << ((double) c.green() / 255.0 )
                        << " " << ((double) c.blue() / 255.0 );
          export_output << " " << n_est[ 0 ] << " " << n_est[ 1 ]
                        << " " << n_est[ 2 ] << std::endl;
        }
      export_output.close();
      trace.endBlock();
    }
  if ( params.normals )
    {
      trace.beginBlock( "Exporting cells normals." );
      std::ostringstream export_sstr;
      export_sstr << fname << "-" << nameEstimator << "-normals-"
                  << estimator.h() << ".txt";
      std::ofstream export_output( export_sstr.str().c_str() );
      export_output << "# kx ky kz sign n_est[0] n_est[1] n_est[2] n_true[0] n_true[1] n_true[2]" << std::endl;
      unsigned int i = 0;
      range = CountedPtr<VisitorRange>( new VisitorRange( new Visitor( surface, *(surface.begin()) )) );
      for ( VisitorConstIterator it = range->begin(), itE = range->end(); it != itE; ++it, ++i )
        {
          Quantity n_est = n_estimations[ i ];
          Quantity n_true_est = n_true_estimations[ i ];
          Surfel s = *it;
          export_output
            << K.sKCoord( s, 0 ) << " " << K.sKCoord( s, 1 ) << " " << K.sKCoord( s, 2 )
            << " " << K.sSign( s )
            << " " << n_est[ 0 ] << " " << n_est[ 1 ] << " " << n_est[ 2 ]
            << " " << n_true_est[ 0 ] << " " << n_true_est[ 1 ] << " " << n_true_est[ 2 ]
            << std::endl;
        }
      export_output.close();
      trace.endBlock();
    }
  if ( params.noff )
    {
      trace.beginBlock( "Exporting NOFF file." );
      std::ostringstream export_sstr;
      export_sstr << fname << "-" << nameEstimator << "-noff-"
                  << estimator.h() << ".off";
      std::ofstream export_output( export_sstr.str().c_str() );
      std::map<Surfel,Quantity> normals;
      unsigned int i = 0;
      range = CountedPtr<VisitorRange>( new VisitorRange( new Visitor( surface, *(surface.begin()) )) );
      for ( VisitorConstIterator it = range->begin(), itE = range->end(); it != itE; ++it, ++i )
        {
          Quantity n_est = n_estimations[ i ];
          normals[ *it ] = n_est;
        }
      CanonicSCellEmbedder<KSpace> surfelEmbedder( K );
      typedef SCellEmbedderWithNormal< CanonicSCellEmbedder<KSpace> > Embedder;
      Embedder embedder( surfelEmbedder, normals );
      surface.exportAs3DNOFF( export_output, embedder );
      export_output.close();
      trace.endBlock();
    }
#ifdef WITH_VISU3D_QGLVIEWER
  if ( params.exportX != "None" )
    {
      typedef typename KSpace::Space Space;
      typedef Viewer3D<Space,KSpace> MyViewever3D;
      typedef Display3DFactory<Space,KSpace> MyDisplay3DFactory;
      int argc = 1;
      char name[] = "Viewer";
      char* argv[ 1 ];
      argv[ 0 ] = name;
      Surfel s;
      QApplication application( argc, argv );
      MyViewever3D viewer( K );
      viewer.show();
      viewer << SetMode3D( s.className(), "Basic" );
      trace.beginBlock( "Viewing surface." );
      bool adev =  params.view == "AngleDeviation";

      unsigned int i = 0;
      range = CountedPtr<VisitorRange>( new VisitorRange( new Visitor( surface, *(surface.begin()) )) );
      for ( VisitorConstIterator it = range->begin(), itE = range->end(); it != itE; ++it, ++i )
        {
          Quantity n_est = n_estimations[ i ];
          Quantity n_true_est = n_true_estimations[ i ];
          Scalar angle_error = acos( n_est.dot( n_true_est ) )*180.0 / 3.14159625;
          s = *it;
          Color c = grad( 0 );
          if ( adev ) c = grad( max( 0.0, min( angle_error, 40.0 ) ) );
          viewer.setFillColor( c );
          MyDisplay3DFactory::drawOrientedSurfelWithNormal( viewer, s, n_est, false );
        }
      trace.endBlock();
      viewer << MyViewever3D::updateDisplay;
      application.exec();
    }
#endif

}

template <typename KSpace,
          typename ImplicitShape,
          typename Surface,
          typename KernelFunction,
          typename PointPredicate>
void chooseEstimator
( const AllParams &params,     //< command-line parameters
  const KSpace& K,                 //< cellular grid space
  const ImplicitShape& shape, //< implicit shape "ground truth"
  const Surface& surface,     //< digital surface approximating shape
  const KernelFunction& chi,  //< the kernel function
  const PointPredicate& ptPred )   //< analysed implicit digital shape as a PointPredicate
{
  string nameEstimator = params.estimator;
  double h = params.gridstep;
  typedef functors::ShapeGeometricFunctors::ShapeNormalVectorFunctor<ImplicitShape> NormalFunctor;
  typedef TrueDigitalSurfaceLocalEstimator<KSpace, ImplicitShape, NormalFunctor> TrueEstimator;
  TrueEstimator true_estimator;
  true_estimator.attach( shape );
  true_estimator.setParams( K, NormalFunctor(), 20, 0.1, 0.01 );
  true_estimator.init( h, surface.begin(), surface.end() );
  if ( nameEstimator == "True" )
    {
      trace.beginBlock( "Chosen estimator is: True." );
      typedef TrueDigitalSurfaceLocalEstimator<KSpace, ImplicitShape, NormalFunctor> Estimator;
      Estimator estimator;
      estimator.attach( shape );
      estimator.setParams( K, NormalFunctor(), 20, 0.1, 0.01 );
      estimator.init( h, surface.begin(), surface.end() );
      trace.endBlock();
      computeEstimation( params, K, shape, surface, true_estimator, estimator );
    }
  else if ( nameEstimator == "VCM" )
    {
      trace.beginBlock( "Chosen estimator is: VCM." );
      typedef typename KSpace::Space Space;
      typedef typename Surface::DigitalSurfaceContainer SurfaceContainer;
      typedef ExactPredicateLpSeparableMetric<Space,2> Metric;
      typedef VoronoiCovarianceMeasureOnDigitalSurface<SurfaceContainer,Metric,
                                                       KernelFunction> VCMOnSurface;
      typedef functors::VCMNormalVectorFunctor<VCMOnSurface> NormalFunctor;
      typedef VCMDigitalSurfaceLocalEstimator<SurfaceContainer,Metric,
                                              KernelFunction, NormalFunctor> VCMNormalEstimator;
      int embedding = params.embedding;
      Surfel2PointEmbedding embType = embedding == 0 ? Pointels :
                                      embedding == 1 ? InnerSpel : OuterSpel;
      double R = params.R_radius;
      double r = params.r_radius;
      double t = params.trivial_radius;
      double alpha = params.alpha;
      if ( alpha != 0.0 ) R *= pow( h, alpha-1.0 );
      if ( alpha != 0.0 ) r *= pow( h, alpha-1.0 );
      trace.info() << "- R=" << R << " r=" << r << " t=" << t << std::endl;
      VCMNormalEstimator estimator;
      estimator.attach( surface );
      estimator.setParams( embType, R, r, chi, t, Metric(), true );
      estimator.init( h, surface.begin(), surface.end() );
      trace.endBlock();
      computeEstimation( params, K, shape, surface, true_estimator, estimator );
    }
  else if ( nameEstimator == "II" )
    {
      trace.beginBlock( "Chosen estimator is: II." );
      typedef typename KSpace::Space Space;
      typedef HyperRectDomain<Space> Domain;
      typedef ImageContainerBySTLVector< Domain, bool> Image;
      typedef typename Domain::ConstIterator DomainConstIterator;
      typedef functors::SimpleThresholdForegroundPredicate<Image> ThresholdedImage;
      typedef functors::IINormalDirectionFunctor<Space> IINormalFunctor;
      typedef IntegralInvariantCovarianceEstimator<KSpace, ThresholdedImage, IINormalFunctor> IINormalEstimator;
      double r = params.r_radius;
      double alpha = params.alpha;
      if ( alpha != 0.0 ) r *= pow( h, alpha-1.0 );
      trace.info() << " r=" << r << std::endl;
      trace.beginBlock( "Preparing characteristic set." );
      Domain domain( K.lowerBound(), K.upperBound() );
      Image image( domain );
      for ( DomainConstIterator it = domain.begin(), itE = domain.end(); it != itE; ++it )
        {
          image.setValue( *it, ptPred( *it ) );
        }
      trace.endBlock();
      trace.beginBlock( "Initialize II estimator." );
      ThresholdedImage thresholdedImage( image, false );
      IINormalEstimator ii_estimator( K, thresholdedImage );
      ii_estimator.setParams( r );
      ii_estimator.init( h, surface.begin(), surface.end() );
      trace.endBlock();
      trace.endBlock();
      computeEstimation( params, K, shape, surface, true_estimator, ii_estimator );
   }

}

template <typename KSpace,
          typename ImplicitShape,
          typename Surface,
          typename PointPredicate>
void chooseKernel
( const AllParams& params,     //< command-line parameters
  const KSpace& K,                 //< cellular grid space
  const ImplicitShape& shape,      //< implicit shape "ground truth"
  const Surface& surface,          //< digital surface approximating shape
  const PointPredicate& ptPred )   //< analysed implicit digital shape as a PointPredicate
{
  string kernel = params.kernel;
  double h = params.gridstep;
  double r = params.r_radius;
  double alpha = params.alpha;
  if ( alpha != 0.0 ) r *= pow( h, alpha-1.0 );
  if ( kernel == "hat" ) {
    typedef typename KSpace::Point Point;
    typedef functors::HatPointFunction<Point,double> KernelFunction;
    KernelFunction chi_r( 1.0, r );
    chooseEstimator( params, K, shape, surface, chi_r, ptPred );
  } else if ( kernel == "ball" ) {
    typedef typename KSpace::Point Point;
    typedef functors::BallConstantPointFunction<Point,double> KernelFunction;
    KernelFunction chi_r( 1.0, r );
    chooseEstimator( params, K, shape, surface, chi_r, ptPred );
  }
}

template <typename KSpace,
          typename ImplicitShape,
          typename ImplicitDigitalShape >
int chooseSurface
( const AllParams &params, //< command-line parameters
  const KSpace& K,                 //< cellular grid space
  const ImplicitShape& shape,      //< implicit shape "ground truth"
  const ImplicitDigitalShape& dshape ) //< analysed implicit digital shape
{
  // Selecting a model of surface depending on noise / not noise.
  typedef double Scalar;
  if ( params.noise == 0.0 )
    { // no noise
      trace.beginBlock( "Make digital surface..." );
      typedef LightImplicitDigitalSurface<KSpace,ImplicitDigitalShape> SurfaceContainer;
      typedef DigitalSurface< SurfaceContainer > Surface;
      typedef typename Surface::Surfel Surfel;
      SurfelAdjacency< KSpace::dimension > surfAdj( true );
      Surfel bel;
      try {
        bel = Surfaces<KSpace>::findABel( K, dshape, 10000 );
      } catch (DGtal::InputException e) {
        trace.error() << "ERROR Unable to find bel." << std::endl;
        return 3;
      }
      SurfaceContainer* surfaceContainer = new SurfaceContainer( K, dshape, surfAdj, bel );
      CountedPtr<Surface> ptrSurface( new Surface( surfaceContainer ) ); // acquired
      trace.info() << "- surface component has " << ptrSurface->size() << " surfels." << std::endl;
      trace.endBlock();
      chooseKernel( params, K, shape, *ptrSurface, dshape );
    }
  else
    { // noise
      trace.beginBlock( "Make digital surface..." );
      typedef typename ImplicitDigitalShape::Domain Domain;
      typedef KanungoNoise< ImplicitDigitalShape, Domain > KanungoPredicate;
      typedef LightImplicitDigitalSurface< KSpace, KanungoPredicate > SurfaceContainer;
      typedef DigitalSurface< SurfaceContainer > Surface;
      typedef typename Surface::Surfel Surfel;
      SurfelAdjacency< KSpace::dimension > surfAdj( true );
      Surfel bel;
      const Domain shapeDomain = dshape.getDomain();
      KanungoPredicate* noisified_dshape = new KanungoPredicate( dshape, shapeDomain, params.noise );
      // We have to search for a big connected component.
      CountedPtr<Surface> ptrSurface;
      double minsize = dshape.getUpperBound()[0] - dshape.getLowerBound()[0];
      unsigned int nb_surfels = 0;
      unsigned int tries = 0;
      do {
        try { // Search initial bel
          bel = Surfaces<KSpace>::findABel( K, *noisified_dshape, 10000 );
        } catch (DGtal::InputException e) {
          trace.error() << "ERROR Unable to find bel." << std::endl;
          return 3;
        }
        SurfaceContainer* surfaceContainer = new SurfaceContainer( K, *noisified_dshape, surfAdj, bel );
        ptrSurface = CountedPtr<Surface>( new Surface( surfaceContainer ) ); // acquired
        nb_surfels = ptrSurface->size();
      } while ( ( nb_surfels < 2 * minsize ) && ( tries++ < 150 ) );
      if( tries >= 150 )
        {
          trace.error() << "ERROR cannot find a proper bel in a big enough component." << std::endl;
          return 4;
        }
      trace.info() << "- surface component has " << nb_surfels << " surfels." << std::endl;
      trace.endBlock();
      chooseKernel( params, K, shape, *ptrSurface, *noisified_dshape );
    }
  return 0;
}


///////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
  
  AllParams allParams;
  
  // parse command line ----------------------------------------------
  CLI::App app;
  std::stringstream ss_descr;
  ss_descr << "Computes a normal vector field over a digitized 3D implicit surface for several estimators (II|VCM|Trivial|True), specified with -e. You may add Kanungo noise with option -N. These estimators are compared with ground truth. You may then: 1) visualize the normals or the angle deviations with -V (if WITH_QGL_VIEWER is enabled), 2) outputs them as a list of cells/estimations with -n, 3) outputs them as a ImaGene file with -O, 4) outputs them as a NOFF file with -O, 5) computes estimation statistics with option -S.";
  
  ss_descr << "Example of use:\n"
       << "./generic3dNormalEstimator -p \"90-3*x^2-2*y^2-z^2\" -o VCM-ellipse -a -10 -A 10 -e VCM -R 3 -r 3 -t 2 -E 0 -x Normals" << endl << endl
       << "Example of implicit surfaces (specified by -p):" << endl
       << " - ellipse  : 90-3*x^2-2*y^2-z^2 " << endl
       << " - torus    : -1*(x^2+y^2+z^2+6*6-2*2)^2+4*6*6*(x^2+y^2) " << endl
       << " - rcube    : 6561-x^4-y^4-z^4" << endl
       << " - goursat  : 8-0.03*x^4-0.03*y^4-0.03*z^4+2*x^2+2*y^2+2*z^2" << endl
       << " - distel   : 10000-(x^2+y^2+z^2+1000*(x^2+y^2)*(x^2+z^2)*(y^2+z^2))" << endl
       << " - leopold  : 100-(x^2*y^2*z^2+4*x^2+4*y^2+3*z^2)" << endl
       << " - diabolo  : x^2-(y^2+z^2)^2" << endl
       << " - heart    : -1*(x^2+2.25*y^2+z^2-1)^3+x^2*z^3+0.1125*y^2*z^3" << endl
       << " - crixxi   : -0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3" << endl << endl;
  ss_descr << "Implemented estimators (specified by -e):" << endl
       << " - True     : supposed to be the ground truth for computations. Of course, it is only approximations." << endl
       << " - VCM      : normal estimator by digital Voronoi Covariance Matrix. Radii parameters are given by -R, -r." << endl
       << " - II       : normal estimator by Integral Invariants. Radius parameter is given by -r." << endl
       << " - Trivial  : the normal obtained by average trivial surfel normals in a ball neighborhood. Radius parameter is given by -r." << endl
       << endl;
  ss_descr << "Note:" << endl
       << "     - This is a normal *direction* evaluator more than a normal vector evaluator. Orientations of normals are deduced from ground truth. This is due to the fact that II and VCM only estimates normal directions." << endl
       << "     - This tool only analyses one surface component, and one that contains at least as many surfels as the width of the digital bounding box. This is required when analysing noisy data, where a lot of the small components are spurious. The drawback is that you cannot analyse the normals on a surface with several components." << endl;
  
  app.description(ss_descr.str());
  
  app.add_option("--polynomial,-p",allParams.polynomials, "the implicit polynomial whose zero-level defines the shape of interest.") ->required();
  app.add_option("--noise,-N", allParams.noise,"the Kanungo noise level l=arg, with l^d the probability that a point at distance d is flipped inside/outside.", true );
  app.add_option("--minAABB,-a", allParams.minAABB, "the min value of the AABB bounding box (domain)", true );
  app.add_option("--maxAABB,-A", allParams.maxAABB, "the max value of the AABB bounding box (domain)", true );
  app.add_option("--gridstep,-g", allParams.gridstep, "the gridstep that defines the digitization (often called h).", true);
  app.add_option("--estimator,-e", allParams.estimator, "the chosen normal estimator: True | VCM | II | Trivial", true)
  -> check(CLI::IsMember({"True", "VCM", "II", "Trivial"}));
  app.add_option("--R-radius,-R", allParams.R_radius, "the constant for parameter R in R(h)=R h^alpha (VCM).",true);
  app.add_option("--r-radius,-r", allParams.r_radius, "the constant for parameter r in r(h)=r h^alpha (VCM,II,Trivial).", true);
  app.add_option("--kernel,-k", allParams.kernel, "the function chi_r, either hat or ball.", true);
  app.add_option("--alpha", allParams.alpha, "the parameter alpha in r(h)=r h^alpha (VCM).", true);
  app.add_option("--trivial-radius,-t", allParams.trivial_radius, "the parameter t defining the radius for the Trivial estimator. Also used for reorienting the VCM.", true);
  app.add_option("--embedding,-E",allParams.embedding, "the surfel -> point embedding for VCM estimator: 0: Pointels, 1: InnerSpel, 2: OuterSpel.", true);
  app.add_option("--output,-o", allParams.output, "the output basename. All generated files will have the form <arg>-*, for instance <arg>-angle-deviation-<gridstep>.txt, <arg>-normals-<gridstep>.txt, <arg>-cells-<gridstep>.txt, <arg>-noff-<gridstep>.off.", true);
  app.add_flag("--angle-deviation-stats,-S", allParams.angle_deviation_stats, "computes angle deviation error and outputs them in file <basename>-angle-deviation-<gridstep>.txt, as specified by -o <basename>.");

  app.add_option("--export,-x",allParams.exportX, "exports surfel normals which can be viewed with ImaGene tool 'viewSetOfSurfels' in file <basename>-cells-<gridstep>.txt, as specified by -o <basename>. Parameter <arg> is None|Normals|AngleDeviation. The color depends on the angle deviation in degree: 0 metallic blue, 5 light cyan, 10 light green, 15 light yellow, 20 yellow, 25 orange, 30 red, 35, dark red, 40- grey", true );
  app.add_flag("--normals,-n", allParams.normals, "outputs every surfel, its estimated normal, and the ground truth normal in file <basename>-normals-<gridstep>.txt, as specified by -o <basename>.");
  app.add_flag("--noff,-O", allParams.noff, "exports the digital surface with normals as NOFF file <basename>-noff-<gridstep>.off, as specified by -o <basename>..");
#ifdef WITH_VISU3D_QGLVIEWER
  app.add_option("--view,-V", allParams.view, "view the digital surface with normals.  Parameter <arg> is None|Normals|AngleDeviation. The color depends on the angle deviation in degree: 0 metallic blue, 5 light cyan, 10 light green, 15 light yellow, 20 yellow, 25 orange, 30 red, 35, dark red, 40- grey.");
#endif

  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  
  trace.beginBlock( "Make implicit shape..." );
  typedef Z3i::Space Space;
  typedef double Scalar;
  typedef MPolynomial< 3, Scalar > Polynomial3;
  typedef MPolynomialReader<3, Scalar> Polynomial3Reader;
  typedef ImplicitPolynomial3Shape<Space> ImplicitShape;
  Polynomial3 poly;
  Polynomial3Reader reader;
  string::const_iterator iter = reader.read( poly, allParams.polynomials.begin(), allParams.polynomials.end() );
  if ( iter != allParams.polynomials.end() )
    {
      trace.error() << "ERROR reading polynomial: I read only <" << allParams.polynomials.substr( 0, iter - allParams.polynomials.begin() )
                    << ">, and I built P=" << poly << std::endl;
      return 2;
    }
  CountedPtr<ImplicitShape> shape( new ImplicitShape( poly ) ); // smart pointer
  trace.endBlock();

  trace.beginBlock( "Make implicit digital shape..." );
  typedef Z3i::KSpace KSpace;
  typedef Space::RealPoint RealPoint;
  typedef GaussDigitizer< Space, ImplicitShape > ImplicitDigitalShape;
  typedef ImplicitDigitalShape::Domain Domain;
  Scalar min_x = allParams.minAABB;
  Scalar max_x = allParams.maxAABB;
  Scalar h = allParams.gridstep;
  RealPoint p1( min_x, min_x, min_x );
  RealPoint p2( max_x, max_x, max_x );
  CountedPtr<ImplicitDigitalShape> dshape( new ImplicitDigitalShape() );
  dshape->attach( *shape );
  dshape->init( p1, p2, h );
  Domain domain = dshape->getDomain();
  KSpace K;
  K.init( domain.lowerBound(), domain.upperBound(), true );
  trace.info() << "- domain is " << domain << std::endl;
  trace.endBlock();

  chooseSurface( allParams, K, *shape, *dshape );

  return 0;
}

