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
#include <string>

#include "DGtal/base/Common.h"
#include "DGtal/base/Clock.h"
#include "DGtal/helpers/StdDefs.h"

#include "CLI11.hpp"


//shapes
#include "DGtal/shapes/implicit/ImplicitBall.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"

//Digitizer
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/geometry/volumes/KanungoNoise.h"
#include "DGtal/topology/CanonicSCellEmbedder.h"


//Estimators
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/geometry/surfaces/FunctorOnCells.h"
#include "DGtal/geometry/volumes/distance/LpMetric.h"

#include "DGtal/geometry/curves/estimation/TrueLocalEstimatorOnPoints.h"
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantCovarianceEstimator.h"

#include "DGtal/geometry/surfaces/estimation/LocalEstimatorFromSurfelFunctorAdapter.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/MongeJetFittingMeanCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/MongeJetFittingGaussianCurvatureEstimator.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/MongeJetFittingPrincipalCurvaturesEstimator.h"

using namespace DGtal;
using namespace functors;



/**
 @page Doc3dLocalEstimators 3dLocalEstimators
 
 @brief Compares local estimators on implicit shapes using DGtal library.

 @b Usage:  3dLocalEstimators [options] --shape <shape> --h <h> --radius <radius> --estimators <binaryWord> --output <output>

    Below are the different available families of estimators:
    - Integral Invariant Mean
    - Integral Invariant Gaussian
    - Monge Jet Fitting Mean
    - Monge Jet Fitting Gaussian

 @b Allowed @b options @b are :
 @code
 -h,--help                             Print this help message and exit
 -s,--shape TEXT REQUIRED              Shape
 -o,--output TEXT=result.dat REQUIRED  Output file
 -r,--radius FLOAT REQUIRED            Kernel radius for IntegralInvariant
 --alpha FLOAT=0.333333                Alpha parameter for Integral Invariant computation
 --h FLOAT REQUIRED                    Grid step
 -a,--minAABB FLOAT=-10                Min value of the AABB bounding box (domain)
 -A,--maxAABB FLOAT=10                 Max value of the AABB bounding box (domain)
 -n,--noise FLOAT=0                    Level of noise to perturb the shape
 -l,--lambda                           Use the shape to get a better approximation of the surface (optional)
 --properties TEXT=110                 the i-th property is disabled iff there is a 0 at position i
 -e,--estimators TEXT=110              the i-th estimator is disabled iff there is a 0 at position i

 @endcode

 @b Example:
     The following example will estimate the curvature of an implicit cone shape and will produce six resulting files (toto_II_gaussian.dat,toto_II_mean.dat, toto_MongeJetFitting_gaussian.dat, toto_MongeJetFitting_mean.dat, toto_True_gaussian.dat, toto_True_mean.dat)
     @code
     ./estimators/3dlocalEstimators --shape "z^2-x^2-y^2" --output result --h 0.4 --radius 1.0

     @endcode
      You can check other example of implicit shapes like the followinf:
     - whitney  : x^2-y*z^2
     - 4lines   : x*y*(y-x)*(y-z*x)
     - cone     : z^2-x^2-y^2
     - simonU   : x^2-z*y^2+x^4+y^4
     - cayley3  : 4*(x^2 + y^2 + z^2) + 16*x*y*z - 1
     - crixxi   : -0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3
**/


typedef std::pair<double,double> PrincipalCurvatures;
template < typename Shape, typename KSpace, typename ConstIterator, typename OutputIterator >
void
estimateTrueMeanCurvatureQuantity( const ConstIterator & it_begin,
                                  const ConstIterator & it_end,
                                  OutputIterator & output,
                                  const KSpace & K,
                                  const double & h,
                                  Shape * aShape )
{
  typedef typename KSpace::Space::RealPoint RealPoint;
  typedef CanonicSCellEmbedder< KSpace > Embedder;
  
  Embedder embedder( K );
  RealPoint currentRealPoint;
  
  for ( ConstIterator it = it_begin; it != it_end; ++it )
  {
    currentRealPoint = embedder( *it_begin ) * h;
    *output = aShape->meanCurvature( currentRealPoint );
    ++output;
  }
}


template < typename Shape, typename KSpace, typename ConstIterator, typename OutputIterator >
void
estimateTrueGaussianCurvatureQuantity( const ConstIterator & it_begin,
                                      const ConstIterator & it_end,
                                      OutputIterator & output,
                                      const KSpace & K,
                                      const double & h,
                                      Shape * aShape )
{
  typedef typename KSpace::Space::RealPoint RealPoint;
  typedef CanonicSCellEmbedder< KSpace > Embedder;
  
  Embedder embedder( K );
  RealPoint currentRealPoint;
  
  for ( ConstIterator it = it_begin; it != it_end; ++it )
  {
    currentRealPoint = embedder( *it_begin ) * h;
    *output = aShape->gaussianCurvature( currentRealPoint );
    ++output;
  }
}

template < typename Shape, typename KSpace, typename ConstIterator, typename OutputIterator >
void
estimateTruePrincipalCurvaturesQuantity( const ConstIterator & it_begin,
                                        const ConstIterator & it_end,
                                        OutputIterator & output,
                                        const KSpace & K,
                                        const double & h,
                                        Shape * aShape )
{
  typedef typename KSpace::Space::RealPoint RealPoint;
  typedef CanonicSCellEmbedder< KSpace > Embedder;
  
  Embedder embedder( K );
  RealPoint currentRealPoint;
  
  for ( ConstIterator it = it_begin; it != it_end; ++it )
  {
    currentRealPoint = embedder( *it_begin ) * h;
    double k1, k2;
    aShape->principalCurvatures( currentRealPoint, k1, k2 );
    PrincipalCurvatures result;
    result.first = k1;
    result.second = k2;
    *output = result;
    ++output;
  }
}

template <typename Space, typename Shape>
bool
compareShapeEstimators( const std::string & filename,
                       const Shape * aShape,
                       const typename Space::RealPoint & border_min,
                       const typename Space::RealPoint & border_max,
                       const double & h,
                       const double & radius_kernel,
                       const double & alpha,
                       const std::string & options,
                       const std::string & properties,
                       const bool & lambda_optimized,
                       double noiseLevel = 0.0 )
{
  typedef typename Space::RealPoint RealPoint;
  typedef GaussDigitizer< Z3i::Space, Shape > DigitalShape;
  typedef Z3i::KSpace KSpace;
  typedef typename KSpace::SCell SCell;
  typedef typename KSpace::Surfel Surfel;
  
  bool withNoise = ( noiseLevel <= 0.0 ) ? false : true;
  
  ASSERT (( noiseLevel < 1.0 ));
  // Digitizer
  DigitalShape dshape;
  dshape.attach( *aShape );
  dshape.init( border_min, border_max, h );
  
  KSpace K;
  if ( ! K.init( dshape.getLowerBound(), dshape.getUpperBound(), true ) )
  {
    std::cerr << "[3dLocalEstimators] error in creating KSpace." << std::endl;
    return false;
  }
  
  try
  {
    if ( withNoise )
    {
      typedef KanungoNoise< DigitalShape, Z3i::Domain > KanungoPredicate;
      typedef LightImplicitDigitalSurface< KSpace, KanungoPredicate > Boundary;
      typedef DigitalSurface< Boundary > MyDigitalSurface;
      typedef typename MyDigitalSurface::ConstIterator ConstIterator;
      
      typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
      typedef GraphVisitorRange< Visitor > VisitorRange;
      typedef typename VisitorRange::ConstIterator VisitorConstIterator;
      
      // typedef PointFunctorFromPointPredicateAndDomain< KanungoPredicate, Z3i::Domain, unsigned int > MyPointFunctor;
      // typedef FunctorOnCells< MyPointFunctor, KSpace > MySpelFunctor;
      
      // Extracts shape boundary
      auto d = dshape.getDomain();
      KanungoPredicate  noisifiedObject ( dshape, d, noiseLevel );
      SCell bel = Surfaces< KSpace >::findABel( K, noisifiedObject, 10000 );
      Boundary * boundary = new Boundary( K, noisifiedObject, SurfelAdjacency< KSpace::dimension >( true ), bel );
      MyDigitalSurface surf ( *boundary );
      
      double minsize = dshape.getUpperBound()[0] - dshape.getLowerBound()[0];
      unsigned int tries = 0;
      while( surf.size() < 2 * minsize || tries > 150 )
      {
        delete boundary;
        bel = Surfaces< KSpace >::findABel( K, noisifiedObject, 10000 );
        boundary = new Boundary( K, noisifiedObject, SurfelAdjacency< KSpace::dimension >( true ), bel );
        surf = MyDigitalSurface( *boundary );
        ++tries;
      }
      
      if( tries > 150 )
      {
        std::cerr << "Can't found a proper bel. So .... I ... just ... kill myself." << std::endl;
        return false;
      }
      
      VisitorRange * range;
      VisitorConstIterator ibegin;
      VisitorConstIterator iend;
      
      // Estimations
      Clock c;
      
      // True
      if( options.at( 0 ) != '0' )
      {
        // True Mean Curvature
        if( properties.at( 0 ) != '0' )
        {
          trace.beginBlock( "True mean curvature" );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_mean.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# True Mean Curvature estimation" << std::endl;
          
          std::ostream_iterator< double > out_it_true_mean( file, "\n" );
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          c.startClock();
          
          estimateTrueMeanCurvatureQuantity( ibegin,
                                            iend,
                                            out_it_true_mean,
                                            K,
                                            h,
                                            aShape );
          
          double TTrueMeanCurv = c.stopClock();
          file << "# time = " << TTrueMeanCurv << std::endl;
          
          file.close();
          delete range;
          
          trace.endBlock();
        }
        
        // True Gaussian Curvature
        if( properties.at( 1 ) != '0' )
        {
          trace.beginBlock( "True Gausian curvature" );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_gaussian.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# True Gaussian Curvature estimation" << std::endl;
          
          std::ostream_iterator< double > out_it_true_gaussian( file, "\n" );
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          c.startClock();
          
          estimateTrueGaussianCurvatureQuantity( ibegin,
                                                iend,
                                                out_it_true_gaussian,
                                                K,
                                                h,
                                                aShape );
          
          double TTrueGaussianCurv = c.stopClock();
          file << "# time = " << TTrueGaussianCurv << std::endl;
          
          file.close();
          delete range;
          
          trace.endBlock();
        }
        
        // True Principal Curvatures
        if( properties.at( 2 ) != '0' )
        {
          trace.beginBlock( "True principal curvatures" );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_principal.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# True Gaussian Curvature estimation" << std::endl;
          
          std::ostream_iterator< std::string > out_it_true_pc( file, "\n" );
          
          std::vector<PrincipalCurvatures> v_results;
          std::back_insert_iterator< std::vector<PrincipalCurvatures> > bkIt(v_results);
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          c.startClock();
          
          estimateTruePrincipalCurvaturesQuantity( ibegin,
                                                  iend,
                                                  bkIt,//out_it_true_pc,
                                                  K,
                                                  h,
                                                  aShape );
          
          for(unsigned int ii = 0; ii < v_results.size(); ++ii )
          {
            std::stringstream ss;
            ss << v_results[ii].first << " " << v_results[ii].second;
            *out_it_true_pc = ss.str();
            ++out_it_true_pc;
          }
          double TTruePrincCurv = c.stopClock();
          file << "# time = " << TTruePrincCurv << std::endl;
          
          file.close();
          
          delete range;
          
          trace.endBlock();
        }
      }
      
      double re = radius_kernel * std::pow( h, alpha ); // to obtains convergence results, re must follow the rule re=kh^(1/3)
      
      // II
      if( options.at( 1 ) != '0' )
      {
        // Integral Invariant Mean Curvature
        if( properties.at( 0 ) != '0' )
        {
          trace.beginBlock( "II mean curvature" );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_II_mean.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Mean Curvature estimation from the Integral Invariant" << std::endl;
          file << "# computed kernel radius = " << re << std::endl;
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          c.startClock();
          
          typedef functors::IIMeanCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
          typedef IntegralInvariantVolumeEstimator< KSpace, KanungoPredicate, MyIICurvatureFunctor > MyIICurvatureEstimator;
          
          MyIICurvatureFunctor curvatureFunctor;
          curvatureFunctor.init( h, re );
          
          MyIICurvatureEstimator curvatureEstimator( curvatureFunctor );
          curvatureEstimator.attach( K, noisifiedObject );
          curvatureEstimator.setParams( re/h );
          curvatureEstimator.init( h, ibegin, iend );
          
          std::ostream_iterator< double > out_it_ii_mean( file, "\n" );
          curvatureEstimator.eval( ibegin, iend, out_it_ii_mean );
          
          double TIIMeanCurv = c.stopClock();
          file << "# time = " << TIIMeanCurv << std::endl;
          file.close();
          
          delete range;
          
          trace.endBlock();
        }
        
        // Integral Invariant Gaussian Curvature
        if( properties.at( 1 ) != '0' )
        {
          trace.beginBlock( "II Gaussian curvature" );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_II_gaussian.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from the Integral Invariant" << std::endl;
          file << "# computed kernel radius = " << re << std::endl;
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          c.startClock();
          
          typedef functors::IIGaussianCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
          typedef IntegralInvariantCovarianceEstimator< KSpace, KanungoPredicate, MyIICurvatureFunctor > MyIICurvatureEstimator;
          
          MyIICurvatureFunctor curvatureFunctor;
          curvatureFunctor.init( h, re );
          
          MyIICurvatureEstimator curvatureEstimator( curvatureFunctor );
          curvatureEstimator.attach( K, noisifiedObject );
          curvatureEstimator.setParams( re/h );
          curvatureEstimator.init( h, ibegin, iend );
          
          std::ostream_iterator< double > out_it_ii_gaussian( file, "\n" );
          curvatureEstimator.eval( ibegin, iend, out_it_ii_gaussian );
          
          double TIIGaussCurv = c.stopClock();
          file << "# time = " << TIIGaussCurv << std::endl;
          file.close();
          
          delete range;
          
          trace.endBlock();
        }
        
        // Integral Invariant Principal Curvatures
        if( properties.at( 2 ) != '0' )
        {
          trace.beginBlock( "II Principal curvatures" );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_II_principal.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from the Integral Invariant" << std::endl;
          file << "# computed kernel radius = " << re << std::endl;
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          c.startClock();
          
          typedef functors::IIPrincipalCurvatures3DFunctor<Z3i::Space> MyIICurvatureFunctor;
          typedef IntegralInvariantCovarianceEstimator< KSpace, KanungoPredicate, MyIICurvatureFunctor > MyIICurvatureEstimator;
          
          MyIICurvatureFunctor curvatureFunctor;
          curvatureFunctor.init( h, re );
          
          MyIICurvatureEstimator curvatureEstimator( curvatureFunctor );
          curvatureEstimator.attach( K, noisifiedObject );
          curvatureEstimator.setParams( re/h );
          curvatureEstimator.init( h, ibegin, iend );
          
          std::vector<PrincipalCurvatures> v_results;
          std::back_insert_iterator< std::vector<PrincipalCurvatures> > bkIt(v_results);
          
          curvatureEstimator.eval( ibegin, iend, bkIt );
          
          std::ostream_iterator< std::string > out_it_ii_principal( file, "\n" );
          for( unsigned int ii = 0; ii < v_results.size(); ++ii )
          {
            std::stringstream ss;
            ss << v_results[ii].first << " " << v_results[ii].second;
            *out_it_ii_principal = ss.str();
            ++out_it_ii_principal;
          }
          
          double TIIGaussCurv = c.stopClock();
          file << "# time = " << TIIGaussCurv << std::endl;
          file.close();
          
          delete range;
          
          trace.endBlock();
        }
      }
      
      // Monge
      if( options.at( 2 ) != '0' )
      {
        // Monge Mean Curvature
        if( properties.at( 0 ) != '0' )
        {
          trace.beginBlock( "Monge mean curvature" );
          typedef LpMetric<Space> L2Metric;
          typedef functors::MongeJetFittingMeanCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorMean;
          typedef functors::ConstValue< double > ConvFunctor;
          typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, L2Metric, FunctorMean, ConvFunctor> ReporterH;
          CanonicSCellEmbedder<KSpace> embedder( K );
          FunctorMean estimatorH( embedder, h );
          ConvFunctor convFunc(1.0);
          L2Metric l2(2.0);
          ReporterH reporterH;
          reporterH.attach( surf );
          reporterH.setParams( l2, estimatorH, convFunc, re/h );
          
          c.startClock();
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          reporterH.init( h , ibegin, iend );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_mean.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Mean Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
          file << "# computed kernel radius = " << re << std::endl;
          std::ostream_iterator< double > out_it_monge_mean( file, "\n" );
          
          delete range;
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          reporterH.eval(ibegin, iend, out_it_monge_mean);
          double TMongeMeanCurv = c.stopClock();
          file << "# time = " << TMongeMeanCurv << std::endl;
          file.close();
          delete range;
          
          
          trace.endBlock();
        }
        
        // Monge Gaussian Curvature
        if( properties.at( 1 ) != '0' )
        {
          trace.beginBlock( "Monge Gaussian curvature" );
          
          typedef functors::MongeJetFittingGaussianCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorGaussian;
          typedef functors::ConstValue< double > ConvFunctor;
          typedef LpMetric<Space> L2Metric;
          typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, L2Metric, FunctorGaussian, ConvFunctor> ReporterK;
          CanonicSCellEmbedder<KSpace> embedder( K );
          FunctorGaussian estimatorK( embedder, h );
          ConvFunctor convFunc(1.0);
          ReporterK reporterK;
          reporterK.attach( surf );
          L2Metric l2(2.0);
          reporterK.setParams( l2, estimatorK, convFunc, re/h );
          c.startClock();
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          reporterK.init( h , ibegin, iend );
          
          delete range;
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_gaussian.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
          file << "# computed kernel radius = " << re << std::endl;
          std::ostream_iterator< double > out_it_monge_gaussian( file, "\n" );
          reporterK.eval(ibegin, iend , out_it_monge_gaussian);
          double TMongeGaussCurv = c.stopClock();
          file << "# time = " << TMongeGaussCurv << std::endl;
          file.close();
          delete range;
          
          
          trace.endBlock();
        }
        
        // Monge Principal Curvatures
        if( properties.at( 2 ) != '0' )
        {
          trace.beginBlock( "Monge Principal Curvature" );
          typedef LpMetric<Space> L2Metric;
          typedef functors::MongeJetFittingPrincipalCurvaturesEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorPrincCurv;
          typedef functors::ConstValue< double > ConvFunctor;
          typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, L2Metric, FunctorPrincCurv, ConvFunctor> ReporterK;
          CanonicSCellEmbedder<KSpace> embedder( K );
          FunctorPrincCurv estimatorK( embedder, h );
          ConvFunctor convFunc(1.0);
          L2Metric l2(2.0);
          ReporterK reporterK;
          reporterK.attach( surf );
          reporterK.setParams( l2, estimatorK, convFunc, re/h );
          
          c.startClock();
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          reporterK.init( h , ibegin, iend  );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_principal.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
          file << "# computed kernel radius = " << re << std::endl;
          std::ostream_iterator< std::string > out_it_monge_principal( file, "\n" );
          
          std::vector<PrincipalCurvatures> v_results;
          std::back_insert_iterator<std::vector<PrincipalCurvatures> > bkIt(v_results);
          
          delete range;
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          reporterK.eval(ibegin, iend , bkIt);//out_it_monge_principal);
          
          for(unsigned int ii = 0; ii < v_results.size(); ++ii )
          {
            std::stringstream ss;
            ss << v_results[ii].first << " " << v_results[ii].second;
            *out_it_monge_principal = ss.str();
            ++out_it_monge_principal;
          }
          
          double TMongeGaussCurv = c.stopClock();
          file << "# time = " << TMongeGaussCurv << std::endl;
          file.close();
          delete range;
          
          
          trace.endBlock();
        }
      }
    }
    else // no noise
    {
      typedef LightImplicitDigitalSurface< KSpace, DigitalShape > Boundary;
      typedef DigitalSurface< Boundary > MyDigitalSurface;
      typedef typename MyDigitalSurface::ConstIterator ConstIterator;
      
      typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
      typedef GraphVisitorRange< Visitor > VisitorRange;
      typedef typename VisitorRange::ConstIterator VisitorConstIterator;
      
      // Extracts shape boundary
      SCell bel = Surfaces<KSpace>::findABel ( K, dshape, 10000 );
      Boundary boundary( K, dshape, SurfelAdjacency< KSpace::dimension >( true ), bel );
      MyDigitalSurface surf ( boundary );
      
      VisitorRange * range;
      VisitorConstIterator ibegin;
      VisitorConstIterator iend;
      
      unsigned int cntIn = 0;
      for( typename Z3i::Domain::ConstIterator it = dshape.getDomain().begin(), ite = dshape.getDomain().end(); it != ite; ++it )
      {
        if( dshape.operator ()(*it))
        {
          ++cntIn;
        }
      }
      
      std::cout << "boundary:" << surf.size() << std::endl;
      std::cout << "full:" << cntIn << std::endl;
      
      // Estimations
      Clock c;
      
      // True
      if( options.at( 0 ) != '0' )
      {
        // True Mean Curvature
        if( properties.at( 0 ) != '0' )
        {
          trace.beginBlock( "True mean curvature" );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_mean.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# True Mean Curvature estimation" << std::endl;
          
          std::ostream_iterator< double > out_it_true_mean( file, "\n" );
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          c.startClock();
          
          estimateTrueMeanCurvatureQuantity( ibegin,
                                            iend,
                                            out_it_true_mean,
                                            K,
                                            h,
                                            aShape );
          
          double TTrueMeanCurv = c.stopClock();
          file << "# time = " << TTrueMeanCurv << std::endl;
          
          file.close();
          delete range;
          
          trace.endBlock();
        }
        
        // True Gaussian Curvature
        if( properties.at( 1 ) != '0' )
        {
          trace.beginBlock( "True Gaussian curvature" );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_gaussian.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# True Gaussian Curvature estimation" << std::endl;
          
          std::ostream_iterator< double > out_it_true_gaussian( file, "\n" );
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          c.startClock();
          
          estimateTrueGaussianCurvatureQuantity( ibegin,
                                                iend,
                                                out_it_true_gaussian,
                                                K,
                                                h,
                                                aShape );
          
          double TTrueGaussianCurv = c.stopClock();
          file << "# time = " << TTrueGaussianCurv << std::endl;
          
          file.close();
          
          delete range;
          
          trace.endBlock();
        }
        
        // True Principal Curvatures
        if( properties.at( 2 ) != '0' )
        {
          trace.beginBlock( "True principal curvatures" );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_True_principal.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# True Gaussian Curvature estimation" << std::endl;
          
          std::ostream_iterator< std::string > out_it_true_pc( file, "\n" );
          
          std::vector<PrincipalCurvatures> v_results;
          std::back_insert_iterator< std::vector<PrincipalCurvatures> > bkIt(v_results);
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          c.startClock();
          
          estimateTruePrincipalCurvaturesQuantity( ibegin,
                                                  iend,
                                                  bkIt,// out_it_true_pc,
                                                  K,
                                                  h,
                                                  aShape );
          
          
          for(unsigned int ii = 0; ii < v_results.size(); ++ii )
          {
            std::stringstream ss;
            ss << v_results[ii].first << " " << v_results[ii].second;
            *out_it_true_pc = ss.str();
            ++out_it_true_pc;
          }
          
          double TTruePrincCurv = c.stopClock();
          file << "# time = " << TTruePrincCurv << std::endl;
          
          file.close();
          
          delete range;
          
          trace.endBlock();
        }
      }
      
      double re = radius_kernel * std::pow( h, alpha ); // to obtains convergence results, re must follow the rule re=kh^(1/3)
      
      // II
      if( options.at( 1 ) != '0' )
      {
        // Integral Invariant Mean Curvature
        if( properties.at( 0 ) != '0' )
        {
          trace.beginBlock( "II mean curvature" );
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_II_mean.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Mean Curvature estimation from the Integral Invariant" << std::endl;
          file << "# computed kernel radius = " << re << std::endl;
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          c.startClock();
          
          typedef functors::IIMeanCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
          typedef IntegralInvariantVolumeEstimator< KSpace, DigitalShape, MyIICurvatureFunctor > MyIICurvatureEstimator;
          
          MyIICurvatureFunctor curvatureFunctor;
          curvatureFunctor.init( h, re );
          
          MyIICurvatureEstimator curvatureEstimator( curvatureFunctor );
          curvatureEstimator.attach( K, dshape );
          curvatureEstimator.setParams( re/h );
          curvatureEstimator.init( h, ibegin, iend );
          
          std::ostream_iterator< double > out_it_ii_mean( file, "\n" );
          curvatureEstimator.eval( ibegin, iend, out_it_ii_mean );
          
          double TIIMeanCurv = c.stopClock();
          file << "# time = " << TIIMeanCurv << std::endl;
          file.close();
          
          delete range;
          
          trace.endBlock();
        }
        
        // Integral Invariant Gaussian Curvature
        if( properties.at( 1 ) != '0' )
        {
          trace.beginBlock( "II Gaussian curvature" );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_II_gaussian.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from the Integral Invariant" << std::endl;
          file << "# computed kernel radius = " << re << std::endl;
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          c.startClock();
          
          typedef functors::IIGaussianCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
          typedef IntegralInvariantCovarianceEstimator< KSpace, DigitalShape, MyIICurvatureFunctor > MyIICurvatureEstimator;
          
          MyIICurvatureFunctor curvatureFunctor;
          curvatureFunctor.init( h, re );
          
          MyIICurvatureEstimator curvatureEstimator( curvatureFunctor );
          curvatureEstimator.attach( K, dshape );
          curvatureEstimator.setParams( re/h );
          curvatureEstimator.init( h, ibegin, iend );
          
          std::ostream_iterator< double > out_it_ii_gaussian( file, "\n" );
          curvatureEstimator.eval( ibegin, iend, out_it_ii_gaussian );
          
          double TIIGaussCurv = c.stopClock();
          file << "# time = " << TIIGaussCurv << std::endl;
          file.close();
          
          delete range;
          
          trace.endBlock();
        }
        
        // Integral Invariant Principal Curvatures
        if( properties.at( 2 ) != '0' )
        {
          trace.beginBlock( "II Principal curvatures" );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_II_principal.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from the Integral Invariant" << std::endl;
          file << "# computed kernel radius = " << re << std::endl;
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          c.startClock();
          
          typedef functors::IIPrincipalCurvatures3DFunctor<Z3i::Space> MyIICurvatureFunctor;
          typedef IntegralInvariantCovarianceEstimator< KSpace, DigitalShape, MyIICurvatureFunctor > MyIICurvatureEstimator;
          
          MyIICurvatureFunctor curvatureFunctor;
          curvatureFunctor.init( h, re );
          
          MyIICurvatureEstimator curvatureEstimator( curvatureFunctor );
          curvatureEstimator.attach( K, dshape );
          curvatureEstimator.setParams( re/h );
          curvatureEstimator.init( h, ibegin, iend );
          
          std::vector<PrincipalCurvatures> v_results;
          std::back_insert_iterator< std::vector<PrincipalCurvatures> > bkIt(v_results);
          
          curvatureEstimator.eval( ibegin, iend, bkIt );
          
          std::ostream_iterator< std::string > out_it_ii_principal( file, "\n" );
          for( unsigned int ii = 0; ii < v_results.size(); ++ii )
          {
            std::stringstream ss;
            ss << v_results[ii].first << " " << v_results[ii].second;
            *out_it_ii_principal = ss.str();
            ++out_it_ii_principal;
          }
          
          double TIIGaussCurv = c.stopClock();
          file << "# time = " << TIIGaussCurv << std::endl;
          file.close();
          
          delete range;
          
          trace.endBlock();
        }
      }
      
      // Monge
      if( options.at( 2 ) != '0' )
      {
        // Monge Mean Curvature
        if( properties.at( 0 ) != '0' )
        {
          trace.beginBlock( "Monge mean curvature" );
          typedef LpMetric<Space> L2Metric;
          typedef functors::MongeJetFittingMeanCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorMean;
          typedef functors::ConstValue< double > ConvFunctor;
          typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, L2Metric, FunctorMean, ConvFunctor> ReporterH;
          CanonicSCellEmbedder<KSpace> embedder( K );
          FunctorMean estimatorH( embedder, h );
          ConvFunctor convFunc(1.0);
          L2Metric l2(2.0);
          ReporterH reporterH;
          reporterH.attach( surf );
          reporterH.setParams( l2, estimatorH, convFunc, re/h );
          c.startClock();
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          reporterH.init( h , ibegin, iend );
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_mean.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Mean Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
          file << "# computed kernel radius = " << re << std::endl;
          std::ostream_iterator< double > out_it_monge_mean( file, "\n" );
          
          delete range;
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          reporterH.eval(ibegin, iend , out_it_monge_mean);
          double TMongeMeanCurv = c.stopClock();
          file << "# time = " << TMongeMeanCurv << std::endl;
          file.close();
          delete range;
          
          
          trace.endBlock();
        }
        
        // Monge Gaussian Curvature
        if( properties.at( 1 ) != '0' )
        {
          trace.beginBlock( "Monge Gaussian curvature" );
          typedef LpMetric<Space> L2Metric;
          typedef functors::MongeJetFittingGaussianCurvatureEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorGaussian;
          typedef functors::ConstValue< double > ConvFunctor;
          typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, L2Metric, FunctorGaussian, ConvFunctor> ReporterK;
          CanonicSCellEmbedder<KSpace> embedder( K );
          FunctorGaussian estimatorK( embedder, h );
          ConvFunctor convFunc(1.0);
          ReporterK reporterK;
          reporterK.attach( surf );
          L2Metric l2(2.0);
          reporterK.setParams( l2, estimatorK, convFunc, re/h );
          
          c.startClock();
          
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          reporterK.init( h , ibegin, iend );
          
          delete range;
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_gaussian.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
          file << "# computed kernel radius = " << re << std::endl;
          std::ostream_iterator< double > out_it_monge_gaussian( file, "\n" );
          reporterK.eval(ibegin, iend , out_it_monge_gaussian);
          double TMongeGaussCurv = c.stopClock();
          file << "# time = " << TMongeGaussCurv << std::endl;
          file.close();
          delete range;
          
          
          trace.endBlock();
        }
        
        // Monge Principal Curvatures
        if( properties.at( 2 ) != '0' )
        {
          trace.beginBlock( "Monge Principal Curvature" );
          typedef LpMetric<Space> L2Metric;
          
          typedef functors::MongeJetFittingPrincipalCurvaturesEstimator<Surfel, CanonicSCellEmbedder<KSpace> > FunctorPrincCurv;
          typedef functors::ConstValue< double > ConvFunctor;
          typedef LocalEstimatorFromSurfelFunctorAdapter<typename MyDigitalSurface::DigitalSurfaceContainer, L2Metric, FunctorPrincCurv, ConvFunctor> ReporterK;
          CanonicSCellEmbedder<KSpace> embedder( K );
          FunctorPrincCurv estimatorK( embedder, h );
          ConvFunctor convFunc(1.0);
          ReporterK reporterK;
          reporterK.attach(surf);
          L2Metric l2(2.0);
          reporterK.setParams(l2, estimatorK, convFunc, re/h);
          c.startClock();
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          reporterK.init( h , ibegin, iend );
          
          delete range;
          range = new VisitorRange( new Visitor( surf, *surf.begin() ));
          ibegin = range->begin();
          iend = range->end();
          
          char full_filename[360];
          sprintf( full_filename, "%s%s", filename.c_str(), "_MongeJetFitting_principal.dat" );
          std::ofstream file( full_filename );
          file << "# h = " << h << std::endl;
          file << "# Gaussian Curvature estimation from CGAL Monge from and Jet Fitting" << std::endl;
          file << "# computed kernel radius = " << re << std::endl;
          std::ostream_iterator< std::string > out_it_monge_principal( file, "\n" );
          
          std::vector<PrincipalCurvatures> v_results;
          std::back_insert_iterator<std::vector<PrincipalCurvatures> > bkIt(v_results);
          reporterK.eval(ibegin, iend , bkIt);//out_it_monge_principal);
          
          for(unsigned int ii = 0; ii < v_results.size(); ++ii )
          {
            std::stringstream ss;
            ss << v_results[ii].first << " " << v_results[ii].second;
            *out_it_monge_principal = ss.str();
            ++out_it_monge_principal;
          }
          
          double TMongeGaussCurv = c.stopClock();
          file << "# time = " << TMongeGaussCurv << std::endl;
          file.close();
          delete range;
          
          
          trace.endBlock();
        }
        
      }
    }
  }
  catch ( InputException e )
  {
    std::cerr << "[estimatorCurvatureComparator3D]"
    << " error."
    << e.what()
    << " "
    << filename << " "
    << border_min[0] << " "
    << border_max[0] << " "
    << h << " "
    << radius_kernel << " "
    << lambda_optimized << " "
    << options << " "
    << alpha << " "
    << std::endl;
    return false;
  }
  return true;
}

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam( std::string param )
{
  trace.error() << " Parameter: " << param << " is required.";
  trace.info() << std::endl;
  exit( 1 );
}


int main( int argc, char** argv )
{
#ifndef WITH_CGAL
#error You need to have activated CGAL (WITH_CGAL) to include this file.
#endif
#ifndef WITH_EIGEN
#error You need to have activated EIGEN (WITH_EIGEN) to include this file.
#endif
  
  // parse command line ----------------------------------------------
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  
  app.description("Compare local estimators on implicit shapes using DGtal library\n Basic usage:\n \t3dlocalEstimators --shape <shape> --h <h> --radius <radius> --estimators <binaryWord> --output <output> \n \n Below are the different available families of estimators: \n \t - Integral Invariant Mean \n \t - Integral Invariant Gaussian\n \t - Monge Jet Fitting Mean\n\t - Monge Jet Fitting Gaussian\n The i-th family of estimators is enabled if the i-th character of the binary word is not 0. The default binary word is '1100'. This means that the first family of estimators, ie. Integral Invariant, is enabled, whereas the next ones are disabled. Below are the different available properties:\n \t - Mean Curvature\n \t - Gaussian Curvature\n \t - k1/k2 \n \n Example: \n ./estimators/3dlocalEstimators --shape \"z^2-x^2-y^2\" --output result --h 0.4 --radius 1.0" );
  
  
  std::string poly_str;
  std::string file_export {"result.dat"};
  std::string properties {"110"};
  std::string options {"110"};
  
  double radius;
  double h;
  double alpha {1.0/3.0};
  double border_minV {-10.0};
  double border_maxV {10.0};
  double noiseLevel {0.0};
  bool lambda_optimized {false};

  app.add_option("--shape,-s", poly_str, "Shape") ->required();
  app.add_option("--output,-o", file_export, "Output file", true) ->required();
  
  app.add_option("--radius,-r", radius,"Kernel radius for IntegralInvariant") ->required();
  app.add_option("--alpha",alpha, "Alpha parameter for Integral Invariant computation", true );
  app.add_option("--h", h, "Grid step") ->required();
  app.add_option("--minAABB,-a", border_minV,  "Min value of the AABB bounding box (domain)", true);
  app.add_option("--maxAABB,-A", border_maxV,  "Max value of the AABB bounding box (domain)", true);
  app.add_option("--noise,-n", noiseLevel, "Level of noise to perturb the shape", true );
  
  app.add_flag("--lambda,-l", lambda_optimized, "Use the shape to get a better approximation of the surface (optional)");
  app.add_option("--properties", properties, "the i-th property is disabled iff there is a 0 at position i", true);
  app.add_option("--estimators,-e", options, "the i-th estimator is disabled iff there is a 0 at position i", true);
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------



  int nb = 3; //number of available properties
  if (properties.size() < nb)
  {
    trace.error() << " At least " << nb
    << " characters are required "
    << " with option --properties.";
    trace.info() << std::endl;
    exit(1);
  }
  
  typedef Z3i::Space::RealPoint RealPoint;
  typedef Z3i::Space::RealPoint::Coordinate Ring;
  
  RealPoint border_min( border_minV, border_minV, border_minV );
  RealPoint border_max( border_maxV, border_maxV, border_maxV );
  
  /// Construction of the polynomial shape
  
  typedef MPolynomial< 3, Ring > Polynomial3;
  typedef MPolynomialReader<3, Ring> Polynomial3Reader;
  typedef ImplicitPolynomial3Shape<Z3i::Space> ImplicitShape;
  
  Polynomial3 poly;
  Polynomial3Reader reader;
  std::string::const_iterator iter = reader.read( poly, poly_str.begin(), poly_str.end() );
  if ( iter != poly_str.end() )
  {
    std::cerr << "ERROR: I read only <"
    << poly_str.substr( 0, iter - poly_str.begin() )
    << ">, and I built P=" << poly << std::endl;
    return 1;
  }
  
  ImplicitShape* shape = new ImplicitShape( poly );
  
  compareShapeEstimators< Z3i::Space, ImplicitShape > (
                                                       file_export,
                                                       shape,
                                                       border_min, border_max,
                                                       h,
                                                       radius,
                                                       alpha,
                                                       options,
                                                       properties,
                                                       lambda_optimized,
                                                       noiseLevel );
  
  delete shape;
}
