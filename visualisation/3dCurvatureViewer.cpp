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
 * @file 3dCurvatureViewer.cpp
 * @ingroup surfaceTools
 * @author Jérémy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2013/01/10
 *
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include <cstring>

#include "CLI11.hpp"


// Shape constructors
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include <DGtal/topology/SetOfSurfels.h>

#include "DGtal/images/ImageHelper.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"

// Integral Invariant includes
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantCovarianceEstimator.h"

// Drawing
#include "DGtal/io/viewers/PolyscopeViewer.h"
#include "DGtal/io/colormaps/GradientColorMap.h"

using namespace DGtal;
using namespace functors;

/**
 @page Doc3DCurvatureViewer 3DCurvatureViewer

 @brief  Computes and visualizes mean or gaussian curvature of binary shapes.

  Vol file viewer, with curvature (mean or Gaussian, see parameters) information on surface.
  Blue color means lowest curvature
  Yellow color means highest curvature
  Red means the in-between

  Uses IntegralInvariantCurvatureEstimation
  @see related article:
        Coeurjolly, D.; Lachaud, J.O; Levallois, J., (2013). Integral based Curvature
        Estimators in Digital Geometry. DGCI 2013. Retrieved from
        https://liris.cnrs.fr/publis/?id=5866





 @b Usage:  3dCurvatureViewer -i file.vol --radius 5 --mode mean


 @b Allowed @b options @b are:

 @code

 Positionals:
   1 TEXT:FILE REQUIRED                  vol file (.vol, .longvol .p3d, .pgm3d and if DGTAL_WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         vol file (.vol, .longvol .p3d, .pgm3d and if DGTAL_WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
   -r,--radius FLOAT REQUIRED            Kernel radius for IntegralInvariant
   -t,--threshold UINT=8                 Min size of SCell boundary of an object
   -l,--minImageThreshold INT=0          set the minimal image threshold to define the image object (object defined by the voxel with intensity belonging to ]minImageThreshold, maxImageThreshold ] ).
   -u,--maxImageThreshold INT=255        set the maximal image threshold to define the image object (object defined by the voxel with intensity belonging to ]minImageThreshold, maxImageThreshold] ).
   -m,--mode TEXT:{mean,gaussian,k1,k2,prindir1,prindir2,normal}=mean
                                         type of output : mean, gaussian, k1, k2, prindir1, prindir2 or normal(default mean)
   -o,--exportOBJ TEXT                   Export the scene to specified OBJ/MTL filename (extensions added).
   -d,--exportDAT TEXT                   Export resulting curvature (for mean, gaussian, k1 or k2 mode) in a simple data file each line representing a surfel.
   --exportOnly                          Used to only export the result without the 3d Visualisation (usefull for scripts).
   -s,--imageScale FLOAT x 3             scaleX, scaleY, scaleZ: re sample the source image according with a grid of size 1.0/scale (usefull to compute curvature on image defined on anisotropic grid). Set by default to 1.0 for the three axis.
   -n,--normalization BOOLEAN            When exporting to OBJ, performs a normalization so that the geometry fits in [-1/2,1/2]^3
   
 @endcode

 Below are the different available modes:

	 - "mean" for the mean curvature
	 - "gaussian" for the Gaussian curvature
	 - "k1" for the first principal curvature
	 - "k2" for the second principal curvature
	 - "prindir1" for the first principal curvature direction
	 - "prindir2" for the second principal curvature direction
	 - "normal" for the normal vector


 @b Example:

 Now we compare the different curvature values from the two shapes:
 @code
   3dCurvatureViewer  $DGtal/examples/samples/lobster.vol -r 10 -l 40 -u 255 -m mean
 @endcode

 You should obtain such a visualisation:
 @image html res3dCurvatureViewer.png "resulting visualisation of mean curvature."


 @see
 @ref 3dCurvatureViewer.cpp,
 @ref Doc3DCurvatureViewerNoise
 */

const Color  AXIS_COLOR_RED( 200, 20, 20, 255 );
const Color  AXIS_COLOR_GREEN( 20, 200, 20, 255 );
const Color  AXIS_COLOR_BLUE( 20, 20, 200, 255 );
const double AXIS_LINESIZE = 0.05;
///////////////////////////////////////////////////////////////////////////////

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam( std::string param )
{
  trace.error() << " Parameter: " << param << " is required.";
  trace.info() << std::endl;
}


int main( int argc, char** argv )
{
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  double re_convolution_kernel;
  double noiseLevel {0.5};
  unsigned int threshold {8};
  int minImageThreshold {0};
  int maxImageThreshold {255};
  std::string mode {"mean"};
  std::string export_dat_filename;
  bool exportOnly {false};
  std::vector< double> vectScale;
  bool normalization {false};
  
  app.description("Visualisation of 3d curvature from .vol file using curvature from Integral Invarian\nBasic usage:\n \t3dCurvatureViewerNoise file.vol --radius 5 --mode mean  \n Below are the different available modes: \n\t - \"mean\" for the mean curvature \n \t - \"mean\" for the mean curvature\n\t - \"gaussian\" for the Gaussian curvature\n\t - \"k1\" for the first principal curvature\n\t - \"k2\" for the second principal curvature\n\t - \"prindir1\" for the first principal curvature direction\n\t - \"prindir2\" for the second principal curvature direction\n\t - \"normal\" for the normal vector");

   
   
   app.add_option("-i,--input,1", inputFileName, "vol file (.vol, .longvol .p3d, .pgm3d and if DGTAL_WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
   ->required()
   ->check(CLI::ExistingFile);
   
   app.add_option("--radius,-r", re_convolution_kernel, "Kernel radius for IntegralInvariant" )
     ->required();
   app.add_option("--threshold,-t", threshold, "Min size of SCell boundary of an object", true);
   app.add_option("--minImageThreshold,-l",minImageThreshold,  "set the minimal image threshold to define the image object (object defined by the voxel with intensity belonging to ]minImageThreshold, maxImageThreshold ] ).", true);
   app.add_option("--maxImageThreshold,-u",maxImageThreshold, "set the maximal image threshold to define the image object (object defined by the voxel with intensity belonging to ]minImageThreshold, maxImageThreshold] ).", true);
   app.add_option("--mode,-m", mode, "type of output : mean, gaussian, k1, k2, prindir1, prindir2 or normal(default mean)", true)
    -> check(CLI::IsMember({"mean","gaussian", "k1", "k2", "prindir1","prindir2", "normal" }));
   app.add_option("--exportDAT,-d",export_dat_filename, "Export resulting curvature (for mean, gaussian, k1 or k2 mode) in a simple data file each line representing a surfel."  );
   app.add_flag("--exportOnly", exportOnly, "Used to only export the result without the 3d Visualisation (usefull for scripts).");
   
   app.add_option("--imageScale,-s", vectScale,  "scaleX, scaleY, scaleZ: re sample the source image according with a grid of size 1.0/scale (usefull to compute curvature on image defined on anisotropic grid). Set by default to 1.0 for the three axis.")
    ->expected(3);
  
   app.add_option("--normalization,-n",normalization, "When exporting to OBJ, performs a normalization so that the geometry fits in [-1/2,1/2]^3") ;
   
   app.get_formatter()->column_width(40);
   CLI11_PARSE(app, argc, argv);
   // END parse command line using CLI ----------------------------------------------

  bool neededArgsGiven=true;

  
    bool enable_visu = !exportOnly; ///<! Default PolyscopeViewer viewer. Disabled if exportOnly is set.
    bool enable_dat = export_dat_filename != ""; ///<! Export to a .dat file.

    if( !enable_visu && !enable_dat )
      {
        trace.error() << "You should specify what you want to export with --export and/or --exportDat, or remove --exportOnly." << std::endl;
        neededArgsGiven = false;
      }

  
  double h = 1.0;
  std::vector<  double > aGridSizeReSample;
  if( vectScale.size() == 3)
     {
       aGridSizeReSample.push_back(1.0/vectScale.at(0));
       aGridSizeReSample.push_back(1.0/vectScale.at(1));
       aGridSizeReSample.push_back(1.0/vectScale.at(2));
     }
   else
     {
       aGridSizeReSample.push_back(1.0);
       aGridSizeReSample.push_back(1.0);
       aGridSizeReSample.push_back(1.0);
     }



  // Construction of the shape from vol file
  typedef Z3i::Space::RealPoint RealPoint;
  typedef Z3i::Point Point;
  typedef ImageSelector< Z3i::Domain,  int>::Type Image;
  typedef DGtal::functors::BasicDomainSubSampler< HyperRectDomain<SpaceND<3, int> >,
                                                  DGtal::int32_t, double >   ReSampler;
  typedef DGtal::ConstImageAdapter<Image, Image::Domain, ReSampler,
                                   Image::Value,  DGtal::functors::Identity >  SamplerImageAdapter;
  typedef IntervalForegroundPredicate< SamplerImageAdapter > ImagePredicate;
  typedef BinaryPointPredicate<DomainPredicate<Image::Domain>, ImagePredicate, AndBoolFct2  > Predicate;
  typedef Z3i::KSpace KSpace;
  typedef KSpace::SCell SCell;
  typedef KSpace::Cell Cell;

  trace.beginBlock("Loading the file");
  Image image = GenericReader<Image>::import( inputFileName );

  PointVector<3,int> shiftVector3D( 0 ,0, 0 );
  DGtal::functors::BasicDomainSubSampler< HyperRectDomain< SpaceND< 3, int > >,
                                          DGtal::int32_t, double > reSampler(image.domain(),
                                                                             aGridSizeReSample,  shiftVector3D);
  const functors::Identity identityFunctor{};
  SamplerImageAdapter sampledImage ( image, reSampler.getSubSampledDomain(), reSampler, identityFunctor );
  ImagePredicate predicateIMG = ImagePredicate( sampledImage,  minImageThreshold, maxImageThreshold );
  DomainPredicate<Z3i::Domain> domainPredicate( sampledImage.domain() );
  AndBoolFct2 andF;
  Predicate predicate(domainPredicate, predicateIMG, andF  );


  Z3i::Domain domain =  sampledImage.domain();
  Z3i::KSpace K;
  bool space_ok = K.init( domain.lowerBound()-Z3i::Domain::Point::diagonal(),
                          domain.upperBound()+Z3i::Domain::Point::diagonal(), true );
  if (!space_ok)
    {
      trace.error() << "Error in the Khalimsky space construction."<<std::endl;
      return 2;
    }
  CanonicSCellEmbedder< KSpace > embedder( K );
  SurfelAdjacency< Z3i::KSpace::dimension > Sadj( true );
  trace.endBlock();
  // Viewer settings


  // Extraction of components
  typedef KSpace::SurfelSet SurfelSet;
  typedef SetOfSurfels< KSpace, SurfelSet > MySetOfSurfels;
  typedef DigitalSurface< MySetOfSurfels > MyDigitalSurface;

  trace.beginBlock("Extracting surfaces");
  std::vector< std::vector<SCell > > vectConnectedSCell;
  Surfaces<KSpace>::extractAllConnectedSCell(vectConnectedSCell,K, Sadj, predicate, false);
  std::ofstream outDat;
  if( enable_dat )
    {
      trace.info() << "Exporting curvature as dat file: "<< export_dat_filename <<std::endl;
      outDat.open( export_dat_filename.c_str() );
      outDat << "# data exported from 3dCurvatureViewer implementing the II curvature estimator (Coeurjolly, D.; Lachaud, J.O; Levallois, J., (2013). Integral based Curvature"
             << "  Estimators in Digital Geometry. DGCI 2013.) " << std::endl;
      outDat << "# format: surfel coordinates (in Khalimsky space) curvature: "<< mode <<  std::endl;
    }

  trace.info()<<"Number of components= "<<vectConnectedSCell.size()<<std::endl;
  trace.endBlock();

  if( vectConnectedSCell.size() == 0 )
    {
      trace.error()<< "No surface component exists. Please check the vol file threshold parameter.";
      trace.info()<<std::endl;
      exit(2);
    }

  typedef PolyscopeViewer<Z3i::Space, Z3i::KSpace> Viewer;
  Viewer viewer( K );
  viewer.allowReuseList = true;

  for( unsigned int i = 0; i<vectConnectedSCell.size(); ++i )
    {
      if( vectConnectedSCell[i].size() <= threshold )
        {
          continue;
        }

      MySetOfSurfels  aSet(K, Sadj);

      for( std::vector<SCell>::const_iterator it = vectConnectedSCell.at(i).begin();
           it != vectConnectedSCell.at(i).end();
           ++it )
        {
          aSet.surfelSet().insert( *it);
        }

      MyDigitalSurface digSurf( aSet );


      typedef DepthFirstVisitor<MyDigitalSurface> Visitor;
      typedef GraphVisitorRange< Visitor > VisitorRange;
      typedef VisitorRange::ConstIterator SurfelConstIterator;
      VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
      SurfelConstIterator abegin = range.begin();
      SurfelConstIterator aend = range.end();

      VisitorRange range2( new Visitor( digSurf, *digSurf.begin() ) );
      SurfelConstIterator abegin2 = range2.begin();

      trace.beginBlock("Curvature computation on a component");
      if( ( mode.compare("gaussian") == 0 ) || ( mode.compare("mean") == 0 )
          || ( mode.compare("k1") == 0 ) || ( mode.compare("k2") == 0 ))
        {
          typedef double Quantity;
          std::vector< Quantity > results;
          std::back_insert_iterator< std::vector< Quantity > > resultsIterator( results );
          if ( mode.compare("mean") == 0 )
            {
              typedef functors::IIMeanCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
              typedef IntegralInvariantVolumeEstimator<Z3i::KSpace, Predicate, MyIICurvatureFunctor> MyIIEstimator;

              MyIICurvatureFunctor functor;
              functor.init( h, re_convolution_kernel );

              MyIIEstimator estimator( functor );
              estimator.attach( K, predicate );
              estimator.setParams( re_convolution_kernel/h );
              estimator.init( h, abegin, aend );

              estimator.eval( abegin, aend, resultsIterator );
            }
          else if ( mode.compare("gaussian") == 0 )
            {
              typedef functors::IIGaussianCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
              typedef IntegralInvariantCovarianceEstimator<Z3i::KSpace, Predicate, MyIICurvatureFunctor> MyIIEstimator;

              MyIICurvatureFunctor functor;
              functor.init( h, re_convolution_kernel );

              MyIIEstimator estimator( functor ); estimator.attach( K,
                                                                    predicate ); estimator.setParams( re_convolution_kernel/h );
              estimator.init( h, abegin, aend );

              estimator.eval( abegin, aend, resultsIterator );
            }
          else if ( mode.compare("k1") == 0 )
            {
              typedef functors::IIFirstPrincipalCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
              typedef IntegralInvariantCovarianceEstimator<Z3i::KSpace, Predicate, MyIICurvatureFunctor> MyIIEstimator;

              MyIICurvatureFunctor functor;
              functor.init( h, re_convolution_kernel );

              MyIIEstimator estimator( functor );
              estimator.attach( K, predicate );
              estimator.setParams( re_convolution_kernel/h );
              estimator.init( h, abegin, aend );

              estimator.eval( abegin, aend, resultsIterator );
            }
          else if ( mode.compare("k2") == 0 )
            {
              typedef functors::IISecondPrincipalCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
              typedef IntegralInvariantCovarianceEstimator<Z3i::KSpace, Predicate, MyIICurvatureFunctor> MyIIEstimator;

              MyIICurvatureFunctor functor;
              functor.init( h, re_convolution_kernel );

              MyIIEstimator estimator( functor );
              estimator.attach( K, predicate );
              estimator.setParams( re_convolution_kernel/h );
              estimator.init( h, abegin, aend );

              estimator.eval( abegin, aend, resultsIterator );
            }
          trace.endBlock();


          // Drawing results
          trace.beginBlock("Visualisation");
          Quantity min = results[ 0 ];
          Quantity max = results[ 0 ];
          for ( unsigned int i = 1; i < results.size(); ++i )
            {
              if ( results[ i ] < min )
                {
                  min = results[ i ];
                }
              else if ( results[ i ] > max )
                {
                  max = results[ i ];
                }
            }
          trace.info() << "Max value= "<<max<<"  min value= "<<min<<std::endl;
          ASSERT( min <= max );
          for ( unsigned int i = 0; i < results.size(); ++i )
            {
              if( enable_visu )
                {
                  viewer << WithQuantity(*abegin2, "curvature", results[i]);
                }

              if( enable_dat )
                {
                  Point kCoords = K.uKCoords(K.unsigns(*abegin2));
                  outDat << kCoords[0] << " " << kCoords[1] << " " << kCoords[2] <<  " " <<  results[i] << std::endl;
                }

              ++abegin2;
            }
        }
      else
        {
          typedef Z3i::Space::RealVector Quantity;
          std::vector< Quantity > results;
          std::back_insert_iterator< std::vector< Quantity > > resultsIterator( results );

          if( mode.compare("prindir1") == 0 )
            {
              typedef functors::IIFirstPrincipalDirectionFunctor<Z3i::Space> MyIICurvatureFunctor;
              typedef IntegralInvariantCovarianceEstimator<Z3i::KSpace, Predicate, MyIICurvatureFunctor> MyIIEstimator;

              MyIICurvatureFunctor functor;
              functor.init( h, re_convolution_kernel );

              MyIIEstimator estimator( functor );
              estimator.attach( K, predicate );
              estimator.setParams( re_convolution_kernel/h );
              estimator.init( h, abegin, aend );

              estimator.eval( abegin, aend, resultsIterator );
            }
          else if( mode.compare("prindir2") == 0 )
            {
              typedef functors::IISecondPrincipalDirectionFunctor<Z3i::Space> MyIICurvatureFunctor;
              typedef IntegralInvariantCovarianceEstimator<Z3i::KSpace, Predicate, MyIICurvatureFunctor> MyIIEstimator;

              MyIICurvatureFunctor functor;
              functor.init( h, re_convolution_kernel );

              MyIIEstimator estimator( functor );
              estimator.attach( K, predicate );
              estimator.setParams( re_convolution_kernel/h );
              estimator.init( h, abegin, aend );

              estimator.eval( abegin, aend, resultsIterator );
            }
          else
            if( mode.compare("normal") == 0 )
              {
                typedef functors::IINormalDirectionFunctor<Z3i::Space> MyIICurvatureFunctor;
                typedef IntegralInvariantCovarianceEstimator<Z3i::KSpace, Predicate, MyIICurvatureFunctor> MyIIEstimator;

                MyIICurvatureFunctor functor;
                functor.init( h, re_convolution_kernel );

                MyIIEstimator estimator( functor );
                estimator.attach( K, predicate );
                estimator.setParams( re_convolution_kernel/h );
                estimator.init( h, abegin, aend );

                estimator.eval( abegin, aend, resultsIterator );
              }



          ///Visualizaton / export
          for ( unsigned int i = 0; i < results.size(); ++i )
            {
              DGtal::Dimension kDim = K.sOrthDir( *abegin2 );
              SCell outer = K.sIndirectIncident( *abegin2, kDim);
              if ( predicate(Z3i::Point(embedder(outer),functors::Round<>()) ))
                {
                  outer = K.sDirectIncident( *abegin2, kDim);
                }

              Cell unsignedSurfel = K.uCell( K.sKCoords(*abegin2) );

              if( enable_visu )
                {
                  viewer << DGtal::Color(255,255,255,255)
                         << unsignedSurfel;
                }

              if( enable_dat )
                {
                  Point kCoords = K.uKCoords(K.unsigns(*abegin2));
                  outDat << kCoords[0] << " " << kCoords[1] << " " << kCoords[2] << " "
                         << results[i][0] << " " << results[i][1] << " " << results[i][2]
                         << std::endl;
                }

              RealPoint center = embedder( outer );

              if( enable_visu )
                {
                  if( mode.compare("prindir1") == 0 )
                    {
                      viewer.drawColor( AXIS_COLOR_BLUE );
                    }
                  else if( mode.compare("prindir2") == 0 )
                    {
                      viewer.drawColor( AXIS_COLOR_RED );
                    }
                  else if( mode.compare("normal") == 0 )
                    {
                      viewer.drawColor( AXIS_COLOR_GREEN );
                    }

                  viewer.drawLine(
                                  RealPoint(
                                            center[0] -  0.5 * results[i][0],
                                            center[1] -  0.5 * results[i][1],
                                            center[2] -  0.5 * results[i][2]
                                            ),
                                  RealPoint(
                                            center[0] +  0.5 * results[i][0],
                                            center[1] +  0.5 * results[i][1],
                                            center[2] +  0.5 * results[i][2]
                                            ));
                }

              ++abegin2;
            }
        }
      trace.endBlock();
    }

  if( enable_visu )
    {
      viewer.show();
    }
  if( enable_dat )
    {
      outDat.close();
    }

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
