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
 * @file 3dCurvatureViewerNoise.cpp
 * @ingroup visualisation
 *
 * @author Jérémy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2013/01/10
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include <cstring>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

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

// Noise
#include "DGtal/geometry/volumes/KanungoNoise.h"

// Integral Invariant includes
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantCovarianceEstimator.h"

// Drawing
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"

#ifdef WITH_VISU3D_QGLVIEWER
#include "DGtal/io/viewers/Viewer3D.h"
#endif

using namespace DGtal;
using namespace functors;

/**
 @page Doc3DCurvatureViewerNoise 3DCurvatureViewerNoise

 @brief  Same as @ref Doc3DCurvatureViewer, but allows to add some noise to objects.

  Vol file viewer, with curvature (mean or Gaussian, see parameters) information on surface.
  Blue color means lowest curvature
  Yellow color means highest curvature
  Red means the in-between

  Uses IntegralInvariantCurvatureEstimation
  @see related article:
        Coeurjolly, D.; Lachaud, J.O; Levallois, J., (2013). Integral based Curvature
        Estimators in Digital Geometry. DGCI 2013. Retrieved from
        https://liris.cnrs.fr/publis/?id=5866




 @b Usage:  3dCurvatureViewerNoise -i file.vol --radius 5 --mode mean --noise 0.5


 @b Allowed @b options @b are:

 @code
  -h [ --help ]                         display this message
  -i [ --input ] arg                    .vol file
  -r [ --radius ] arg                   Kernel radius for IntegralInvariant
  -k [ --noise ] arg (=0.5)             Level of Kanungo noise ]0;1[
  -t [ --threshold ] arg (=8)           Min size of SCell boundary of an object
  -l [ --minImageThreshold ] arg (=0)   set the minimal image threshold to
                                        define the image object (object defined
                                        by the voxel with intensity belonging
                                        to ]minImageThreshold,
                                        maxImageThreshold ] ).
  -u [ --maxImageThreshold ] arg (=255) set the minimal image threshold to
                                        define the image object (object defined
                                        by the voxel with intensity belonging
                                        to ]minImageThreshold,
                                        maxImageThreshold] ).
  -m [ --mode ] arg (=mean)             type of output : mean, gaussian, k1,
                                        k2, prindir1, prindir2 or
                                        normal(default mean)
  -o [ --exportOBJ ] arg                Export the scene to specified OBJ/MTL
                                        filename (extensions added).
  -d [ --exportDAT ] arg                Export resulting curvature (for mean,
                                        gaussian, k1 or k2 mode) in a simple
                                        data file each line representing a
                                        surfel.
  --exportOnly                          Used to only export the result without
                                        the 3d Visualisation (usefull for
                                        scripts).
  -s [ --imageScale ] arg               scaleX, scaleY, scaleZ: re sample the
                                        source image according with a grid of
                                        size 1.0/scale (usefull to compute
                                        curvature on image defined on
                                        anisotropic grid). Set by default to
                                        1.0 for the three axis.
  -n [ --normalization ]                When exporting to OBJ, performs a
                                        normalization so that the geometry fits
                                        in [-1/2,1/2]^3


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

 @code
$ 3dCurvatureViewerNoise -i $DGtal/examples/samples/Al.100.vol -r 15 -m mean -k 0.5
 @endcode

 You should obtain such a visualisation:
 @image html res3dCurvatureViewerNoise.png "resulting visualisation of mean curvature."


 @see
 @ref 3dCurvatureViewerNoise.cpp,
 @ref Doc3DCurvatureViewer
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

namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value< std::string >(), ".vol file")
    ("radius,r",  po::value< double >(), "Kernel radius for IntegralInvariant" )
    ("noise,k",  po::value< double >()->default_value(0.5), "Level of Kanungo noise ]0;1[" )
    ("threshold,t",  po::value< unsigned int >()->default_value(8), "Min size of SCell boundary of an object" )
    ("minImageThreshold,l",  po::value<  int >()->default_value(0), "set the minimal image threshold to define the image object (object defined by the voxel with intensity belonging to ]minImageThreshold, maxImageThreshold ] )." )
    ("maxImageThreshold,u",  po::value<  int >()->default_value(255), "set the minimal image threshold to define the image object (object defined by the voxel with intensity belonging to ]minImageThreshold, maxImageThreshold] )." )
    ("mode,m", po::value< std::string >()->default_value("mean"), "type of output : mean, gaussian, k1, k2, prindir1, prindir2 or normal(default mean)")
    ("exportOBJ,o", po::value< std::string >(), "Export the scene to specified OBJ/MTL filename (extensions added)." )
    ("exportDAT,d", po::value<std::string>(), "Export resulting curvature (for mean, gaussian, k1 or k2 mode) in a simple data file each line representing a surfel. ")
    ("exportOnly", "Used to only export the result without the 3d Visualisation (usefull for scripts)." )
    ("imageScale,s", po::value<std::vector<double> >()->multitoken(), "scaleX, scaleY, scaleZ: re sample the source image according with a grid of size 1.0/scale (usefull to compute curvature on image defined on anisotropic grid). Set by default to 1.0 for the three axis.  ")
    ("normalization,n", "When exporting to OBJ, performs a normalization so that the geometry fits in [-1/2,1/2]^3") ;

  bool parseOK = true;
  po::variables_map vm;
  try
    {
      po::store( po::parse_command_line( argc, argv, general_opt ), vm );
    }
  catch( const std::exception & ex )
    {
      parseOK = false;
      trace.error() << " Error checking program options: " << ex.what() << std::endl;
    }
  bool neededArgsGiven=true;

  if (parseOK && !(vm.count("input"))){
    missingParam("--input");
    neededArgsGiven=false;
  }
  if (parseOK && !(vm.count("radius"))){
    missingParam("--radius");
    neededArgsGiven=false;
  }

  bool normalization = false;
  if (parseOK && vm.count("normalization"))
    normalization = true;

  std::string mode;
  if( parseOK )
    mode =  vm["mode"].as< std::string >();
  if ( parseOK && ( mode.compare("gaussian") != 0 ) && ( mode.compare("mean") != 0 ) &&
       ( mode.compare("k1") != 0 ) && ( mode.compare("k2") != 0 ) &&
       ( mode.compare("prindir1") != 0 ) && ( mode.compare("prindir2") != 0 )&& ( mode.compare("normal") != 0 ))
    {
      parseOK = false;
      trace.error() << " The selected mode ("<<mode << ") is not defined."<<std::endl;
    }

  double noiseLevel = vm["noise"].as< double >();
  if( noiseLevel < 0.0 || noiseLevel > 1.0 )
    {
      parseOK = false;
      trace.error() << "The noise level should be in the interval: ]0, 1["<< std::endl;
    }

#ifndef WITH_VISU3D_QGLVIEWER
  bool enable_visu = false;
#else
  bool enable_visu = !vm.count("exportOnly"); ///<! Default QGLViewer viewer. Disabled if exportOnly is set.
#endif
  bool enable_obj = vm.count("exportOBJ"); ///<! Export to a .obj file.
  bool enable_dat = vm.count("exportDAT"); ///<! Export to a .dat file.

  if( !enable_visu && !enable_obj && !enable_dat )
    {
#ifndef WITH_VISU3D_QGLVIEWER
      trace.error() << "You should specify what you want to export with --export and/or --exportDat." << std::endl;
#else
      trace.error() << "You should specify what you want to export with --export and/or --exportDat, or remove --exportOnly." << std::endl;
#endif
      neededArgsGiven = false;
    }

  if(!neededArgsGiven || !parseOK || vm.count("help") || argc <= 1 )
    {
      trace.info()<< "Visualisation of 3d curvature from .vol file using curvature from Integral Invariant" <<std::endl
                  << general_opt << "\n"
                  << "Basic usage: "<<std::endl
                  << "\t3dCurvatureViewerNoise -i file.vol --radius 5 --mode mean --noise 0.5"<<std::endl
                  << std::endl
                  << "Below are the different available modes: " << std::endl
                  << "\t - \"mean\" for the mean curvature" << std::endl
                  << "\t - \"gaussian\" for the Gaussian curvature" << std::endl
                  << "\t - \"k1\" for the first principal curvature" << std::endl
                  << "\t - \"k2\" for the second principal curvature" << std::endl
                  << "\t - \"prindir1\" for the first principal curvature direction" << std::endl
                  << "\t - \"prindir2\" for the second principal curvature direction" << std::endl
                  << "\t - \"normal\" for the normal vector" << std::endl
                  << std::endl;
      return 0;
    }
  unsigned int threshold = vm["threshold"].as< unsigned int >();
  int minImageThreshold =  vm["minImageThreshold"].as<  int >();
  int maxImageThreshold =  vm["maxImageThreshold"].as<  int >();

  double h = 1.0;

  std::string export_obj_filename;
  std::string export_dat_filename;

  if( enable_obj )
    {
      export_obj_filename = vm["exportOBJ"].as< std::string >();
      if( export_obj_filename.find(".obj") == std::string::npos )
        {
          std::ostringstream oss;
          oss << export_obj_filename << ".obj" << std::endl;
          export_obj_filename = oss.str();
        }
    }


  if( enable_dat )
    {
      export_dat_filename = vm["exportDAT"].as<std::string>();
    }

  double re_convolution_kernel = vm["radius"].as< double >();


  std::vector<  double > aGridSizeReSample;
  if( vm.count( "imageScale" ))
    {
      std::vector< double> vectScale = vm["imageScale"].as<std::vector<double > >();
      if( vectScale.size() != 3 )
        {
          trace.error() << "The grid size should contains 3 elements" << std::endl;
          return 0;
        }
      else
        {
          aGridSizeReSample.push_back(1.0/vectScale.at(0));
          aGridSizeReSample.push_back(1.0/vectScale.at(1));
          aGridSizeReSample.push_back(1.0/vectScale.at(2));
        }
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
  typedef KanungoNoise< ImagePredicate, Z3i::Domain > KanungoPredicate;
  typedef BinaryPointPredicate<DomainPredicate<Image::Domain>, KanungoPredicate, AndBoolFct2  > Predicate;
  typedef Z3i::KSpace KSpace;
  typedef KSpace::SCell SCell;
  typedef KSpace::Cell Cell;

  trace.beginBlock("Loading the file");
  std::string filename = vm["input"].as< std::string >();
  Image image = GenericReader<Image>::import( filename );

  PointVector<3,int> shiftVector3D( 0 ,0, 0 );
  DGtal::functors::BasicDomainSubSampler< HyperRectDomain< SpaceND< 3, int > >,
                                          DGtal::int32_t, double > reSampler(image.domain(),
                                                                             aGridSizeReSample,  shiftVector3D);
  const functors::Identity identityFunctor{};
  SamplerImageAdapter sampledImage ( image, reSampler.getSubSampledDomain(), reSampler, identityFunctor );
  ImagePredicate predicateIMG = ImagePredicate( sampledImage,  minImageThreshold, maxImageThreshold );
  KanungoPredicate noisifiedPredicateIMG( predicateIMG, sampledImage.domain(), noiseLevel );
  DomainPredicate<Z3i::Domain> domainPredicate( sampledImage.domain() );
  AndBoolFct2 andF;
  Predicate predicate(domainPredicate, noisifiedPredicateIMG, andF  );


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

#ifdef WITH_VISU3D_QGLVIEWER
  QApplication application( argc, argv );
  typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;
#endif
  typedef Board3D<Z3i::Space, Z3i::KSpace> Board;

#ifdef WITH_VISU3D_QGLVIEWER
  Viewer viewer( K );
#endif
  Board board( K );

#ifdef WITH_VISU3D_QGLVIEWER
  if( enable_visu )
    {
      viewer.show();
    }
#endif

  unsigned int i = 0;
  unsigned int max_size = 0;
  for( unsigned int ii = 0; ii<vectConnectedSCell.size(); ++ii )
    {
      if( vectConnectedSCell[ii].size() <= threshold )
        {
          continue;
        }
      if( vectConnectedSCell[ii].size() > max_size )
        {
          max_size = vectConnectedSCell[ii].size();
          i = ii;
        }
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
      typedef GradientColorMap< Quantity > Gradient;
      Gradient cmap_grad( min, (max==min)? max+1: max );
      cmap_grad.addColor( Color( 50, 50, 255 ) );
      cmap_grad.addColor( Color( 255, 0, 0 ) );
      cmap_grad.addColor( Color( 255, 255, 10 ) );

#ifdef WITH_VISU3D_QGLVIEWER
      if( enable_visu )
        {
          viewer << SetMode3D((*abegin2).className(), "Basic" );
        }
#endif
      if( enable_obj )
        {
          board << SetMode3D((K.unsigns(*abegin2)).className(), "Basic" );
        }


      for ( unsigned int i = 0; i < results.size(); ++i )
        {
#ifdef WITH_VISU3D_QGLVIEWER
          if( enable_visu )
            {
              viewer << CustomColors3D( Color::Black, cmap_grad( results[ i ] ));
              viewer << *abegin2;
            }
#endif

          if( enable_obj )
            {
              board << CustomColors3D( Color::Black, cmap_grad( results[ i ] ));
              board      << K.unsigns(*abegin2);
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
        }else if( mode.compare("normal") == 0 )
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

      //Visualization + export

#ifdef WITH_VISU3D_QGLVIEWER
      if( enable_visu )
        {
          viewer << SetMode3D(K.uCell( K.sKCoords(*abegin2) ).className(), "Basic" );
        }
#endif

      if( enable_obj )
        {
          board << SetMode3D(K.uCell( K.sKCoords(*abegin2) ).className(), "Basic" );
        }

      for ( unsigned int i = 0; i < results.size(); ++i )
        {
          DGtal::Dimension kDim = K.sOrthDir( *abegin2 );
          SCell outer = K.sIndirectIncident( *abegin2, kDim);
          if ( predicate(embedder(outer)) )
            {
              outer = K.sDirectIncident( *abegin2, kDim);
            }

          Cell unsignedSurfel = K.uCell( K.sKCoords(*abegin2) );

#ifdef WITH_VISU3D_QGLVIEWER
          if( enable_visu )
            {
              viewer << CustomColors3D( DGtal::Color(255,255,255,255),
                                        DGtal::Color(255,255,255,255))
                     << unsignedSurfel;
            }
#endif

          if( enable_obj )
            {
              board << CustomColors3D( DGtal::Color(255,255,255,255),
                                       DGtal::Color(255,255,255,255))
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

#ifdef WITH_VISU3D_QGLVIEWER
          if( enable_visu )
            {
              if( mode.compare("prindir1") == 0 )
                {
                  viewer.setLineColor( AXIS_COLOR_BLUE );
                }
              else if( mode.compare("prindir2") == 0 )
                {
                  viewer.setLineColor( AXIS_COLOR_RED );
                }
              else if( mode.compare("normal") == 0 )
                {
                  viewer.setLineColor( AXIS_COLOR_GREEN );
                }

              viewer.addLine (
                              RealPoint(
                                        center[0] -  0.5 * results[i][0],
                                        center[1] -  0.5 * results[i][1],
                                        center[2] -  0.5 * results[i][2]
                                        ),
                              RealPoint(
                                        center[0] +  0.5 * results[i][0],
                                        center[1] +  0.5 * results[i][1],
                                        center[2] +  0.5 * results[i][2]
                                        ),
                              AXIS_LINESIZE );
            }
#endif

          if( enable_obj )
            {
              if( mode.compare("prindir1") == 0 )
                {
                  board.setFillColor( AXIS_COLOR_BLUE );
                }
              else if( mode.compare("prindir2") == 0 )
                {
                  board.setFillColor( AXIS_COLOR_RED );
                }
              else if( mode.compare("normal") == 0 )
                {
                  board.setFillColor( AXIS_COLOR_GREEN );
                }

              board.addCylinder (
                                 RealPoint(
                                           center[0] -  0.5 * results[i][0],
                                           center[1] -  0.5 * results[i][1],
                                           center[2] -  0.5 * results[i][2]),
                                 RealPoint(
                                           center[0] +  0.5 * results[i][0],
                                           center[1] +  0.5 * results[i][1],
                                           center[2] +  0.5 * results[i][2]),
                                 0.2 );
            }

          ++abegin2;
        }
    }
  trace.endBlock();

#ifdef WITH_VISU3D_QGLVIEWER
  if( enable_visu )
    {
      viewer << Viewer3D<>::updateDisplay;
    }
#endif
  if( enable_obj )
    {
      trace.info()<< "Exporting object: " << export_obj_filename << " ...";
      board.saveOBJ(export_obj_filename,normalization);
      trace.info() << "[done]" << std::endl;
    }
  if( enable_dat )
    {
      outDat.close();
    }

#ifdef WITH_VISU3D_QGLVIEWER
  if( enable_visu )
    {
      return application.exec();
    }
#endif

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
