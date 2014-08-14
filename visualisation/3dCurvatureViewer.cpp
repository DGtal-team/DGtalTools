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
 * Vol file viewer, with curvature (mean or Gaussian, see parameters) information on surface.
 * Blue color means lowest curvature
 * Yellow color means highest curvature
 * Red means the in-between
 *
 * Uses IntegralInvariantCurvatureEstimation
 * @see related article:
 *       Coeurjolly, D.; Lachaud, J.O; Levallois, J., (2013). Integral based Curvature
 *       Estimators in Digital Geometry. DGCI 2013. Retrieved from
 *       https://liris.cnrs.fr/publis/?id=5866
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

// Integral Invariant includes
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantCovarianceEstimator.h"

// Drawing
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include <QtGui/QApplication>

 using namespace DGtal;

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
  ("threshold,t",  po::value< unsigned int >()->default_value(8), "Min size of SCell boundary of an object" )
  ("minImageThreshold,l",  po::value<  int >()->default_value(0), "set the minimal image threshold to define the image object (object defined by the voxel with intensity belonging to ]minImageThreshold, maxImageThreshold ] )." )
  ("maxImageThreshold,u",  po::value<  int >()->default_value(1), "set the minimal image threshold to define the image object (object defined by the voxel with intensity belonging to ]minImageThreshold, maxImageThreshold] )." )  
("mode,m", po::value< std::string >()->default_value("mean"), "type of output : mean, gaussian, k1, k2, prindir1 or prindir2 (default mean)")
  ("export,e", po::value< std::string >(), "Export the scene to specified OBJ filename." )
  ("exportDat,E", po::value<std::string>(), "Export resulting curvature (for mean, gaussian, k1 or k2 mode) in a simple data file each line representing a surfel. ")
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
                  ( mode.compare("prindir1") != 0 ) && ( mode.compare("prindir2") != 0 ))
  {
    parseOK = false;
    trace.error() << " The selected mode ("<<mode << ") is not defined."<<std::endl;
  }


  if(!neededArgsGiven || !parseOK || vm.count("help") || argc <= 1 )
  {
    trace.info()<< "Visualisation of 3d curvature from .vol file using curvature from Integral Invariant" <<std::endl
    << general_opt << "\n"
    << "Basic usage: "<<std::endl
    << "\t3dCurvatureViewer -i file.vol --radius 5 --mode mean"<<std::endl
    << std::endl
    << "Below are the different available modes: " << std::endl
    << "\t - \"mean\" for the mean curvature" << std::endl
    << "\t - \"gaussian\" for the Gaussian curvature" << std::endl
    << "\t - \"k1\" for the first principal curvature" << std::endl
    << "\t - \"k2\" for the second principal curvature" << std::endl
    << "\t - \"prindir1\" for the first principal curvature direction" << std::endl
    << "\t - \"prindir2\" for the second principal curvature direction" << std::endl
    << std::endl;
    return 0;
  }
  unsigned int threshold = vm["threshold"].as< unsigned int >();
  int minImageThreshold =  vm["minImageThreshold"].as<  int >();
  int maxImageThreshold =  vm["maxImageThreshold"].as<  int >();

  bool exportOnly = vm.count("exportOnly");

  double h = 1.0;
 

  std::string export_path;
  bool myexport = false;
  bool myexportDat = false;
  string exportDatFilename;

  if(vm.count("export")){
    export_path = vm["export"].as< std::string >();
    if( export_path.find(".obj") == std::string::npos )
    {
      std::ostringstream oss; 
      oss << export_path << ".obj" << endl; 
      export_path = oss.str();
    } 
    myexport=true;
  }


  if(vm.count("exportDat")){
    exportDatFilename = vm["exportDat"].as<std::string>();
    myexportDat = true;
  }
 
  double re_convolution_kernel = vm["radius"].as< double >();


  std::vector<  double > aGridSizeReSample;
  if(vm.count("imageScale")){
    std::vector< double> vectScale = vm["imageScale"].as<std::vector<double > >();
    if(vectScale.size()!=3){
      trace.error() << "The grid size should contains 3 elements" << std::endl;
      return 0;
    }else{
      aGridSizeReSample.push_back(1.0/vectScale.at(0));
      aGridSizeReSample.push_back(1.0/vectScale.at(1));
      aGridSizeReSample.push_back(1.0/vectScale.at(2));
    }
  }else{
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
  typedef BinaryPointPredicate<DomainPredicate<Image::Domain>, ImagePredicate, DGtal::AndBoolFct2  > Predicate;
  typedef Z3i::KSpace KSpace;
  typedef KSpace::SCell SCell;
  typedef KSpace::Cell Cell;
  typedef KSpace::Surfel Surfel;

  std::string filename = vm["input"].as< std::string >();
  Image image = GenericReader<Image>::import( filename );
  
  PointVector<3,int> shiftVector3D(0 ,0, 0);      
  DGtal::functors::BasicDomainSubSampler< HyperRectDomain<SpaceND<3, int> >,  
                                          DGtal::int32_t, double > reSampler(image.domain(),
                                                                             aGridSizeReSample,  shiftVector3D);  
  SamplerImageAdapter sampledImage (image, reSampler.getSubSampledDomain(), reSampler, functors::Identity());
  ImagePredicate predicateIMG = ImagePredicate( sampledImage,  minImageThreshold, maxImageThreshold );
  DomainPredicate<Z3i::Domain> domainPredicate( sampledImage.domain() );
  DGtal::AndBoolFct2 andF;
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

  // Viewer settings

  QApplication application( argc, argv );
  typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;

  Viewer viewer( K );
  if(!exportOnly){
    viewer.show();
  }
  // Extraction of components
  typedef KSpace::SurfelSet SurfelSet;
  typedef SetOfSurfels< KSpace, SurfelSet > MySetOfSurfels;
  typedef DigitalSurface< MySetOfSurfels > MyDigitalSurface;
  
  
  typedef Board3D<Z3i::Space, Z3i::KSpace> Board;
  Board board( K );


  std::vector< std::vector<SCell > > vectConnectedSCell;
  Surfaces<KSpace>::extractAllConnectedSCell(vectConnectedSCell,K, Sadj, predicate, false);
  std::ofstream outDat;
  if(myexportDat){
    trace.info() << "Exporting curvature as dat file: "<< exportDatFilename <<std::endl;
    outDat.open(exportDatFilename.c_str());
    outDat << "# data exported from 3dCurvatureViewer implementing the II curvature estimator (Coeurjolly, D.; Lachaud, J.O; Levallois, J., (2013). Integral based Curvature"
           << "  Estimators in Digital Geometry. DGCI 2013.) " << std::endl;
    outDat << "# format: surfel coordinates (in Khalimsky space) curvature: "<< mode <<  std::endl;
  }
  
  for(unsigned int i = 0; i<vectConnectedSCell.size(); i++)
  {
    if( vectConnectedSCell[i].size() <= threshold )
    {
      continue;
    }
    
    MySetOfSurfels  aSet(K, Sadj);

    for(std::vector<SCell>::const_iterator it= vectConnectedSCell.at(i).begin(); it != vectConnectedSCell.at(i).end(); ++it)
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

    trace.beginBlock("curvature computation");
    if( ( mode.compare("gaussian") == 0 ) || ( mode.compare("mean") == 0 ) || ( mode.compare("k1") == 0 ) || ( mode.compare("k2") == 0 ))
    {
      typedef double Quantity;
      std::vector< Quantity > results;
      back_insert_iterator< std::vector< Quantity > > resultsIterator( results );

      if ( ( mode.compare("mean") == 0 ) )
      {
        typedef functors::IIGeometricFunctors::IIMeanCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
        typedef IntegralInvariantVolumeEstimator<Z3i::KSpace, Predicate, MyIICurvatureFunctor> MyIIEstimator;

        MyIICurvatureFunctor functor;
        functor.init( h, re_convolution_kernel );

        MyIIEstimator estimator( functor );
        estimator.attach( K, predicate );
        estimator.setParams( re_convolution_kernel/h );
        estimator.init( h, abegin, aend );

        estimator.eval( abegin, aend, resultsIterator );
      }
      else if ( ( mode.compare("gaussian") == 0 ) )
      {
        typedef functors::IIGeometricFunctors::IIGaussianCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
        typedef IntegralInvariantCovarianceEstimator<Z3i::KSpace, Predicate, MyIICurvatureFunctor> MyIIEstimator;

        MyIICurvatureFunctor functor;
        functor.init( h, re_convolution_kernel );

        MyIIEstimator estimator( functor );
        estimator.attach( K, predicate );
        estimator.setParams( re_convolution_kernel/h );
        estimator.init( h, abegin, aend );

        estimator.eval( abegin, aend, resultsIterator );
      }
      else if ( ( mode.compare("k1") == 0 ) )
      {
        typedef functors::IIGeometricFunctors::IIFirstPrincipalCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
        typedef IntegralInvariantCovarianceEstimator<Z3i::KSpace, Predicate, MyIICurvatureFunctor> MyIIEstimator;

        MyIICurvatureFunctor functor;
        functor.init( h, re_convolution_kernel );

        MyIIEstimator estimator( functor );
        estimator.attach( K, predicate );
        estimator.setParams( re_convolution_kernel/h );
        estimator.init( h, abegin, aend );

        estimator.eval( abegin, aend, resultsIterator );
      }
      else if ( ( mode.compare("k2") == 0 ) )
      {
        typedef functors::IIGeometricFunctors::IISecondPrincipalCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
        typedef IntegralInvariantCovarianceEstimator<Z3i::KSpace, Predicate, MyIICurvatureFunctor> MyIIEstimator;

        MyIICurvatureFunctor functor;
        functor.init( h, re_convolution_kernel );

        MyIIEstimator estimator( functor );
        estimator.attach( K, predicate );
        estimator.setParams( re_convolution_kernel/h );
        estimator.init( h, abegin, aend );

        estimator.eval( abegin, aend, resultsIterator );
      }

      // Drawing results
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

      ASSERT( min <= max );
      typedef GradientColorMap< Quantity > Gradient;
      Gradient cmap_grad( min, (max==min)? max+1: max );
      cmap_grad.addColor( Color( 50, 50, 255 ) );
      cmap_grad.addColor( Color( 255, 0, 0 ) );
      cmap_grad.addColor( Color( 255, 255, 10 ) );

      viewer << SetMode3D((*abegin2).className(), "Basic" );
      if( myexportDat )
      {
        board << SetMode3D((K.unsigns(*abegin2)).className(), "Basic" );
      }
      

      for ( unsigned int i = 0; i < results.size(); ++i )
      {
        viewer.setFillColor(cmap_grad( results[ i ] ));
        viewer << *abegin2;

        if (myexport)
        {
          board << CustomColors3D( Color::Black, cmap_grad( results[ i ] ))
          << K.unsigns(*abegin2);
        }
        if(myexportDat){
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
      back_insert_iterator< std::vector< Quantity > > resultsIterator( results );

      if( mode.compare("prindir1") == 0 )
      {
        typedef functors::IIGeometricFunctors::IIFirstPrincipalDirectionFunctor<Z3i::Space> MyIICurvatureFunctor;
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
        typedef functors::IIGeometricFunctors::IISecondPrincipalDirectionFunctor<Z3i::Space> MyIICurvatureFunctor;
        typedef IntegralInvariantCovarianceEstimator<Z3i::KSpace, Predicate, MyIICurvatureFunctor> MyIIEstimator;

        MyIICurvatureFunctor functor;
        functor.init( h, re_convolution_kernel );

        MyIIEstimator estimator( functor );
        estimator.attach( K, predicate );
        estimator.setParams( re_convolution_kernel/h );
        estimator.init( h, abegin, aend );

        estimator.eval( abegin, aend, resultsIterator );
      }

      viewer << SetMode3D(K.uCell( K.sKCoords(*abegin2) ).className(), "Basic" );

      if( myexport )
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
        viewer.setFillColor(DGtal::Color(255,255,255,255));
        viewer << unsignedSurfel;
        if (myexport)
        {
          board << CustomColors3D( DGtal::Color(255,255,255,255),
           DGtal::Color(255,255,255,255))
          << unsignedSurfel;
        }

        RealPoint center = embedder( outer );

        if( mode.compare("prindir1") == 0 )
        {
          viewer.setLineColor(AXIS_COLOR_BLUE);
        }
        else if( mode.compare("prindir2") == 0 )
        {
          viewer.setLineColor(AXIS_COLOR_RED); 
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
        if( myexport )
        {
          board.setFillColor(AXIS_COLOR_BLUE);
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
  }

  viewer << Viewer3D<>::updateDisplay;

  if (myexport)
  {
    trace.info()<< "Exporting object: " << export_path << " ...";
    board.saveOBJ(export_path,normalization);
    trace.info() << "[done]" << std::endl;
  }
  if(myexportDat){
      outDat.close();
  }

  if(!exportOnly){
    return application.exec();
  }else{
    return 0;
  
  }
}

///////////////////////////////////////////////////////////////////////////////
