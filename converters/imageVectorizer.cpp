#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/readers/PNMReader.h"

#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/Image.h"

#include "DGtal/geometry/helpers/ContourHelper.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/geometry/curves/GridCurve.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "DGtal/geometry/curves/FrechetShortcut.h"
#include "DGtal/io/boards/Board2D.h"

#include "DGtal/geometry/curves/CForwardSegmentComputer.h"
#include "DGtal/geometry/curves/GreedySegmentation.h"
#include "DGtal/geometry/curves/GeometricalDSS.h"

#include <vector>
#include <string>

using namespace DGtal;
using namespace Z2i;



///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;




std::vector<LibBoard::Point> basicSampleContour(const std::vector<Point> &aContour, unsigned int rate){
  std::vector<LibBoard::Point> result;
  for(unsigned int i=0; i< aContour.size(); i++){
    if((i%rate)==0){
      result.push_back(LibBoard::Point(aContour.at(i)[0],aContour.at(i)[1]));
    } 
  }
  return result;
}



std::vector<LibBoard::Point> frechetSampleContour(const GridCurve<Z2i::K2> &curve, double errMax){
  std::vector<LibBoard::Point> result;
  
  typedef Curve::PointsRange::ConstIterator Iterator;
  typedef FrechetShortcut<Iterator,int> SegmentComputer;
  typedef GreedySegmentation<SegmentComputer> Segmentation;
  typedef GridCurve<Z2i::K2>::PointsRange Range; //range
  Range r = curve.getPointsRange(); //range
  
  typedef GreedySegmentation<SegmentComputer> Segmentation;  
  Segmentation theSegmentation( r.begin(), r.end(), SegmentComputer(errMax) );  
  Segmentation::SegmentComputerIterator it = theSegmentation.begin();
  Segmentation::SegmentComputerIterator itEnd = theSegmentation.end();
  for ( ; it != itEnd; ++it) {
    Point p = *(*it).begin();
    result.push_back(LibBoard::Point(p[0],p[1])); 

    }
  return result;
}



std::vector<LibBoard::Point> DSSSampleContour(const GridCurve<Z2i::K2> &curve){
  std::vector<LibBoard::Point> result;  
  typedef Curve::PointsRange::ConstIterator Iterator;
  typedef ArithmeticalDSS<Iterator,int,8> SegmentComputer; 
  typedef GreedySegmentation<SegmentComputer> Segmentation;
  typedef GridCurve<Z2i::K2>::PointsRange Range; //range
  Range r = curve.getPointsRange(); //range
  SegmentComputer recognitionAlgorithm;
  Segmentation theSegmentation( r.begin(), r.end(), recognitionAlgorithm );  
  Segmentation::SegmentComputerIterator it = theSegmentation.begin();
  Segmentation::SegmentComputerIterator itEnd = theSegmentation.end();
  for ( ; it != itEnd; ++it) {
    Point p = *(*it).begin();
    result.push_back(LibBoard::Point(p[0],p[1])); 

    }
  return result;
}



  
  





int main( int argc, char** argv )
{
  
  typedef ImageContainerBySTLVector< Z2i::Domain, unsigned int> Image;
  typedef IntervalThresholder<Image::Value> Binarizer; 

  
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("pgmImage,i", po::value<std::string>(), "pgm image file name")
    ("output,o", po::value<std::string>(), "output image file name")
    ("step,s", po::value<int>(), "step value of threshold (default value 10)")
    ("samplingFrechet,f", po::value<double>(), "samplingFrechet (default value 1.0)")
    ("samplingDSS,a", "samplingDSS based sampling")
    ("samplingBasic,p", po::value<int>(), "samplingBasic value (default value 10)");
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
    parseOK=false;
  }
  po::notify(vm);    
  if(vm.count("help")||argc<=1|| !parseOK)
    {
      trace.info()<< "Vectorialise  an pgm image in eps format" <<std::endl << "Basic usage: "<<std::endl
		  << "\t pgm2eps [options] --image <imageName> -step 50 > resultVecto.eos "<<std::endl
		  << general_opt << "\n";
      return 0;
    }

  unsigned int step =10;
  unsigned int samplingFreq =10;
  double frechetErr = 1.0;

  //Parse options
  if (!(vm.count("pgmImage"))){
    trace.info() << "Image file name needed"<< endl;
    return 0;
  } 
  
  string imageFileName = vm["pgmImage"].as<std::string>();
  string imageOutput = vm["output"].as<std::string>();
  Image imageSRC = PNMReader<Image>::importPGM( imageFileName ); 
  
  Z2i::KSpace ks;
  Point ptLow=imageSRC.domain().lowerBound();
  Point ptUp=imageSRC.domain().upperBound();
  
  Image image( Z2i::Domain(ptLow- Vector().diagonal(1),ptUp+ Vector().diagonal(1)));
  
  
  for(Image::Domain::ConstIterator it = imageSRC.domain().begin(),
	itend = imageSRC.domain().end(); it != itend; ++it)
    image.setValue( *it , imageSRC(*it));
  
  trace.error() <<  "lower bound " << image.domain().lowerBound() << endl;
  trace.error() <<  "upper bound " << image.domain().upperBound() << endl;
  
  if(! ks.init( image.domain().lowerBound(), 
		image.domain().upperBound(), true )){
    trace.error() << "Problem in KSpace initialisation"<< endl;
  }
  
  if(vm.count("step")){
    step= vm["step"].as<int>();
  } 
  if(vm.count("samplingBasic")){
    samplingFreq= vm["samplingBasic"].as<int>();
  } 
  if(vm.count("samplingFrechet")){
    frechetErr= vm["samplingFrechet"].as<double>();
  } 


  Board2D exportVecto;  
  SurfelAdjacency<2> sAdj( true );

  for( unsigned int i= step; i<= 255; i+=step){
    Binarizer b(i, 255); 
    PointFunctorPredicate<Image,Binarizer> predicate(image, b);
    trace.info()<< "preocessing step " << i << endl;
    std::vector< std::vector< Point >  >  vectContours;
    Surfaces<Z2i::KSpace>::extractAllPointContours4C( vectContours,
						      ks, predicate, sAdj );  
    
    for(unsigned int k=0; k< vectContours.size(); k++){

      vector<Z2i::Point> aContour = vectContours.at(k);
      GridCurve<Z2i::K2> aCurve; //grid curve
      aCurve.initFromVector(aContour);
      vector<LibBoard::Point> aContourSampled;
      if(vm.count("samplingFrechet")){
	aContourSampled  = frechetSampleContour(aCurve, frechetErr);
      }else if(vm.count("samplingDSS")){
	aContourSampled = DSSSampleContour(aCurve);
      }else{
	aContourSampled = basicSampleContour(aContour, samplingFreq);
      }
      
      if(aContourSampled.size()>0){
	exportVecto.setPenColor(DGtal::Color(i,i,i));
	exportVecto.fillPolyline(aContourSampled);
      }
    }

      

      
      
  }
  

  
  string extension = imageOutput.substr(imageOutput.find_last_of(".") + 1);
  if(extension=="eps"){
    exportVecto.saveEPS(imageOutput.c_str());  
  }else if(extension=="fig"){
    exportVecto.saveFIG(imageOutput.c_str());  
  }else if(extension=="svg"){
    exportVecto.saveFIG(imageOutput.c_str());  
  }else if(extension=="tikz"){
    exportVecto.saveTikZ(imageOutput.c_str());  
  }

  
      
  
}

