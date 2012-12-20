#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/readers/PNMReader.h"

#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/Image.h"

#include "DGtal/geometry/helpers/ContourHelper.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/io/boards/Board2D.h"


#include <vector>
#include <string>

using namespace DGtal;
using namespace Z2i;



///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;



int main( int argc, char** argv )
{
  
  typedef  ImageContainerBySTLVector< Z2i::Domain, unsigned int> Image;
  typedef IntervalThresholder<Image::Value> Binarizer; 
  
  
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("image,i", po::value<std::string>(), "image file name")
    ("output,o", po::value<std::string>(), "output image file name")
    ("step,s", po::value<int>(), "step value of threshold (default value 10)")
    ("sampling,p", po::value<int>(), "sampling value (default value 10)");
  
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
  unsigned int sampling =10;

  //Parse options
  if (!(vm.count("image"))){
    trace.info() << "Image file name needed"<< endl;
    return 0;
  } 
  
  string imageFileName = vm["image"].as<std::string>();
  string imageOutput = vm["output"].as<std::string>();
  Image imageSRC = PNMReader<Image>::importPGM( imageFileName ); 
  
  Z2i::KSpace ks;
  Point ptLow=imageSRC.domain().lowerBound();
  Point ptUp=imageSRC.domain().upperBound();
  
  Image image( Z2i::Domain(ptLow- Vector().diagonal(1),ptUp+ Vector().diagonal(1)));
  
  
  for(Image::Domain::ConstIterator it = imageSRC.domain().begin(),
	itend = imageSRC.domain().end(); it != itend; ++it)
    image.setValue( *it , imageSRC(*it));
  
  if(! ks.init( image.domain().lowerBound(), 
		image.domain().upperBound(), true )){
    trace.error() << "Problem in KSpace initialisation"<< endl;
  }
  
  if(vm.count("step")){
    step= vm["step"].as<int>();
  } 
  if(vm.count("sampling")){
    sampling= vm["sampling"].as<int>();
  } 


  Board2D exportVecto;
  
  SurfelAdjacency<2> sAdj( true );

  for( unsigned int i= step; i<= 255; i+=step){
    Binarizer b(0, i); 
    PointFunctorPredicate<Image,Binarizer> predicate(image, b);
    trace.info()<< "preocessing step " << i << endl;
    std::vector< std::vector< Z2i::Point >  >  vectContours;
    Surfaces<Z2i::KSpace>::extractAllPointContours4C( vectContours,
						      ks, predicate, sAdj );  
    
    for(unsigned int k=0; k< vectContours.size(); k++){
      vector<Z2i::Point> aContour = vectContours.at(k);
      vector<LibBoard::Point> aPolygon;
      
      for(unsigned int l=0; l < aContour.size(); l++){
	if(l%sampling==0){
	  aPolygon.push_back(LibBoard::Point(aContour.at(l)[0],aContour.at(l)[1]));
	}
      }
      
      unsigned int val = i;
      exportVecto.setPenColor(DGtal::Color(val,val,val));
      exportVecto.fillPolyline(aPolygon, 256-i);
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

