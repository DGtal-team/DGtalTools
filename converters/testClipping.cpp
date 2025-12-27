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






int main( int argc, char** argv )
{
  
  typedef ImageContainerBySTLVector< Z2i::Domain, unsigned int> Image;
  typedef IntervalThresholder<Image::Value> Binarizer; 

  
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("output,o", po::value<std::string>(), "output image file name");
  
  
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

  
  string imageOutput = vm["output"].as<std::string>();
  

  Board2D exportVecto;
  std::vector<LibBoard::Point> clippingPath;
  clippingPath.push_back(LibBoard::Point(0,0));
  clippingPath.push_back(LibBoard::Point(0,5));
  clippingPath.push_back(LibBoard::Point(5,5));
  clippingPath.push_back(LibBoard::Point(5,0));
  clippingPath.push_back(LibBoard::Point(5,0));
  clippingPath.push_back(LibBoard::Point(0,0));


  clippingPath.push_back(LibBoard::Point(1,1));
  clippingPath.push_back(LibBoard::Point(1,4));
  clippingPath.push_back(LibBoard::Point(4,4));
  clippingPath.push_back(LibBoard::Point(4,1));
  clippingPath.push_back(LibBoard::Point(4,1));
  clippingPath.push_back(LibBoard::Point(1,1));
  
  exportVecto.setClippingPath (clippingPath);


  std::vector<LibBoard::Point> polygon;
  polygon.push_back(LibBoard::Point(4,4));
  polygon.push_back(LibBoard::Point(8,6));
  polygon.push_back(LibBoard::Point(5,8));
  polygon.push_back(LibBoard::Point(2,3));
  
  exportVecto.setPenColor(DGtal::Color(255,0,0));
  exportVecto.fillPolyline (polygon);

  exportVecto.setPenColor(DGtal::Color(0,0,255));
  exportVecto.drawPolyline (clippingPath);
  
  
  
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

