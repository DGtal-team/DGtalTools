#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/readers/PNMReader.h"
#include "DGtal/geometry/curves/representation/FreemanChain.h"
#include "DGtal/geometry/helpers/ContourHelper.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <vector>
#include <string>

using namespace DGtal;




///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("image,i", po::value<std::string>(), "image file name")
    ("min,m", po::value<int>(), "min image threshold value (default 128)")
    ("max,M", po::value<int>(), "max image threshold value (default 255)")
    
    ("minSize,s", po::value<int>(), "minSize of the extracted freeman chain (default 0)")
    ("contourSelect,s", po::value<vector <int> >()->multitoken(), 
     "Select contour according reference point and maximal distance:  ex. --contourSelect X Y distanceMax")
    ("thresholdRange,R", po::value<vector <int> >()->multitoken(), 
     "use a range interval as threshold : --thresholdRange min increment max : for each possible i, it define a digital sets [min+(i*increment),min+((i+1)*increment)] and extract their boundary. ");
  
  
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  po::notify(vm);    
  if(vm.count("help")||argc<=1)
    {
      trace.info()<< "Extract FreemanChains from thresholded image" <<std::endl << "Basic usage: "<<std::endl
      << "\t image2freeman [options] --image <imageName> -min 128 -max 255 > contours.fc"<<std::endl
      << general_opt << "\n";
      return 0;
    }
  
  
  double minThreshold = 128;
  double maxThreshold = 255;
  unsigned int minSize =0;
  bool select=false;
  bool thresholdRange=vm.count("thresholdRange");
  Z2i::Point selectCenter;
  unsigned int selectDistanceMax = 0; 
  
  
  //Parse options
  if (!(vm.count("image"))){
    trace.info() << "Image file name needed"<< endl;
    return 0;
  } 
  
  if(vm.count("min")){
    minThreshold= vm["min"].as<int>();
  } 
  if(vm.count("max")){
    maxThreshold= vm["max"].as<int>();
  } 
  if(vm.count("minSize")){
    minSize = vm["minSize"].as<int>();
  } 
  if(vm.count("contourSelect")){
    select=true;
    vector<int> cntConstraints= vm["contourSelect"].as<vector <int> >();
    if(cntConstraints.size()!=3){
      trace.info() << "Incomplete option \"--contourSelect\""<< endl;
      return 0;
    }
    selectCenter[0]= cntConstraints.at(0);
    selectCenter[1]= cntConstraints.at(1);
    selectDistanceMax= (unsigned int) cntConstraints.at(2);
  }
  
  int min, max, increment;
  if(! thresholdRange){
    min=(int)minThreshold;
    max= (int)maxThreshold;
    increment =  (int)(maxThreshold- minThreshold);
  }else{
    vector<int> vectRange= vm["thresholdRange"].as<vector <int> >();
    if(vectRange.size()!=3){
      trace.info() << "Incomplete option \"--thresholdRange\""<< endl;
      return 0;
    }
    min=vectRange.at(0);
    increment=vectRange.at(1);
    max = vectRange.at(2);
  }

  typedef ImageSelector < Z2i::Domain, unsigned char>::Type Image;
  typedef IntervalThresholder<Image::Value> Binarizer; 





  string imageFileName = vm["image"].as<std::string>();
  Image image = PNMReader<Image>::importPGM( imageFileName ); 
  Z2i::KSpace ks;
  if(! ks.init( image.domain().lowerBound(), 
		image.domain().upperBound(), true )){
    trace.error() << "Problem in KSpace initialisation"<< endl;
  }

  for(int i=0; minThreshold+i*increment< maxThreshold; i++){
    min = (int)(minThreshold+i*increment);
    max = (int)(minThreshold+(i+1)*increment);
    Binarizer b(min, max); 
    PointFunctorPredicate<Image,Binarizer> predicate(image, b); 
    
    trace.info() << "DGtal contour extraction from thresholds ["<<  min << "," << max << "]" ;
    
    SurfelAdjacency<2> sAdj( true );
    std::vector< std::vector< Z2i::Point >  >  vectContoursBdryPointels;
    Surfaces<Z2i::KSpace>::extractAllPointContours4C( vectContoursBdryPointels,
						      ks, predicate, sAdj );  
    trace.info() << " [done]" << endl;
    for(unsigned int k=0; k<vectContoursBdryPointels.size(); k++){
      if(vectContoursBdryPointels.at(k).size()>minSize){
	if(select){
	  Z2i::Point ptMean = ContourHelper::getMeanPoint(vectContoursBdryPointels.at(k));
	  unsigned int distance = (unsigned int)ceil(sqrt((double)(ptMean[0]-selectCenter[0])*(ptMean[0]-selectCenter[0])+
							  (ptMean[1]-selectCenter[1])*(ptMean[1]-selectCenter[1])));
	  if(distance<=selectDistanceMax){
	    FreemanChain<Z2i::Integer> fc (vectContoursBdryPointels.at(k));    
	    cout << fc.x0 << " " << fc.y0   << " " << fc.chain << endl; 
	  }
	}else{
	  FreemanChain<Z2i::Integer> fc (vectContoursBdryPointels.at(k));    
	  cout << fc.x0 << " " << fc.y0   << " " << fc.chain << endl; 
	}
      }
      
    }
  
    

  }

    
  
}

