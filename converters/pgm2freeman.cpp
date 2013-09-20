#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/readers/GenericReader.h"


#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageSelector.h"

#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/geometry/helpers/ContourHelper.h"

#include "DGtal/topology/helpers/Surfaces.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <vector>
#include <string>

using namespace DGtal;




///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;


typedef ImageSelector < Z2i::Domain, unsigned char>::Type Image;



std::vector<unsigned int> getHistoFromImage(const Image &image){
  const Image::Domain &imgDom = image.domain();
  std::vector<unsigned int> vectHisto(imgDom.size());
  for(Image::Domain::ConstIterator it=imgDom.begin(); it!= imgDom.end(); ++it){
    vectHisto[image(*it)]++;
  }
  return vectHisto;
}



unsigned int 
getOtsuThreshold(const Image &image){
  std::vector<unsigned int> histo = getHistoFromImage(image);
  unsigned int imageSize = image.domain().size();
  unsigned int sumA = 0;
  unsigned int sumB = imageSize;
  unsigned int muA=0;
  unsigned int muB=0;
  unsigned int sumMuAll= 0;
  for( unsigned int t=0; t< histo.size();t++){
    sumMuAll+=histo[t]*t;
  }
  
  unsigned int thresholdRes=0;
  double valMax=0.0;
  for( unsigned int t=0; t< histo.size(); t++){
    sumA+=histo[t];
    if(sumA==0)
      continue; 
    sumB=imageSize-sumA;
    if(sumB==0){
      break;
    }
    
    muA+=histo[t]*t;
    muB=sumMuAll-muA;
    double muAr=muA/(double)sumA;
    double muBr=muB/(double)sumB;
    double sigma=  (double)sumA*(double)sumB*(muAr-muBr)*(muAr-muBr);
    if(valMax<=sigma){
      valMax=sigma;
      thresholdRes=t;
    }
  }
  return thresholdRes;
}




void saveAllContoursAsFc(std::vector< std::vector< Z2i::Point >  >  vectContoursBdryPointels, unsigned int minSize){
  for(unsigned int k=0; k<vectContoursBdryPointels.size(); k++){
    if(vectContoursBdryPointels.at(k).size()>minSize){
      	  FreemanChain<Z2i::Integer> fc (vectContoursBdryPointels.at(k));    
	  std::cout << fc.x0 << " " << fc.y0   << " " << fc.chain << std::endl; 
	  
    }
  }
}


void saveSelContoursAsFC(std::vector< std::vector< Z2i::Point >  >  vectContoursBdryPointels, 
			 unsigned int minSize, Z2i::Point refPoint, double selectDistanceMax){

  for(unsigned int k=0; k<vectContoursBdryPointels.size(); k++){
    if(vectContoursBdryPointels.at(k).size()>minSize){
      Z2i::Point ptMean = ContourHelper::getMeanPoint(vectContoursBdryPointels.at(k));
      unsigned int distance = (unsigned int)ceil(sqrt((double)(ptMean[0]-refPoint[0])*(ptMean[0]-refPoint[0])+
						      (ptMean[1]-refPoint[1])*(ptMean[1]-refPoint[1])));
      if(distance<=selectDistanceMax){
	FreemanChain<Z2i::Integer> fc (vectContoursBdryPointels.at(k));    
	std::cout << fc.x0 << " " << fc.y0   << " " << fc.chain << std::endl; 
      }      
    }    
  }
}




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
    ("contourSelect,c", po::value<std::vector <int> >()->multitoken(), 
     "Select contour according reference point and maximal distance:  ex. --contourSelect X Y distanceMax")
    ("thresholdRangeMin,r", po::value<std::vector <int> >()->multitoken(), 
     "use a range interval as threshold (from min) : --thresholdRangeMin min increment max : for each possible i, it define a digital sets [min, min+((i+1)*increment)] such that min+((i+1)*increment)< max  and extract their boundary. ")
    ("thresholdRangeMax,R", po::value<std::vector <int> >()->multitoken(), 
     "use a range interval as threshold (from max) : --thresholdRangeMax min increment max : for each possible i, it define a digital sets [ max-((i)*increment), max] such that max-((i)*increment)>min  and extract their boundary. ");
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
    parseOK=false;
  }
  po::notify(vm);    
  if(vm.count("help")||argc<=1|| !parseOK)
    {
      trace.info()<< "Extract FreemanChains from thresholded image" <<std::endl << "Basic usage: "<<std::endl
      << "\t image2freeman [options] --image <imageName> -min 128 -max 255 > contours.fc"<<std::endl
      << "Note that if you don't specify any threshold a threshold threshold max is automatically defined from the Otsu algorithm with min=0. "<<std::endl
      << general_opt << "\n";
      return 0;
    }
  
  
  double minThreshold = 128;
  double maxThreshold = 255;
  unsigned int minSize =0;
  bool select=false;
  bool thresholdRange=vm.count("thresholdRangeMin")||vm.count("thresholdRangeMax");
  Z2i::Point selectCenter;
  unsigned int selectDistanceMax = 0; 
 
  typedef IntervalThresholder<Image::Value> Binarizer; 
  std::string imageFileName = vm["image"].as<std::string>();

  Image image = GenericReader<Image>::import( imageFileName ); 
  

  //Parse options
  if (!(vm.count("image"))){
    trace.info() << "Image file name needed"<< std::endl;
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
    std::vector<int> cntConstraints= vm["contourSelect"].as<std::vector <int> >();
    if(cntConstraints.size()!=3){
      trace.info() << "Incomplete option \"--contourSelect\""<< std::endl;
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
    increment =  (int)(maxThreshold - minThreshold);
    if(!vm.count("min")&&!vm.count("max")) {
      min=0;
      trace.info() << "Min/Max threshold values not specified, set min to 0 and computing max with the otsu algorithm...";     
      max = getOtsuThreshold(image);
      trace.info() << "[done] (max= " << max << ") "<< std::endl;
    }
    
  }else{
    std::vector<int> vectRange;
    if ( vm.count("thresholdRangeMax")){
      vectRange= vm["thresholdRangeMax"].as<std::vector <int> >();
    }else{
      vectRange= vm["thresholdRangeMin"].as<std::vector <int> >();
    }
    if(vectRange.size()!=3){
      trace.info() << "Incomplete option \"--thresholdRange\""<< std::endl;
      return 0;
    }
    min=vectRange.at(0);
    increment=vectRange.at(1);
    max = vectRange.at(2);
    minThreshold=min;
    maxThreshold=max;
  }

 

 
  Z2i::KSpace ks;
  if(! ks.init( image.domain().lowerBound(), 
		image.domain().upperBound(), true )){
    trace.error() << "Problem in KSpace initialisation"<< std::endl;
  }
  
  
  if (!thresholdRange){
    Binarizer b(min, max); 
    PointFunctorPredicate<Image,Binarizer> predicate(image, b); 
    trace.info() << "DGtal contour extraction from thresholds ["<<  min << "," << max << "]" ;
    SurfelAdjacency<2> sAdj( true );
    std::vector< std::vector< Z2i::Point >  >  vectContoursBdryPointels;
    Surfaces<Z2i::KSpace>::extractAllPointContours4C( vectContoursBdryPointels,
						      ks, predicate, sAdj );  
    if(select){
      saveSelContoursAsFC(vectContoursBdryPointels,  minSize, selectCenter,  selectDistanceMax);
    }else{
      saveAllContoursAsFc(vectContoursBdryPointels,  minSize); 
    }
  }else{
    for(int i=0; minThreshold+i*increment< maxThreshold; i++){
      if(vm.count("thresholdRangeMin")){
	min = (int)(minThreshold+(i)*increment);
      }
      if(vm.count("thresholdRangeMax")){
	max = (int)(maxThreshold-(i)*increment);
      }
      Binarizer b(min, max); 
      PointFunctorPredicate<Image,Binarizer> predicate(image, b); 
      
      trace.info() << "DGtal contour extraction from thresholds ["<<  min << "," << max << "]" ;
      SurfelAdjacency<2> sAdj( true );
      std::vector< std::vector< Z2i::Point >  >  vectContoursBdryPointels;
      Surfaces<Z2i::KSpace>::extractAllPointContours4C( vectContoursBdryPointels,
							ks, predicate, sAdj );  
      if(select){
	saveSelContoursAsFC(vectContoursBdryPointels,  minSize, selectCenter,  selectDistanceMax);
      }else{
	saveAllContoursAsFc(vectContoursBdryPointels,  minSize); 
      }
      trace.info() << " [done]" << std::endl;
    }
  }

    
  
}

