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
 * @file volCompare.cpp
 *
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2013/07/20
 *
 *
 * This file is part of the DGtal library.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/geometry/volumes/distance/DistanceTransformation.h>
#include <DGtal/math/Statistic.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;

typedef ImageContainerBySTLVector < Z3i::Domain,  int > Image3D;
typedef ImageContainerBySTLVector < Z2i::Domain,  int > Image2D;


bool 
isDiff(const Image3D &refImage, const Image3D &compImage, 
       int refMin, int refMax,  int compMin, int compMax, const Point &pt){
  bool isRefOn = (refImage(pt)<= refMax) && (refImage(pt) >= refMin);
  bool isCompOn = (compImage(pt)<= compMax) && (compImage(pt) >= compMin);
  return ((!isRefOn && isCompOn) || (isRefOn && !isCompOn));
}



void
getStatsFromDistanceMap(Statistic<double> & stats, const Image3D &refImage, const Image3D &compImage, 
			int refMin, int refMax,  int compMin, int compMax, 
			bool statOnFalsePositiveOnly=false){

  // Get the digital set from ref image by computing the surface (use -1 and +1 since the interval of append function are open)
  Z3i::DigitalSet set3dRef (refImage.domain()); 
  SetFromImage<Z3i::DigitalSet>::append<Image3D>(set3dRef, refImage, refMin-1,refMax);  
  typedef NotPointPredicate<Z3i::DigitalSet> NegPredicate;
  

  // Applying the distance transform on the digital surface of the set: 
  typedef  DistanceTransformation<Z3i::Space, NegPredicate, Z3i::L2Metric> DTL2;   
  DTL2 dtL2(&(refImage.domain()), NegPredicate(set3dRef), &Z3i::l2Metric);

  // Get the set of point of compImage: (use -1 and +1 since the interval of append function are open)
  Z3i::DigitalSet set3dComp (compImage.domain()); 
  SetFromImage<Z3i::DigitalSet>::append<Image3D>(set3dComp, compImage, compMin-1, compMax);

  unsigned int nbAdded=0;
  //Applying stats from the set to be compared (from compImage)
  for(Z3i::DigitalSet::ConstIterator it= set3dComp.begin();  it!= set3dComp.end(); ++it){
    if((!statOnFalsePositiveOnly) || (isDiff(refImage, compImage, refMin, refMax, compMin, compMax, *it))){
      DTL2::Value distance = dtL2(*it);   
      stats.addValue(distance);
      nbAdded++;
    }
  }
  
  if(nbAdded==0)
    trace.error() << "No point added to statistics, will failed..." << endl;
}






// total ref: True Positive, True Negative, False Positive, False Negative
std::vector<int> getVoxelsStats(const Image3D &refImage,  int refMin, int refMax, const Image3D &compImage, 
				int compMin, int compMax){
  int truePos = 0; 
  int trueNeg = 0;
  int falsePos = 0;
  int falseNeg = 0;
  int totPosRef=0;
  int totPosComp=0;
  int totNegRef=0;
  int totNegComp=0;

  for(Image3D::Domain::ConstIterator it = refImage.domain().begin(); it!=refImage.domain().end(); it++){
    if(refImage(*it) <= refMax && refImage(*it) >= refMin){
      totPosRef++;
      if(compImage(*it) <= compMax && compImage(*it) >= compMin){
	totPosComp++;
	truePos++;
      }else{
	totNegComp++;
	falseNeg++;
      }
    }else{
      totNegRef++;
      if(compImage(*it) <= compMax && compImage(*it) >= compMin){
	falsePos++;
	totPosComp++;
      }else{
	totNegComp++;
	trueNeg++;
      }      
    }
  }
  std::vector<int> res;
  res.push_back(truePos);
  res.push_back(trueNeg);
  res.push_back(falsePos);
  res.push_back(falseNeg);
  res.push_back(totPosRef);
  res.push_back(totPosComp);
  res.push_back(totNegRef);
  res.push_back(totNegComp);
  
  return res;
}




int main(int argc, char**argv)
{

  
  // parse command line ----------------------------------------------
  // A = ref
  // B= comp
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "volA,a", po::value<std::string>(), "Input filename of volume A." )
    ( "volB,b", po::value<std::string>(), "Input filename of volume B." )
    ( "aMin", po::value<int>()->default_value(0), "min threshold for a voxel to be considered as belonging to the object of volume A. (default 0)" )
    ( "aMax", po::value<int>()->default_value(128), "max threshold for a voxel to be considered as belonging to the object of volume A. (default 128)" )
    ( "bMin", po::value<int>()->default_value(0), "min threshold for a voxel to be considered as belonging to the object of volume B. (default 0)" )
    ( "bMax", po::value<int>()->default_value(128), "max threshold for a voxel to be considered as belonging to the object of volume B. (default 128)" )
    ("statsFromBnotInAOnly" , "apply distance map stats only for voxels of B which are not in A (else compute stats from all distances of the object B).");
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);    
  
  if ( vm.count ( "help" ) || ! vm.count("volA")||! vm.count("volB") )
    {
      trace.info() << "apply basic comparaisons (Number of voxels (B-A), (A-B), ...etc,  statistics on distances) between two volumetric images A and B (shape defined from thresholds). Usefull to determine classical statistics like false positive related stats."<<std::endl
		   << std::endl << "Basic usage: "<<std::endl
		   << "\t volCompare --volA <volAFilename> --volB <volBFilename> "<<std::endl
		   << general_opt << "\n"
		   << "Typical use :\n  volCompare -a imageA.pgm3d --aMin 128 --aMax 255 -b imageB.pgm3d --bMin 128 --bMax 255 --statsFromBnotInAOnly \n" ;

      return 0;
    }

  if(! vm.count("volA")||! vm.count("volB"))
    {
      trace.error() << " Reference and compared volume filename are needed to be defined" << endl;      
      return 0;
    }
 
  std::string volAFilename = vm["volA"].as<std::string>();
  std::string volBFilename = vm["volB"].as<std::string>();
 
  int aMin =  vm["aMin"].as<int>();
  int aMax =  vm["aMax"].as<int>();
  int bMin = vm["bMin"].as<int>();
  int bMax = vm["bMax"].as<int>();
 
  Image3D imageA = GenericReader<Image3D>::import(volAFilename);
  Image3D imageB = GenericReader<Image3D>::import(volBFilename);
 
  std::vector<int> vectStats = getVoxelsStats(imageA, aMin, aMax,  imageB,  bMin, bMax);
  std::cout << "True Positives:" << vectStats.at(0)<< std::endl;
  std::cout << "True Negatives:" << vectStats.at(1)<< std::endl;
  std::cout << "False Positives:" << vectStats.at(2)<< std::endl;
  std::cout << "False Negatives:" << vectStats.at(3)<< std::endl;
  std::cout << "Tot Positives ref=:" << vectStats.at(4)<< std::endl; 
  std::cout << "Tot Positives comp=:" << vectStats.at(5)<< std::endl; 
  std::cout << "Tot Negatives ref=:" << vectStats.at(6)<< std::endl; 
  std::cout << "Tot Negatives comp=:" << vectStats.at(7)<< std::endl; 
 


  trace.info() << "Computing Distance Map stats ...";
  Statistic<double> statDistances(true); 
  getStatsFromDistanceMap(statDistances, imageRef, imageComp, refMin, refMax, compMin, compMax, vm.count("statsFromFalsePosOnly") );
     
  trace.info() << " [done] " << std::endl;
  std::cout << "distance max= " << statDistances.max() << std::endl;
  std::cout << "distance mean= " << statDistances.mean() << std::endl; 
  std::cout << "distance variance= " << statDistances.variance() << std::endl; 
  std::cout << "distance mediane= " << statDistances.median() << std::endl; 
  return 1;
}
