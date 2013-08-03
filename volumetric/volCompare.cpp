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

typedef ImageContainerBySTLVector < Z3i::Domain, unsigned int > Image3D;
typedef ImageContainerBySTLVector < Z2i::Domain, unsigned int > Image2D;





void
getStatsFromDistanceMap(Statistic<int> & stats, const Image3D &refImage, const Image3D &compImage, 
			int refMin, int refMax,  int compMin, int compMax, 
			bool statOnFalsePositiveOnly=false){

  // Get the digital set from ref image by computing the surface
  Z3i::DigitalSet set3dRef (refImage.domain()); 
  SetFromImage<Z3i::DigitalSet>::append<Image3D>(set3dRef, refImage, refMin,refMax);  
  typedef NotPointPredicate<Z3i::DigitalSet> NegPredicate;
  

  // Applying the distance transform on the digital surface of the set: 
  typedef  DistanceTransformation<Z3i::Space, NegPredicate, Z3i::L2Metric> DTL2;   
  DTL2 dtL2(&(refImage.domain()), NegPredicate(set3dRef), &Z3i::l2Metric);

  // Get the set of point of compImage:
  Z3i::DigitalSet set3dComp (compImage.domain()); 
  SetFromImage<Z3i::DigitalSet>::append<Image3D>(set3dComp, compImage, compMin, compMax);

  int distanceMax=0;
  //Applying stats from the set to be compared (from compImage)
  for(Z3i::DigitalSet::ConstIterator it= set3dComp.begin();  it!= set3dComp.end(); ++it){
    if(!statOnFalsePositiveOnly || (refImage(*it)!=compImage(*it))){
      int distance = dtL2(*it);   
      stats.addValue(distance);
    }
  }

}






// total ref: True Positive, True Negative, False Positive, False Negative
std::vector<int> getTPTNFPFNVoxelsStats(const Image3D &refImage, const Image3D &compImage, int refMin, int refMax, 
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
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "referenceVol,r", po::value<std::string>(), "Input reference vol filename." )
    ( "refMin", po::value<int>()->default_value(0), "min threshold for a voxel to be considered as belonging to the compared object. (default 0)" )
    ( "refMax", po::value<int>()->default_value(128), "max threshold for a voxel to be considered as belonging to the compared object. (default 128)" )
    ( "compMin", po::value<int>()->default_value(0), "min threshold for a voxel to be considered as belonging to the compared object. (default 0)" )
    ( "compMax", po::value<int>()->default_value(128), "max threshold for a voxel to be considered as belonging to the compared object. (default 128)" )
    ( "comparedVol,c", po::value<string>(),"Compared vol filename." )
    ("statsFromFalsePosOnly,f" , "apply distance map stats from all false positive voxels (else compute stats from all distances of the compared object).");
  
     bool parseOK=true;
     po::variables_map vm;
     try{
       po::store(po::parse_command_line(argc, argv, general_opt), vm);  
     }catch(const std::exception& ex){
       parseOK=false;
       trace.info()<< "Error checking program options: "<< ex.what()<< endl;
     }
     po::notify(vm);    
  
     if ( vm.count ( "help" ) || ! vm.count("referenceVol")||! vm.count("comparedVol") )
       {
	 trace.info() << "apply basic comparaisons (true/false positive count, statistics on distances) between two volumetric images (shape defined from thresholds)"<<std::endl
		      << std::endl << "Basic usage: "<<std::endl
		      << "\t volCompare --referenceVol <volReferenceFilename> --comparedVol <volComparedFilename> "<<std::endl
		      << general_opt << "\n"
		      << "Typical use :\n  volCompare -r imageRef.pgm3d --refMin 128 --refMax 255 -c imageDiffExt.pgm3d --compMin 128 --compMax 255 -f \n" ;


	 return 0;
       }

     if(! vm.count("referenceVol")||! vm.count("comparedVol"))
       {
	 trace.error() << " Reference and compared volume filename are needed to be defined" << endl;      
	 return 0;
       }
 
     std::string referenceVolFilename = vm["referenceVol"].as<std::string>();
     std::string compareVolFilename = vm["comparedVol"].as<std::string>();
 
     int refMin =  vm["refMin"].as<int>();
     int refMax =  vm["refMax"].as<int>();
     int compMin = vm["compMin"].as<int>();
     int compMax = vm["compMax"].as<int>();
 
     Image3D imageRef = GenericReader<Image3D>::import(referenceVolFilename);
     Image3D imageComp = GenericReader<Image3D>::import(compareVolFilename);
 
    std:vector<int> vectStats = getTPTNFPFNVoxelsStats(imageRef, imageComp, refMin, refMax, compMin, compMax);
     std::cout << "True Positives:" << vectStats.at(0)<< std::endl;
     std::cout << "True Negatives:" << vectStats.at(1)<< std::endl;
     std::cout << "False Positives:" << vectStats.at(2)<< std::endl;
     std::cout << "False Negatives:" << vectStats.at(3)<< std::endl;
     std::cout << "Tot Positives ref=:" << vectStats.at(4)<< std::endl; 
     std::cout << "Tot Positives comp=:" << vectStats.at(5)<< std::endl; 
     std::cout << "Tot Negatives ref=:" << vectStats.at(6)<< std::endl; 
     std::cout << "Tot Negatives comp=:" << vectStats.at(7)<< std::endl; 
 


     trace.info() << "Computing Distance Map stats ...";
     Statistic<int> statDistances(true); 
     getStatsFromDistanceMap(statDistances, imageRef, imageComp, refMin, refMax, compMin, compMax, vm.count("statsFromFalsePosOnly") );
     trace.info() << " [done] " << std::endl;
     std::cout << "distance max= " << statDistances.max() << std::endl;
     std::cout << "distance mean= " << statDistances.mean() << std::endl; 
     std::cout << "distance variance= " << statDistances.variance() << std::endl; 
     std::cout << "distance mediane= " << statDistances.median() << std::endl; 
     return 1;
     }
