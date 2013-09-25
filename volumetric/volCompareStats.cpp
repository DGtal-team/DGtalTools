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
#include <limits>


using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;

typedef ImageContainerBySTLVector < Z3i::Domain,  int > Image3D;
typedef ImageContainerBySTLVector < Z2i::Domain,  int > Image2D;


double
getRMSE(const Image3D & imageA, const Image3D &imageB){
  double sumDiff=0;
  for(Image3D::Domain::ConstIterator it = imageA.domain().begin(); it!=imageA.domain().end(); it++){
    sumDiff+=(imageA(*it)-imageB(*it))*(imageA(*it)-imageB(*it));
  }
  return sqrt(sumDiff/imageA.domain().size());
}


double
getPSNR(const Image3D & imageA, const Image3D &imageB, double rmsd){
  unsigned long long int d =  std::numeric_limits<Image3D::Value>::max();
  return 10.0*log10(d*d/rmsd);
}



void
getVoxelsStats(const Image3D &imageA,  int aMin, int aMax, const Image3D &imageB, 
	       int bMin, int bMax, bool exportStatVoxels,  std::vector<Point> &vectPtBinA,  
	       std::vector<Point> &vectPtCompBinCompA,  std::vector<Point> &vectPtBnotInA, 
	       std::vector<Point> &vectPtnotBInA, bool precisionRecallFMean ){
  int numBinA = 0; // true positif with A as ref shape.
  int numCompBinCompA = 0; // true neg with A as ref shape.
  int numBnotInA = 0; // false pos with A as ref shape
  int numNotBinA = 0; // false neg with A as ref shape
  int numTotalInA = 0; // total pos in reference shape
  int numTotalInB = 0; // total pos in compared shape
  int numTotalinCompA = 0; //total in complement A 
  int numTotalinCompB = 0; //total in complement B

  for(Image3D::Domain::ConstIterator it = imageA.domain().begin(); it!=imageA.domain().end(); it++){
    // voxels in A
    if(imageA(*it) <= aMax && imageA(*it) >= aMin){
      numTotalInA++;
      //voxels in B 
      if(imageB(*it) <= bMax && imageB(*it) >= bMin){
	numTotalInB++;
	numBinA++;
	if(exportStatVoxels) vectPtBinA.push_back(*it);
      }else{
	numTotalinCompB++;
	numNotBinA++;
	if(exportStatVoxels) vectPtnotBInA.push_back(*it);
      }
      // voxels outside A
    }else{
      numTotalinCompA++;
      // voxels in B
      if(imageB(*it) <= bMax && imageB(*it) >= bMin){
	numBnotInA++;
	numTotalInB++;
	if(exportStatVoxels) vectPtBnotInA.push_back(*it);
      }else{
	numTotalinCompB++;
	numCompBinCompA++;
	if(exportStatVoxels) vectPtCompBinCompA.push_back(*it);
      }      
    }
  }

  std::cout << numBinA << " " << numCompBinCompA << " " << numBnotInA 
	    << " " << numNotBinA  << " " << numTotalInA << " " 
	    << numTotalInB << " " << numTotalinCompA << " " << numTotalinCompB;
  if(precisionRecallFMean){
    double precision = (double)numBinA/(numBinA + numBnotInA);
    double recall = (double)numBinA/(numBinA+numNotBinA);
    double fmean = (2.0*precision*recall)/(precision+recall);
    std::cout << " " << precision<<  " " << recall << " " << fmean ; 
  }
}




// total ref: True Positive, True Negative, False Positive, False Negative
void
getVoxelsStats(const Image3D &imageA,  int aMin, int aMax, const Image3D &imageB, 
				int bMin, int bMax, bool precisionRecallFMean ){
  std::vector<Point> v1, v2, v3, v4;
  return getVoxelsStats(imageA, aMin, aMax, imageB, bMin, bMax, false, v1, v2, v3, v4, precisionRecallFMean);
}


void 
exportSetofPoints(string filename, std::vector<Point> aVectPoint){
  std::ofstream ofs;
  ofs.open(filename.c_str(), std::ofstream::out );
  ofs<< "# Set of 3d points with format: x y z" << std::endl;
  for (unsigned int i =0; i< aVectPoint.size(); i++){
    ofs << aVectPoint.at(i)[0] << " " << aVectPoint.at(i)[1] << " "<< aVectPoint.at(i)[2] << std::endl;
  }
  ofs.close();
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
    ("displayTFstats", "Change the comparison diplay by using the  true/false/positive/negative notation and considering the shape A as reference. It also display precision/recall/f-mean statistics.")
    ("exportSDP", "Export voxels belonging to each categorie (voxels of ( B in A) , (NOT in B and NOT in A),   (B and NOT in A) and (Voxels of NOT in B and in A)). ") ;
  
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
      trace.info() << "apply basic comparaisons (Number of voxels (B-A), (A-B), ...etc,  RMSE, PSNR) between two volumetric images A and B (shape defined from thresholds). Usefull to determine classical statistics like false positive related stats."<<std::endl
		   << std::endl << "Basic usage: "<<std::endl
		   << "\t volCompareStats --volA <volAFilename> --volB <volBFilename> "<<std::endl
		   << general_opt << "\n"
		   << "Typical use :\n  volCompareStats -a imageA.pgm3d --aMin 128 --aMax 255 -b imageB.pgm3d --bMin 128 --bMax 255 \n" ;

      return 0;
    }

  if(! vm.count("volA")||! vm.count("volB"))
    {
      trace.error() << " The two volume filename are needed to be defined" << endl;      
      return 0;
    }
 
  std::string volAFilename = vm["volA"].as<std::string>();
  std::string volBFilename = vm["volB"].as<std::string>();
 
  int aMin = vm["aMin"].as<int>();
  int aMax = vm["aMax"].as<int>();
  int bMin = vm["bMin"].as<int>();
  int bMax = vm["bMax"].as<int>();
 
  Image3D imageA = GenericReader<Image3D>::import(volAFilename);
  Image3D imageB = GenericReader<Image3D>::import(volBFilename);
 
  if(vm.count("displayTFstats")){
    std::cout << "# Statistics given with the reference shape A: "<< volAFilename<< " (defined with threshold min: " << aMin << " and max: " << aMax << " )"<< endl;
    std::cout << "# and with the compared shape B: "<< volBFilename << "  (defined with threshold min: " << bMin << " and max: " << bMax << " )"<< endl;
    std::cout << "# #True_Positive #TrueNegative #FalsePositive #FalseNegative  #TotalinA #TotalInB #TotalComplementOfRef #TotalComplementOfComp Precision Recall F-Mean RMSE PSNR "<< endl;    
    
    if(vm.count("exportSDP")){
      std::vector<Point> voxelsBinA, voxelsNotInBNotInA, voxelsBNotInA, voxelsNotInBInA; 
      getVoxelsStats(imageA, aMin, aMax,  imageB,  bMin, bMax, true, voxelsBinA, voxelsNotInBNotInA, voxelsBNotInA, voxelsNotInBInA, true);    
      exportSetofPoints("truePos.sdp",  voxelsBinA);
      exportSetofPoints("trueNeg.sdp",  voxelsNotInBNotInA);
      exportSetofPoints("falsePos.sdp", voxelsBNotInA);
      exportSetofPoints("falseNeg.sdp", voxelsNotInBInA);
    }else{
      getVoxelsStats(imageA, aMin, aMax,  imageB,  bMin, bMax, true);
    }
    
  }else{ 
    std::cout << "# Statistics given with the shape A: "<< volAFilename<< " (defined with threshold min: " << aMin << " and max: " << aMax << " )"<< endl;
    std::cout << "# and shape B: "<< volBFilename << "  (defined with threshold min: " << bMin << " and max: " << bMax << " )"<< endl;
    std::cout << "# #(Voxels of B in A) #(Voxels of NOT in B and NOT in A) #(Voxels of B and NOT in A)  #(Voxels of NOT in B and in A) #(Voxels in A) #(Voxels in B) #(Voxels not in A) #(Voxels not in B)  RMSE PSNR "<< endl;    

    if(vm.count("exportSDP")){
      std::vector<Point> voxelsBinA, voxelsNotInBNotInA, voxelsBNotInA, voxelsNotInBInA; 
      getVoxelsStats(imageA, aMin, aMax,  imageB,  bMin, bMax, true, voxelsBinA, voxelsNotInBNotInA, voxelsBNotInA, voxelsNotInBInA, false);    
      exportSetofPoints("inBinA.sdp",  voxelsBinA);
      exportSetofPoints("notinBnotinA.sdp",  voxelsNotInBNotInA);
      exportSetofPoints("inBnotinA.sdp", voxelsBNotInA);
      exportSetofPoints("notinBinA.sdp", voxelsNotInBInA);
    }else{
      getVoxelsStats(imageA, aMin, aMax,  imageB,  bMin, bMax, false);
    }
  }  

  double rmse= getRMSE(imageA, imageB);
  double psnr= getPSNR(imageA, imageB, rmse);
  
  std::cout << " " << rmse << " " << psnr << endl;

  return 1;
}

