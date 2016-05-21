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
 * @file volShapeMetrics.cpp
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



/**
 @page volShapeMetrics volShapeMetrics
 
 @brief Applies shape measures for comparing two volumetric images A and B (shape defined from thresholds).

 Usefull to determine classical statistics like false positive related stats.
 


 @b Usage: volShapeMetrics --volA \<volAFilename\> --volB \<volBFilename\> 


 @b Allowed @b options @b are : 
 @code
  -h [ --help ]               display this message.
  -a [ --volA ] arg           Input filename of volume A (vol format, and other
                              pgm3d can also be used).
  -b [ --volB ] arg           Input filename of volume B (vol format, and other
                              pgm3d can also be used).
  --aMin arg (=0)             min threshold for a voxel to be considered as 
                              belonging to the object of volume A. (default 0)
  --aMax arg (=128)           max threshold for a voxel to be considered as 
                              belonging to the object of volume A. (default 
                              128)
  --bMin arg (=0)             min threshold for a voxel to be considered as 
                              belonging to the object of volume B. (default 0)
  --bMax arg (=128)           max threshold for a voxel to be considered as 
                              belonging to the object of volume B. (default 
                              128)
  --noDistanceComparisons     to avoid to apply distance map computation if the
                              distance comparaison are not needed.
  --distancesFromBnotInAOnly  apply distance map measures only for voxels of B 
                              which are not in A (else the measure are given 
                              from all distances of the object B).
  --displayTFstats            Change the comparison diplay by using the  
                              true/false/positive/negative notation and 
                              considering the shape A as reference. It also 
                              display precision/recall/f-mean statistics.
  --exportSDP                 Export voxels belonging to each categorie (voxels
                              of ( B in A) , (NOT in B and NOT in A),   (B and 
                              NOT in A) and (Voxels of NOT in B and in A)). 

 @endcode

 @b Example: 

To test this tool, we need to generate a volumetric file to be compared to an original one:
 @code
# generation of the file "eroded.vol" from DGtal examples:
$ $DGtal/build/examples/tutorial-examples/FMMErosion
 @endcode

Then we can apply comparisons of the two shapes:

@code
$ volShapeMetrics -a eroded.vol --aMin 1 --aMax 255 -b $DGtal/examples/samples/cat10.vol --bMin 1 --bMax 255 --displayTFstats --exportSDP 
@endcode

You should obtain different comparison measures and you can display the set of voxels associated to the false positive ( falsePos.sdp ):
@code 
$ 3dSDPViewer -i falsePos.sdp -c 250 40 40 5
@endcode

 @image html resVolShapeMetrics.png "Resulting False positive set of voxels."
 
 @see
 @ref volShapeMetrics.cpp

 */



typedef ImageContainerBySTLVector < Z3i::Domain,  int > Image3D;
typedef ImageContainerBySTLVector < Z2i::Domain,  int > Image2D;


bool 
isDiff(const Image3D &imageA, int aMin, int aMax, const Image3D &imageB, 
       int bMin, int bMax, const Point &pt){
  bool isRefOn = (imageA(pt)<= aMax) && (imageA(pt) >= aMin);
  bool isCompOn = (imageB(pt)<= bMax) && (imageB(pt) >= bMin);
  return ((!isRefOn && isCompOn) || (isRefOn && !isCompOn));
}



void
getStatsFromDistanceMap(Statistic<double> & stats, const Image3D &imageA, int aMin, int aMax,
			const Image3D & imageB,  int bMin, int bMax, 
			bool statOnFalsePositiveOnly=false, Point *ptMax=0){

  // Get the digital set from ref image by computing the surface (use -1 and +1 since the interval of append function are open)
  Z3i::DigitalSet set3dRef (imageA.domain()); 
  SetFromImage<Z3i::DigitalSet>::append<Image3D>(set3dRef, imageA, aMin-1,aMax);  
  typedef functors::NotPointPredicate<Z3i::DigitalSet> NegPredicate;
  

  // Applying the distance transform on the digital surface of the set: 
  typedef  DistanceTransformation<Z3i::Space, NegPredicate, Z3i::L2Metric> DTL2;
  const NegPredicate aPredicate( set3dRef );
  DTL2 dtL2( imageA.domain(), aPredicate, Z3i::l2Metric );

  // Get the set of point of imageB: (use -1 and +1 since the interval of append function are open)
  Z3i::DigitalSet set3dComp (imageB.domain()); 
  SetFromImage<Z3i::DigitalSet>::append<Image3D>(set3dComp, imageB, bMin-1, bMax);

  unsigned int nbAdded=0;
  double maxDist=0;
  //Applying stats from the set to be compared (from imageB)
  for(Z3i::DigitalSet::ConstIterator it= set3dComp.begin();  it!= set3dComp.end(); ++it){
    if((!statOnFalsePositiveOnly) || (isDiff(imageA, aMin, aMax, imageB, bMin, bMax, *it))){
      DTL2::Value distance = dtL2(*it);   
      stats.addValue(distance);
      nbAdded++;
      if(maxDist<distance){
	maxDist=distance;
	if(ptMax!=0){
	  (*ptMax)[0]=(*it)[0]; 
	  (*ptMax)[1]=(*it)[1];
	  (*ptMax)[2]=(*it)[2];
	}
      }
    }
  }
  
  if(nbAdded==0)
    trace.error() << "No point added to statistics, will failed..." << endl;
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
    ( "volA,a", po::value<std::string>(), "Input filename of volume A (vol format, and other pgm3d can also be used)." )
    ( "volB,b", po::value<std::string>(), "Input filename of volume B (vol format, and other pgm3d can also be used)." )
    ( "aMin", po::value<int>()->default_value(0), "min threshold for a voxel to be considered as belonging to the object of volume A. (default 0)" )
    ( "aMax", po::value<int>()->default_value(128), "max threshold for a voxel to be considered as belonging to the object of volume A. (default 128)" )
    ( "bMin", po::value<int>()->default_value(0), "min threshold for a voxel to be considered as belonging to the object of volume B. (default 0)" )
    ( "bMax", po::value<int>()->default_value(128), "max threshold for a voxel to be considered as belonging to the object of volume B. (default 128)" )
    ( "noDistanceComparisons", "to avoid to apply distance map computation if the distance comparaison are not needed.") 
    ("distancesFromBnotInAOnly", "apply distance map measures only for voxels of B which are not in A (else the measure are given from all distances of the object B).")
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
  
  if ( !parseOK || vm.count ( "help" ) || ! vm.count("volA")||! vm.count("volB") )
    {
      trace.info() << "Apply shape measures for comparing two volumetric images A and B (shape defined from thresholds). " << std::endl
                   << "It can compute:" << std::endl
                   << "      - voxel count from voxel partition (number of voxel from (B-A), (A-B) ...): usefull to determine classical statistics like false positive related stats."<<std::endl
		   << "      - euclidean distance between two volumetric images A and B " << std::endl
		   << std::endl << "Basic usage: "<< std::endl
		   << "\t volShapeMetrics --volA <volAFilename> --volB <volBFilename> "<< std::endl
		   << general_opt << "\n"
		   << "Typical use :\n  volShapeMetrics -a imageA.vol --aMin 128 --aMax 255 -b imageB.vol --bMin 128 --bMax 255 --distancesFromBnotInAOnly \n" ;
 
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

  
  std::cout << "# Shape comparisons (generated with volShapeMetrics) given with the reference shape A: "<< volAFilename
	    << " (defined with threshold min: " << aMin << " and max: " << aMax << " )"<< endl
	    << "# and with the compared shape B: "<< volBFilename << "  (defined with threshold min: " << bMin << " and max: " << bMax << " )"<< endl;
  if(vm.count("displayTFstats")){
    std::cout << "# #True_Positive #TrueNegative #FalsePositive #FalseNegative  #TotalinA #TotalInB #TotalComplementOfRef #TotalComplementOfComp Precision Recall F-Mean  ";
  }else{
    std::cout << "# #(Voxels of B in A) #(Voxels of NOT in B and NOT in A) #(Voxels of B and NOT in A)  #(Voxels of NOT in B and in A) #(Voxels in A) #(Voxels in B) #(Voxels not in A) #(Voxels not in B) ";
  }
  
  
  if(!vm.count("noDistanceComparisons")){
    std::cout << " Max(MinDistance(shape B to shape A) Mean(MinDistance(shape B to shape A) Variance(MinDistance(shape B to shape A)) "
	      << " Mediane(MinDistance(shape B to shape A)  Farthest point of B to A ";
    if(vm.count("distancesFromBnotInAOnly")){
      std::cout << "*** for parts of B which are not in A only ***" ;
    }
  }
  std::cout << std::endl; 

  
  if(vm.count("exportSDP")){
    std::vector<Point> voxelsBinA, voxelsNotInBNotInA, voxelsBNotInA, voxelsNotInBInA; 
    if(vm.count("displayTFstats")){
      getVoxelsStats(imageA, aMin, aMax,  imageB,  bMin, bMax, true, voxelsBinA, voxelsNotInBNotInA, voxelsBNotInA, voxelsNotInBInA, true);    
      exportSetofPoints("truePos.sdp",  voxelsBinA);
      exportSetofPoints("trueNeg.sdp",  voxelsNotInBNotInA);
      exportSetofPoints("falsePos.sdp", voxelsBNotInA);
      exportSetofPoints("falseNeg.sdp", voxelsNotInBInA);
    }else{
      std::vector<Point> voxelsBinA, voxelsNotInBNotInA, voxelsBNotInA, voxelsNotInBInA; 
      getVoxelsStats(imageA, aMin, aMax,  imageB,  bMin, bMax, true, voxelsBinA, voxelsNotInBNotInA, voxelsBNotInA, voxelsNotInBInA, false);    
      exportSetofPoints("inBinA.sdp",  voxelsBinA);
      exportSetofPoints("notinBnotinA.sdp",  voxelsNotInBNotInA);
      exportSetofPoints("inBnotinA.sdp", voxelsBNotInA);
      exportSetofPoints("notinBinA.sdp", voxelsNotInBInA);
    }
  }else{
    getVoxelsStats(imageA, aMin, aMax,  imageB,  bMin, bMax, vm.count("displayTFstats"));
  }
  
  
  if(!vm.count("noDistanceComparisons")){
    trace.info() << "Computing Distance Map stats ...";
    Statistic<double> statDistances(true); 
    Point ptMax;
    getStatsFromDistanceMap(statDistances, imageA, aMin, aMax, imageB, bMin, bMax, vm.count("statsFromFalsePosOnly"), &ptMax );    
    std::cout << " " << statDistances.max() << " " << statDistances.mean()  << " " << statDistances.variance()  << "  "<< statDistances.median() << " " << ptMax[0] << " " << ptMax[1] << " " << ptMax[2] << std::endl; 

    trace.info() << " [done] " << std::endl;
  }    
   
   return 1;
      
}
