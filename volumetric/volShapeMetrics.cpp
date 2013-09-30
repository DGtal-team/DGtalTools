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
  typedef NotPointPredicate<Z3i::DigitalSet> NegPredicate;
  

  // Applying the distance transform on the digital surface of the set: 
  typedef  DistanceTransformation<Z3i::Space, NegPredicate, Z3i::L2Metric> DTL2;   
  DTL2 dtL2(&(imageA.domain()), NegPredicate(set3dRef), &Z3i::l2Metric);

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
    ("statsFromBnotInAOnly", "apply distance map stats only for voxels of B which are not in A (else compute stats from all distances of the object B).");
  
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
      trace.info() << "apply comparaisons based on distance measures between two volumetric images A and B (shape defined from thresholds)" <<std::endl
		   << std::endl << "Basic usage: "<<std::endl
		   << "\t volCompareDistances --volA <volAFilename> --volB <volBFilename> "<<std::endl
		   << general_opt << "\n"
		   << "Typical use :\n  volCompareDistances -a imageA.pgm3d --aMin 128 --aMax 255 -b imageB.pgm3d --bMin 128 --bMax 255 --statsFromBnotInAOnly \n" ;

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

  
  trace.info() << "Computing Distance Map stats ...";
  Statistic<double> statDistances(true); 
  Point ptMax;
  getStatsFromDistanceMap(statDistances, imageA, aMin, aMax, imageB, bMin, bMax, vm.count("statsFromFalsePosOnly"), &ptMax );    
  trace.info() << " [done] " << std::endl;
  std::cout << "# Statistics from distances " << std::endl;
  std::cout << "# Max(MinDistance(shape B to shape A) Mean(MinDistance(shape B to shape A) Variance(MinDistance(shape B to shape A) Mediane(MinDistance(shape B to shape A)  Farthest point of B to A " << std::endl;
  std::cout << statDistances.max() << " " << statDistances.mean()  << " " << statDistances.variance()  << "  "<< statDistances.median() << " " << ptMax[0] << " " << ptMax[1] << " " << ptMax[2] << std::endl; 
  
  
  return 1;
}

