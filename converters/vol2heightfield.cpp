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
 * @file vol2heightfield.cpp
 * @ingroup Converters
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2015/03/18
 *
 * 
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;


/**
 @page vol2heightfield vol2heightfield
 @brief  Converts volumetric  file into a projected 2D image given from a normal direction N and from a starting point P.

 The 3D volume is scanned in this normal direction N starting from P with a step 1. If the intensity of the 3d point is inside the given thresholds its 2D gray values are set to the current scan number.


@b Usage: vol2heightfield [input] [output]

@b Allowed @b options @b are:

@code
  -h [ --help ]                    display this message
  -i [ --input ] arg               vol file (.vol, .longvol .p3d, .pgm3d and if
                                   WITH_ITK is selected: dicom, dcm, mha, mhd).
                                   For longvol, dicom, dcm, mha or mhd formats,
                                   the input values are linearly scaled between
                                   0 and 255.
  -o [ --output ] arg              sequence of discrete point file (.sdp) 
  -m [ --thresholdMin ] arg (=128) min threshold (default 128)
  -M [ --thresholdMax ] arg (=255) max threshold (default 255)
  --nx arg (=0)                    set the x component of the projection 
                                   direction.
  --ny arg (=0)                    set the y component of the projection 
                                   direction.
  --nz arg (=1)                    set the z component of the projection 
                                   direction.
  -x [ --centerX ] arg (=0)        choose x center of the projected image.
  -y [ --centerY ] arg (=0)        choose y center of the projected image.
  -z [ --centerZ ] arg (=1)        choose z center of the projected image.
  --width arg (=100)               set the width of the resulting height Field 
                                   image.
  --height arg (=100)              set the height of the resulting height Field
                                   image.
  --heightFieldMaxScan arg (=255)  set the maximal scan deep.
  --setBackgroundLastDepth         change the default background (black with 
                                   the last filled intensity).
  --rescaleInputMin arg (=0)       min value used to rescale the input 
                                   intensity (to avoid basic cast into 8  bits 
                                   image).
  --rescaleInputMax arg (=255)     max value used to rescale the input 
                                   intensity (to avoid basic cast into 8 bits 
                                   image).


@endcode

@b Example:
@code 
$ vol2heightfield -i ${DGtal}/examples/samples/lobster.vol -m 60 -M 500  --nx 0 --ny 0.7 --nz -1 -x 150 -y 0 -z 150 --width 300 --height 300 --heightFieldMaxScan 350  -o resultingHeightMap.pgm 
@endcode

You should obtain such a resulting image:
@image html resVol2heightfield.png "resulting image."
@see
@ref vol2heightfield.cpp

*/


int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char> Image2D;
  typedef DGtal::ConstImageAdapter<Image3D, Z2i::Domain, DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>,
                                   Image3D::Value,  DGtal::functors::Identity >  ImageAdapterExtractor;

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
     ("input,i", po::value<std::string>(), "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
    ("output,o", po::value<std::string>(), "sequence of discrete point file (.sdp) ") 
    ("thresholdMin,m", po::value<int>()->default_value(128), "min threshold (default 128)" )
    ("thresholdMax,M", po::value<int>()->default_value(255), "max threshold (default 255)" )
    ("nx", po::value<double>()->default_value(0), "set the x component of the projection direction." )
    ("ny", po::value<double>()->default_value(0), "set the y component of the projection direction." )
    ("nz", po::value<double>()->default_value(1), "set the z component of the projection direction." )
    ("centerX,x", po::value<unsigned int>()->default_value(0), "choose x center of the projected image." )
    ("centerY,y", po::value<unsigned int>()->default_value(0), "choose y center of the projected image." )
    ("centerZ,z", po::value<unsigned int>()->default_value(1), "choose z center of the projected image." )
    ("width", po::value<unsigned int>()->default_value(100), "set the width of the resulting height Field image." )
    ("height", po::value<unsigned int>()->default_value(100), "set the height of the resulting height Field image." )
    ("heightFieldMaxScan", po::value<unsigned int>()->default_value(255), "set the maximal scan deep." )
    ("setBackgroundLastDepth", "change the default background (black with the last filled intensity).")
    ("rescaleInputMin", po::value<DGtal::int64_t>()->default_value(0), "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).")
    ("rescaleInputMax", po::value<DGtal::int64_t>()->default_value(255), "max value used to rescale the input intensity (to avoid basic cast into 8 bits image).");
  
  
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);    
  if( !parseOK || vm.count("help")||argc<=1)
    {
      std::cout << "Usage: " << argv[0] << " [input] [output]\n"
		<< "Convert volumetric  file into a projected 2D image given from a normal direction N and from a starting point P. The 3D volume is scanned in this normal direction N starting from P with a step 1. If the intensity of the 3d point is inside the given thresholds its 2D gray values are set to the current scan number."
		<< general_opt << "\n";
      std::cout << "Example:\n"
		<< "vol2heightfield -i ${DGtal}/examples/samples/lobster.vol -m 60 -M 500  --nx 0 --ny 0.7 --nz -1 -x 150 -y 0 -z 150 --width 300 --height 300 --heightFieldMaxScan 350  -o resultingHeightMap.pgm \n";
      return 0;
    }
  
  if(! vm.count("input") ||! vm.count("output"))
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;      
      return 0;
    }
  
  string inputFilename = vm["input"].as<std::string>();
  string outputFilename = vm["output"].as<std::string>();
  DGtal::int64_t rescaleInputMin = vm["rescaleInputMin"].as<DGtal::int64_t>();
  DGtal::int64_t rescaleInputMax = vm["rescaleInputMax"].as<DGtal::int64_t>();

  trace.info() << "Reading input file " << inputFilename ; 

  typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
  Image3D inputImage =  GenericReader< Image3D >::importWithValueFunctor( inputFilename,RescalFCT(rescaleInputMin,
                                                                                                  rescaleInputMax,
                                                                                                  0, 255) );


  trace.info() << " [done] " << std::endl ; 
  
  std::ofstream outStream;
  outStream.open(outputFilename.c_str());
  int minTh = vm["thresholdMin"].as<int>();
  int maxTh = vm["thresholdMax"].as<int>();
  
  trace.info() << "Processing image to output file " << outputFilename ; 
  
  unsigned int widthImageScan = vm["height"].as<unsigned int>();
  unsigned int heightImageScan = vm["width"].as<unsigned int>();
  unsigned int maxScan = vm["heightFieldMaxScan"].as<unsigned int>();
  if(maxScan > std::numeric_limits<Image2D::Value>::max()){
    maxScan = std::numeric_limits<Image2D::Value>::max();
    trace.warning()<< "value --setBackgroundLastDepth outside mox value of image. Set to max value:" << maxScan << std::endl; 
  }
  
  unsigned int centerX = vm["centerX"].as<unsigned int>();
  unsigned int centerY = vm["centerY"].as<unsigned int>();
  unsigned int centerZ = vm["centerZ"].as<unsigned int>();

  double nx = vm["nx"].as<double>();
  double ny = vm["ny"].as<double>();
  double nz = vm["nz"].as<double>();
  
  
  Image2D::Domain aDomain2D(DGtal::Z2i::Point(0,0), 
                          DGtal::Z2i::Point(widthImageScan, heightImageScan));
  Z3i::Point ptCenter (centerX, centerY, centerZ);
  Z3i::RealPoint normalDir (nx, ny, nz);
  Image2D resultingImage(aDomain2D);
  
  for(Image2D::Domain::ConstIterator it = resultingImage.domain().begin(); 
      it != resultingImage.domain().end(); it++){
    resultingImage.setValue(*it, 0);
  }
  DGtal::functors::Identity idV;
  
  unsigned int maxDepthFound = 0;
  for(unsigned int k=0; k < maxScan; k++){
    DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedder(inputImage.domain(), 
                                                                        ptCenter+normalDir*k,
                                                                        normalDir,
                                                                        widthImageScan);
    ImageAdapterExtractor extractedImage(inputImage, aDomain2D, embedder, idV);
    for(Image2D::Domain::ConstIterator it = extractedImage.domain().begin(); 
        it != extractedImage.domain().end(); it++){
      if( resultingImage(*it)== 0 &&  extractedImage(*it) < maxTh &&
          extractedImage(*it) > minTh){
        maxDepthFound = k;
        resultingImage.setValue(*it, maxScan-k);
      }
    }    
  }
  if (vm.count("setBackgroundLastDepth")){
    for(Image2D::Domain::ConstIterator it = resultingImage.domain().begin(); 
        it != resultingImage.domain().end(); it++){
      if( resultingImage(*it)== 0 ){
        resultingImage.setValue(*it, maxScan-maxDepthFound);
      }
    }
  }   
  
  resultingImage >> outputFilename;

  trace.info() << " [done] " << std::endl ;   
  return 0;  
}




