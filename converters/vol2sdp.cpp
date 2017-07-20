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
 * @file vol2sdp.cpp
 * @ingroup Converters
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2013/05/07
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
#include "DGtal/io/readers/GenericReader.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;


/**
   @page vol2sdp vol2sdp
   @brief  Extracts digital points from 3d vol files.

   @b Usage: vol2sdp [input] [output]

   @b Allowed @b options @b are:

   @code
   -h [ --help ]                    display this message
   -i [ --input ] arg               vol file (.vol, .longvol .p3d, .pgm3d and if
                                    WITH_ITK is selected: dicom, dcm, mha, mhd) 
                                    or sdp (sequence of discrete points). For 
                                    longvol, dicom, dcm, mha or mhd formats, the
                                    input values are linearly scaled between 0 
                                    and 255.
   -o [ --output ] arg              sequence of discrete point file (.sdp) 
   -e [ --exportImageValues ]       option to export also the image value of the
   voxel in a fourth field.
   -m [ --thresholdMin ] arg (=128) min threshold (default 128)
   -M [ --thresholdMax ] arg (=255) max threshold (default 255)
   --rescaleInputMin arg (=0)       min value used to rescale the input 
   intensity (to avoid basic cast into 8  bits 
   image).
   --rescaleInputMax arg (=255)     max value used to rescale the input 
   intensity (to avoid basic cast into 8 bits 
   image).


   @endcode

   @b Example:
   @code 
   $ vol2sdp -i ${DGtal}/examples/samples/lobster.vol -o volumeList.sdp -m 70
   # Visualisation:
   $ 3dSDPViewer -i volumeList.sdp
   @endcode

   You should obtain such a visualization:
   @image html resVol2sdp.png "resulting visualisation."

   @see
   @ref vol2sdp.cpp

*/

int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd) or sdp (sequence of discrete points). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
    ("output,o", po::value<std::string>(), "sequence of discrete point file (.sdp) " )
    ("exportImageValues,e","option to export also the image value of the voxel in a fourth field.")
    ("thresholdMin,m", po::value<int>()->default_value(128), "min threshold (default 128)" )
    ("thresholdMax,M", po::value<int>()->default_value(255), "max threshold (default 255)" )
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
		<< "Convert volumetric  file into a digital set of points from a given threshold."
		<< general_opt << "\n";
      std::cout << "Example:\n"
		<< "vol2sdp -i ${DGtal}/examples/samples/lobster.vol -o volumeList.sdp \n";
      return 0;
    }
  
  if(! vm.count("input") ||! vm.count("output"))
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;      
      return 0;
    }

  
  string inputFilename = vm["input"].as<std::string>();
  string outputFilename = vm["output"].as<std::string>();
  
  trace.info() << "Reading input file " << inputFilename ; 
  DGtal::int64_t rescaleInputMin = vm["rescaleInputMin"].as<DGtal::int64_t>();
  DGtal::int64_t rescaleInputMax = vm["rescaleInputMax"].as<DGtal::int64_t>();

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
  //Processing all points
  outStream << "# sdp file generate by vol2sdp with source vol:" << inputFilename << " and threshold min: " <<  minTh << " max:" << maxTh << std::endl;
  outStream << "# format: x y z ";
  if(vm.count("exportImageValues")){
    outStream << " image_value";
  }
  outStream << std::endl;
  
  
  for(Image3D::Domain::ConstIterator it=inputImage.domain().begin(); it != inputImage.domain().end(); ++it){
    if(inputImage(*it) >= minTh && inputImage(*it) <= maxTh ){
      outStream << (*it)[0] << " " << (*it)[1] << " " << (*it)[2];
      if(vm.count("exportImageValues")){
        outStream << " " << (unsigned int) inputImage(*it);
      }
      
      outStream << std::endl;
    }
  }
  outStream.close();

  trace.info() << " [done] " << std::endl ;   


  return 0;
  
}




