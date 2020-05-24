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

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;



/**
   @page vol2sdp vol2sdp
   @brief  Extracts digital points from 3d vol files.

   @b Usage: vol2sdp [input] [output]

   @b Allowed @b options @b are:

   @code
   
   vol2sdp -i ${DGtal}/examples/samples/lobster.vol -o volumeList.sdp
   Usage: ./converters/vol2sdp [OPTIONS] 1 [2]

   Positionals:
   1 TEXT:FILE REQUIRED        vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd) or sdp (sequence of discrete points). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
   2 TEXT=result.sdp           sequence of discrete point file (.sdp)

   Options:
   -h,--help                   Print this help message and exit
   -i,--input TEXT:FILE REQUIRED
                              vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd) or sdp (sequence of discrete points). 
                              For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
  -o,--output TEXT=result.sdp sequence of discrete point file (.sdp)
  -e,--exportImageValues      option to export also the image value of the voxel in a fourth field.
  -m,--thresholdMin INT=128   threshold min (excluded) to define binary shape.
  -M,--thresholdMax INT=255   threshold max (included) to define binary shape.
  --rescaleInputMin INT=0     min value used to rescale the input intensity (to avoid basic cast into 8  bits image).
  --rescaleInputMax INT=255   max value used to rescale the input intensity (to avoid basic cast into 8  bits image).

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
  CLI::App app;
  std::string inputFilename;
  std::string outputFilename {"result.sdp"};
  bool exportImageValues {false};
  int thresholdMin {128};
  int thresholdMax {255};  
  DGtal::int64_t rescaleInputMin {0};
  DGtal::int64_t rescaleInputMax {255};
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  
  // parse command line using CLI ----------------------------------------------
  app.description("Convert volumetric  file into a digital set of points from a given threshold.\n vol2sdp -i ${DGtal}/examples/samples/lobster.vol -o volumeList.sdp");
  app.add_option("-i,--input,1", inputFilename,  "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd) or sdp (sequence of discrete points). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("--output,-o,2", outputFilename, "sequence of discrete point file (.sdp)", true);
  app.add_flag("--exportImageValues,-e",exportImageValues, "option to export also the image value of the voxel in a fourth field.");
  app.add_option("--thresholdMin,-m", thresholdMin, "threshold min (excluded) to define binary shape.", true);
  app.add_option("--thresholdMax,-M", thresholdMax, "threshold max (included) to define binary shape.", true);
  app.add_option("--rescaleInputMin", rescaleInputMin, "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);
  app.add_option("--rescaleInputMax", rescaleInputMax, "max value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
  Image3D inputImage =  GenericReader< Image3D >::importWithValueFunctor( inputFilename,RescalFCT(rescaleInputMin,
                                                                                                  rescaleInputMax,
                                                                                                  0, 255) );
  trace.info() << " [done] " << std::endl ; 
  std::ofstream outStream;
  outStream.open(outputFilename.c_str());

  trace.info() << "Processing image to output file " << outputFilename ; 
  //Processing all points
  outStream << "# sdp file generate by vol2sdp with source vol:" << inputFilename
            << " and threshold min: " <<  thresholdMin << " max:" << thresholdMax << std::endl;
  outStream << "# format: x y z ";
  if(exportImageValues){
    outStream << " image_value";
  }
  outStream << std::endl;
  

  
  for(Image3D::Domain::ConstIterator it=inputImage.domain().begin(); it != inputImage.domain().end(); ++it){
    if(inputImage(*it) >= thresholdMin && inputImage(*it) <= thresholdMax ){
      outStream << (*it)[0] << " " << (*it)[1] << " " << (*it)[2];
      if(exportImageValues){
        outStream << " " << (unsigned int) inputImage(*it);
      }
      
      outStream << std::endl;
    }
  }
  outStream.close();
  trace.info() << " [done] " << std::endl ;   
  return EXIT_SUCCESS;
  
}

