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
 * @file vol2raw.cpp
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/05/01
 *
 *
 * This file is part of the DGtal library.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/RawWriter.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
 @page vol2raw vol2raw
 @brief Converts a vol to a 8-bit raw file.

@b Usage: vol2raw [input] [output]

@b Allowed @b options @b are:

@code

Usage: ./converters/vol2raw [OPTIONS] 1 [2]

Positionals:
  1 TEXT:FILE REQUIRED                  vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
  2 TEXT                                output file (.raw).

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
  -o,--output TEXT                      output file (.raw).
  --rescaleInputMin INT=0               min value used to rescale the input intensity (to avoid basic cast into 8  bits image).
  --rescaleInputMax INT=255             max value used to rescale the input intensity (to avoid basic cast into 8  bits image).

@endcode

@b Example:
@code 
   vol2raw  input.vol  out.raw
@endcode

@see
@ref vol2raw.cpp

*/


/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam ( std::string param )
{
  trace.error() <<" Parameter: "<<param<<" is required..";
  trace.info() <<std::endl;
  exit ( 1 );
}


int main(int argc, char**argv)
{

  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.raw"};
  DGtal::int64_t rescaleInputMin {0};
  DGtal::int64_t rescaleInputMax {255};
  
  app.description("Convert a vol to a 8-bit raw file.\n Example: vol2raw  ${DGtal}/examples/samples/lobster.vol res.raw \n");
  app.add_option("-i,--input,1", inputFileName, "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("--output,-o,2",outputFileName ,"output file (.raw).");
  app.add_option("--rescaleInputMin", rescaleInputMin, "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);
  app.add_option("--rescaleInputMax", rescaleInputMax, "max value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);
 
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;
  typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
  MyImageC imageC =  GenericReader< MyImageC >::importWithValueFunctor( inputFileName ,RescalFCT(rescaleInputMin,
                                                                                                 rescaleInputMax,
                                                                                                 0, 255) );

  bool res =  RawWriter< MyImageC >::exportRaw8(outputFileName, imageC);
  trace.info() << "Raw export done, image dimensions: "  << imageC.domain().upperBound()[0]-imageC.domain().lowerBound()[0]+1
               << " " << imageC.domain().upperBound()[1]-imageC.domain().lowerBound()[1]+1
               << " " << imageC.domain().upperBound()[2]-imageC.domain().lowerBound()[2]+1 << std::endl;
    
  if (res)
     return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}
