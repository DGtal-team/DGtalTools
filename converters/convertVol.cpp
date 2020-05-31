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
 * @file convertVol.cpp
 * @ingroup surfaceTools
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
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;

/**
 @page convertVol convertVol
 @brief Converts volumetric file into volumetric file from different formats (pgm3d, vol, longvol). This tool can also be used to upgrade a Version-2 Vol or Longvol file to the new (compressed) Version-3.

@b Usage: convertVol [input] [output]

@b Allowed @b options @b are:

@code
Positionals:
  1 TEXT:FILE REQUIRED                  volumetric file (.pgm3d, .vol, .longvol).

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         volumetric file (.pgm3d, .vol, .longvol).
  -o,--output TEXT                      volumetric file (.pgm3d, .vol, .longvol)

@endcode

@b Examples:
@code 
$ convertVol -i ${DGtal}/examples/samples/lobster.vol -o convertedVol.p3d 
@endcode

 To upgrade a "Version-2" vol file to a "Version-3" (default Vol writer):
 @code
 $ convertVol -i ${DGtal}/examples/samples/lobster.vol -o ${DGtal}/examples/samples/lobster.vol
 @endcode
 

@see convertVol.cpp

*/


int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char> Image3D;

  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.vol"};
  
  app.description("Convert volumetric file into volumetric file from different formats (pgm3d, vol, longvol)\n ");
  app.add_option("-i,--input,1", inputFileName, "volumetric file (.pgm3d, .vol, .longvol)." )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("-o,--output,2", outputFileName, "volumetric file (.pgm3d, .vol, .longvol)", true);

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  trace.info() << "Reading input file " << inputFileName ; 
  Image3D inputImage = DGtal::GenericReader<Image3D>::import(inputFileName);
  trace.info() << " [done] " << std::endl ; 
  trace.info() << "Writing output file " << outputFileName ; 
  DGtal::GenericWriter<Image3D>::exportFile(outputFileName,  inputImage);
  trace.info() << " [done] " << std::endl ;   

  return EXIT_SUCCESS;  
}




