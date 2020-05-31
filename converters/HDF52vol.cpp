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
 * @file HDF52vol.cpp
 * @author Martial Tola (\c martial.tola@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2013/09/11
 *
 *
 * This file is part of the DGtal library.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/HDF5Reader.h>
#include <DGtal/io/writers/VolWriter.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
 @page HDF52vol HDF52vol
 @brief Converts a 3D 8-bit HDF5 file to vol.

@b Usage: HDF52vol [input] [output]

@b Allowed @b options @b are:

@code

Positionals:
  1 TEXT:FILE REQUIRED                  volumetric file (.pgm3d, .vol, .longvol).

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         the input FreemanChain file name
  -o,--output TEXT                      the output filename
@endcode

@b Example:
@code
  $HDF52vol -i ${DGtal}/tests/samples/ex_image2.h5 -o out.vol
@endcode

@see HDF52vol.cpp

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
   std::string outputFileName {"result.vol"};

   app.description("Convert a 3D 8-bit HDF5 file to vol.");
   app.add_option("-i,--input,1", inputFileName, "Input HDF5 file." )
    ->required()
    ->check(CLI::ExistingFile);
   app.add_option("-o,--output,2", outputFileName, "Output vol filename.", true );

   app.get_formatter()->column_width(40);
   CLI11_PARSE(app, argc, argv);
   // END parse command line using CLI ----------------------------------------------
   
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;

  MyImageC imageC = HDF5Reader< MyImageC >::importHDF5_3D( inputFileName, "/UInt8Array3D" );
  bool res = VolWriter< MyImageC>::exportVol(outputFileName, imageC);

  if (res)
    return EXIT_SUCCESS;
  else 
    {
      trace.error()<< "Error while exporting the volume."<<std::endl;
      return EXIT_FAILURE;
    }
}
