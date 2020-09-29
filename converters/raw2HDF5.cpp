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
 * @file raw2HDF5.cpp
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
#include <DGtal/io/readers/RawReader.h>
#include <DGtal/io/writers/HDF5Writer.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;
using namespace Z3i;



/**
 @page raw2HDF5 raw2HDF5
 @brief  Converts a 3D 8-bit raw file to HDF5.

@b Usage: raw2HDF5 [input] [output]

@b Allowed @b options @b are:

@code

Positionals:
  1 TEXT:FILE REQUIRED                  Input raw file.

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         Input raw file.
  -o,--output TEXT=result.hdf5           Output hdf5 filename.
  --x UINT REQUIRED                     x extent.
  --y UINT REQUIRED                     y extent.
  --z UINT REQUIRED                     z extent.

@endcode

@b Example:
@code
  $ raw2HDF5 -x 128 -y 128 -z 128 -i input.raw -o output.hd
@endcode


@see raw2HDF5.cpp

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
  std::string outputFileName {"result.hdf5"};
  unsigned int x, y, z;
   
  app.description("Converts a 3D 8-bit raw file to HDF5.\n Basic usage \n \traw2HDF5 -x 128 -y 128 -z 128 --input <RawFileName> --output <HDF5OutputFileName>");
   
  app.add_option("-i,--input,1", inputFileName, "Input raw file." )
    ->required()
    ->check(CLI::ExistingFile);
   app.add_option("-o,--output,2",outputFileName,"Output hdf5 filename.", true);
   app.add_option("--x,-x", x, "x extent." )
   ->required();
   app.add_option("--y,-y", y, "y extent." )
     ->required();
   app.add_option("--z,-z", z, "z extent." )
     ->required(); 
   
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------


  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;

  MyImageC  imageC = RawReader< MyImageC >::importRaw8 ( inputFileName, Z3i::Vector(x,y,z) );
  bool res =  HDF5Writer< MyImageC>::exportHDF5_3D(outputFileName, imageC, "/UInt8Array3D");


  if (res)
     return EXIT_SUCCESS;
  else 
    {
      trace.error()<< "Error while exporting the volume."<<std::endl;
      return EXIT_FAILURE;
    }

}
