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
 * @file raw2vol.cpp
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
#include <DGtal/io/readers/RawReader.h>
#include <DGtal/io/writers/VolWriter.h>
#include <DGtal/io/readers/VolReader.h>
#include <DGtal/io/writers/RawWriter.h>
#include <DGtal/helpers/StdDefs.h>
#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
 @page raw2vol raw2vol
 @brief  Converts a  8-bit raw file to vol.

@b Usage: raw2vol [input] [output]

@b Allowed @b options @b are:

@code
Allowed options are: :

Positionals:
  1 TEXT:FILE REQUIRED                  Input raw file.

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         Input raw file.
  -o,--output TEXT=result.vol           Output vol filename.
  --x UINT REQUIRED                     x extent.
  --y UINT REQUIRED                     y extent.
  --z UINT REQUIRED                     z extent.
@endcode

@b Example:
@code
  $ raw2vol -x 128 -y 128 -z 128  -i input.raw  -o output.vol 

@endcode

@see raw2vol.cpp

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
   unsigned int x, y, z;
   app.description("Converts a  8-bit raw file to  vol.\n Basic example:\n \t raw2vol -x 128 -y 128 -z 128 --input <RawFileName> --output <VolOutputFileName>");
   app.add_option("-i,--input,1", inputFileName, "Input raw file." )
    ->required()
    ->check(CLI::ExistingFile);
   app.add_option("-o,--output,2",outputFileName,"Output vol filename.", true);
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
  MyImageC  imageC = RawReader< MyImageC >::importRaw8 ( inputFileName, Z3i::Vector(x,y,z));
  bool res =  VolWriter< MyImageC>::exportVol(outputFileName, imageC);

  if (res)
     return EXIT_SUCCESS;
  else 
    {
      trace.error()<< "Error while exporting the volume."<<std::endl;
      return EXIT_FAILURE;
    }
}
