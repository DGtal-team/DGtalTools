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
 * @file visuDistanceTransform.cpp
 * @ingroup surfaceTools
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2012/07/08
 *
 *
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"

#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/Display3D.h"
#include "DGtal/shapes/Mesh.h"

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;

/**
 @page ofs2off ofs2off
 @brief  Convert OFS file into OFF mesh format.

@b Usage: ofs2off [input] [output]

@b Allowed @b options @b are:

@code
Positionals:
  1 TEXT:FILE REQUIRED                  ofs file (.ofs).

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         ofs file (.ofs).
  -o,--output TEXT                      ofs file (.off)
@endcode

@b Example:
@code
  $ ofs2off -i input.ofs -o output.off 
@endcode

@see ofs2off.cpp

*/



int main( int argc, char** argv )
{
// parse command line using CLI ----------------------------------------------
   CLI::App app;
   std::string inputFileName;
   std::string outputFileName {"result.off"};

   app.description("Convert OFS file into OFF mesh format.");
   app.add_option("-i,--input,1", inputFileName, "ofs file (.ofs)." )
    ->required()
    ->check(CLI::ExistingFile);
   app.add_option("-o,--output,2", outputFileName, "ofs file (.off)");
   
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------
   
   // We store the colors
   Mesh<Z3i::RealPoint> anImportedMesh(true);
   bool import = anImportedMesh << inputFileName;
   bool exported = anImportedMesh >> outputFileName;
   if(!import || !exported){
     trace.info() << "Conversion failed: " <<  (exported? "Reading OFS failed. ":  "Export OFF failed. ") << std::endl;
     return 0;
   }
  trace.info() << "[done]. "<< std::endl;
  return 0;
}
