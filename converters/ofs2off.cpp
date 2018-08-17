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

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;

/**
 @page ofs2off ofs2off
 @brief  Convert OFS file into OFF mesh format.

@b Usage: ofs2off [input] [output]

@b Allowed @b options @b are:

@code
  -h [ --help ]         display this message
  -i [ --input ] arg    ofs file (.ofs) 
  -o [ --output ] arg   ofs file (.off) 
@endcode

@b Example:
@code
  $ ofs2off -i input.ofs -o output.off 
@endcode

@see ofs2off.cpp

*/


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "ofs file (.ofs) " )
    ("output,o", po::value<std::string>(), "ofs file (.off) " );


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
		<< "Convert OFS file into OFF mesh format"
		<< general_opt << "\n";
      return 0;
    }

  if(! vm.count("input")||! vm.count("output"))
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;
      return 0;
    }


  string inputFilename = vm["input"].as<std::string>();
  string outputFilename = vm["output"].as<std::string>();

  // We store the colors
  Mesh<Display3D<>::BallD3D> anImportedMesh(true);
  bool import = anImportedMesh << inputFilename;
  bool exported = anImportedMesh >> outputFilename;
  if(!import || !exported){
    trace.info() << "Conversion failed: " <<  (exported? "Reading OFS failed. ":  "Export OFF failed. ") << std::endl;
    return 0;
  }



  trace.info() << "[done]. "<< std::endl;

  return 0;

}
