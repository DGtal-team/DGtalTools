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


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;

/**
 @page convertVol convertVol
 @brief Converts volumetric file into volumetric file from different formats (pgm3d, vol, longvol). This tool can also be used to upgrade a Version-2 Vol or Longvol file to the new (compressed) Version-3.


@b Usage: convertVol [input] [output]

@b Allowed @b options @b are:

@code
  -h [ --help ]         display this message
  -i [ --input ] arg    volumetric file (.pgm3d, .vol, .longvol) 
  -o [ --output ] arg   volumetric file (.pgm3d, .vol, .longvol) 
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


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char> Image3D;

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "volumetric file (.pgm3d, .vol, .longvol) " )
    ("output,o", po::value<std::string>(), "volumetric file (.pgm3d, .vol, .longvol) " );
    
  
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
		<< "Convert volumetric file into volumetric file from different formats (pgm3d, vol, longvol) "
		<< general_opt << "\n";
      std::cout << "Example:\n"
		<< "convertVol -i ${DGtal}/examples/samples/lobster.vol -o convertedVol.p3d \n";
      return 0;
    }
  
  if(! vm.count("input")||! vm.count("output"))
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;      
      return 0;
    }

  
  string inputFilename = vm["input"].as<std::string>();
  string outputFilename = vm["output"].as<std::string>();
  
  trace.info() << "Reading input file " << inputFilename ; 
  Image3D inputImage = DGtal::GenericReader<Image3D>::import(inputFilename);
  trace.info() << " [done] " << std::endl ; 
  trace.info() << "Writing output file " << outputFilename ; 
  DGtal::GenericWriter<Image3D>::exportFile(outputFilename,  inputImage);
  trace.info() << " [done] " << std::endl ;   


  return 0;
  
}




