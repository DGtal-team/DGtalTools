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

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;

/**
 @page raw2HDF5 raw2HDF5
 @brief  Converts a 3D 8-bit raw file to HDF5.

@b Usage: raw2HDF5 [input] [output]

@b Allowed @b options @b are:

@code
  -h [ --help ]             display this message
  -i [ --input ] arg        Input raw file.
  -o [ --output ] arg       Output HDF5 filename.
  -x [ --x ] arg            x extent.
  -y [ --y ] arg            y extent.
  -z [ --z ] arg            z extent.

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

  // parse command line ----------------------------------------------
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>(), "Input raw file." )
    ( "output,o", po::value<string>(),"Output HDF5 filename." )
    ( "x,x", po::value<unsigned int>(),"x extent." )
    ( "y,y", po::value<unsigned int >(),"y extent." )
    ( "z,z", po::value<unsigned int>(),"z extent." );
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }

  po::notify ( vm );
  if (!parseOK || vm.count ( "help" ) ||argc<=1 )
    {
      trace.info() << "Converts a 3D 8-bit raw file to HDF5."<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\traw2HDF5 -x 128 -y 128 -z 128 --input <RawFileName> --output <HDF5OutputFileName> "<<std::endl
                   << general_opt << "\n";
      return 0;
    }

  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  std::string filename = vm["input"].as<std::string>();
  if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
  std::string outputFileName = vm["output"].as<std::string>();
 if ( ! ( vm.count ( "x" ) ) ) missingParam ( "--x" );
  unsigned int x =  vm["x"].as<unsigned int>();
 if ( ! ( vm.count ( "y" ) ) ) missingParam ( "--y" );
  unsigned int y =  vm["y"].as<unsigned int>();
 if ( ! ( vm.count ( "z" ) ) ) missingParam ( "--z" );
  unsigned int z =  vm["z"].as<unsigned int>();
 

  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;

  MyImageC  imageC = RawReader< MyImageC >::importRaw8 ( filename, Z3i::Vector(x,y,z) );
  bool res =  HDF5Writer< MyImageC>::exportHDF5_3D(outputFileName, imageC, "/UInt8Array3D");

  if (res)
    return 0;
  else 
    {
      trace.error()<< "Error while exporting the HDF5 file."<<std::endl;
      return 1;
    }
}
