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

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;

/**
 @page HDF52vol HDF52vol
 @brief Converts a 3D 8-bit HDF5 file to vol.

@b Usage: HDF52vol [input] [output]

@b Allowed @b options @b are:

@code
  -h [ --help ]                     display this message
  -i [ --input ] arg                the input FreemanChain file name
  -o [ --output ] arg (=result.pgm)  the output filename
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

  // parse command line ----------------------------------------------
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>(), "Input HDF5 file." )
    ( "output,o", po::value<string>(),"Output vol filename." );
  
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
      trace.info() << "Convert a 3D 8-bit HDF5 file to vol."<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\tHDF52vol --input <HDF5FileName> --output <VolOutputFileName> "<<std::endl
                   << general_opt << "\n";
      return 0;
    }

  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  std::string filename = vm["input"].as<std::string>();
  if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
  std::string outputFileName = vm["output"].as<std::string>();
 

  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;

  MyImageC imageC = HDF5Reader< MyImageC >::importHDF5_3D( filename, "/UInt8Array3D" );
  bool res = VolWriter< MyImageC>::exportVol(outputFileName, imageC);

  if (res)
    return 0;
  else 
    {
      trace.error()<< "Error while exporting the volume."<<std::endl;
      return 1;
    }
}
