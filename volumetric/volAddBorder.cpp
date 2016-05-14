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
 * @file voAddBorder.cpp
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
#include <DGtal/io/readers/VolReader.h>
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
 @page volAddBorder volAddBorder
 
 @brief Adds a border of one voxel with value 0 around a vol file.

 @b Usage: 	volAddBorder --input \<volFileName\> --o \<volOutputFileName\> 


 @b Allowed @b options @b are : 
 @code
  -h [ --help ]         display this message.
  -i [ --input ] arg    Input vol file.
  -o [ --output ] arg   Output filename.

 @endcode

 @b Example: 

 @code
 $  volAddBorder -i $DGtal/examples/samples/Al.100.vol -o Al.100border.vol
 @endcode
 

 You should obtain an resulting image with domain:
-1, -1, -1 x 100, 100, 100 instead  0, 0, 0 x 99, 99, 99

 
 @see
 @ref volAddBorder.cpp

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
  po::options_description general_opt ( "Allowed options are " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>(), "Input vol file." )
    ( "output,o", po::value<string>(),"Output filename." );
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify ( vm );
  if ( !parseOK || vm.count ( "help" ) ||argc<=1 )
    {
      trace.info() << "Add a border of one voxel with value 0 around a vol file."<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\tvolAddBorder --input <volFileName> --o <volOutputFileName> "<<std::endl
                   << general_opt << "\n"
                   << "Example: \n \t volAddBorder -i $DGtal/examples/samples/Al.100.vol -o Al.100border.vol";
      return 0;
    }

  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  std::string filename = vm["input"].as<std::string>();
  if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
  std::string outputFileName = vm["output"].as<std::string>();

  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;

  MyImageC  imageC = VolReader< MyImageC >::importVol ( filename );
  MyImageC  outputImage( Z3i::Domain( imageC.domain().lowerBound() - Vector().diagonal(1),
                                      imageC.domain().upperBound() + Vector().diagonal(1)));
  
  //Fast Copy
  for(MyImageC::Domain::ConstIterator it = imageC.domain().begin(),
        itend = imageC.domain().end(); it != itend; ++it)
    outputImage.setValue( *it , imageC(*it));
  

  bool res =  VolWriter< MyImageC>::exportVol(outputFileName, outputImage);
  if (res) return 0; else return 1;
}
