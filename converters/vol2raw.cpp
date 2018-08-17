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
 * @file vol2raw.cpp
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
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/RawWriter.h>
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
 @page vol2raw vol2raw
 @brief Converts a vol to a 8-bit raw file.

@b Usage: vol2raw [input] [output]

@b Allowed @b options @b are:

@code
    -h [ --help ]         display this message.
    -i [ --input ] arg    vol file (.vol, .longvol .p3d, .pgm3d and if 
                          WITH_ITK is selected: dicom, dcm, mha, mhd). For
                          longvol, dicom, dcm, mha or mhd formats, the 
                          input values are linearly scaled between 0 and 
                          255.
    -o [ --output ] arg   Output filename.
    --rescaleInputMin arg (=0)   min value used to rescale the input intensity 
                          (to avoid basic cast into 8  bits image).
    --rescaleInputMax arg (=255) max value used to rescale the input intensity 
                          (to avoid basic cast into 8 bits image).

@endcode

@b Example:
@code 
   vol2raw -i input.vol -o out.raw
@endcode

@see
@ref vol2raw.cpp

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
    ("input,i", po::value<std::string>(), "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
    ( "output,o", po::value<string>(),"Output filename." )
    ("rescaleInputMin", po::value<DGtal::int64_t>()->default_value(0), "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).")
    ("rescaleInputMax", po::value<DGtal::int64_t>()->default_value(255), "max value used to rescale the input intensity (to avoid basic cast into 8 bits image).");
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
      trace.info() << "Convert a vol to a 8-bit raw file."<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\tvol2raw --input <volFileName> --o <RawOutputFileName> "<<std::endl
                   << general_opt << "\n";
      return 0;
    }

  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  std::string filename = vm["input"].as<std::string>();
  if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
  std::string outputFileName = vm["output"].as<std::string>();
  DGtal::int64_t rescaleInputMin = vm["rescaleInputMin"].as<DGtal::int64_t>();
  DGtal::int64_t rescaleInputMax = vm["rescaleInputMax"].as<DGtal::int64_t>();
  
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;
  typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
  MyImageC imageC =  GenericReader< MyImageC >::importWithValueFunctor( filename ,RescalFCT(rescaleInputMin,
                                                                                               rescaleInputMax,
                                                                                               0, 255) );

  
  
  bool res =  RawWriter< MyImageC >::exportRaw8(outputFileName, imageC);
  trace.info() << "Raw export done, image dimensions: "  << imageC.domain().upperBound()[0]-imageC.domain().lowerBound()[0]+1
               << " " << imageC.domain().upperBound()[1]-imageC.domain().lowerBound()[1]+1
               << " " << imageC.domain().upperBound()[2]-imageC.domain().lowerBound()[2]+1 << std::endl;
    
  if (res) return 0; else return 1;
}
