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
 * @file volIntensityScale.cpp
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2015/01/12
 *
 * 
 *
 * This file is part of the DGtalTools.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/VolReader.h>
#include <DGtal/io/writers/VolWriter.h>
#include "DGtal/io/writers/GenericWriter.h"
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ConstImageAdapter.h>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;


/**
 @page volIntensityScale volIntensityScale
 
 @brief Applies a linear rescaling of the image intensity from an input intensity interval [InMin, InMax] into an output interval [OutMin, OutMax].

 @b Usage: 	 volIntensityScale --input \<volFileName\> --output \<volOutputFileName\>  (both files can be independently in vol, pgm3D, p3d format)



 @b Allowed @b options @b are : 
 @code
  -h [ --help ]             display this message.
  -i [ --input ] arg        Input vol file.
  -o [ --output ] arg       volumetric output file (.vol, .pgm, .pgm3d, 
                            .longvol) 
  -m [ --inMin ] arg (=0)   the min value of the input image.
  -M [ --inMax ] arg (=255) the max value of the input image.
  --outMin arg (=0)         the min value of the output image.
  --outMax arg (=255)       the max value of the output image.

 @endcode

 @b Example: 

 In this example, we select all intensities included in [0,100] and scale them into [0, 255]: 

 @code
  $ volIntensityScale -i  ${DGtal}/examples/samples/lobster.vol  --inMin 0 --inMax 100 -o lobster0-100.vol 

 @endcode


 You should obtain such a result:
 @image html resVolIntensityScale.png "(a) source image, (b) intensity scaled image."
 
 @see
 @ref volIntensityScale.cpp

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
}


int main(int argc, char**argv)
{

  // parse command line ----------------------------------------------
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>(), "Input vol file." )
    ( "output,o", po::value<std::string>(), "volumetric output file (.vol, .pgm, .pgm3d, .longvol) " )
    ( "inMin,m", po::value<int>()->default_value(0), "the min value of the input image." )
    ( "inMax,M", po::value<int>()->default_value(255), "the max value of the input image." )
    ( "outMin", po::value<int>()->default_value(0), "the min value of the output image." )
    ( "outMax", po::value<int>()->default_value(255), "the max value of the output image." );
  
  bool parseOK=true;
  po::variables_map vm;
  
  try{
    po::store ( po::parse_command_line ( argc, argv, general_opt ), vm );
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
    //Parse options
  if (parseOK &&  ! ( vm.count ( "input" ) ) ) {parseOK=false; missingParam ( "--input" );};
  if (parseOK && ! ( vm.count ( "output" ) ) ) {parseOK=false; missingParam ( "--output" );};
  

  po::notify ( vm );
  if (!parseOK || vm.count ( "help" ))
    {
      trace.info() << "Apply a linear rescaling of the image intensity from an input intensity interval [InMin, InMax] into an output interval [OutMin, OutMax]." <<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\t volIntensityScale --input <volFileName> --output <volOutputFileName>  (both files can be independently in vol, pgm3D, p3d format)"<<std::endl
                   << general_opt << "\n";
      std::cout << "Example:\n"
		<< "volIntensityScale -i  ${DGtal}/examples/samples/lobster.vol  --inMin 0 --inMax 100 -o lobster0-100.vol \n";
      return 0;
    }



  std::string filename = vm["input"].as<std::string>();
  std::string outputFileName = vm["output"].as<std::string>();
  
  int inMin = vm["inMin"].as<int>();
  int inMax = vm["inMax"].as<int>();
  int outMin = vm["outMin"].as<int>();
  int outMax = vm["outMax"].as<int>();
    
  trace.beginBlock("Loading file");

  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;
  typedef DGtal::functors::Rescaling<int ,unsigned char > RescalFCT;
  MyImageC image = VolReader< MyImageC,  RescalFCT  >::importVol( filename, 
                                                                  RescalFCT(inMin,
                                                                            inMax,
                                                                            outMin, outMax));  
  trace.endBlock();

  
  trace.beginBlock("Exporting...");
  bool res =  GenericWriter<MyImageC>::exportFile(outputFileName, image);
  trace.endBlock();
  if (res) return 0; else return 1;   
}




