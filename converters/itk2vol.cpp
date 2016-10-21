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
 * @file itk2vol.cpp
 * @ingroup converters
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2015/06/06
 *
 * 
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/readers/ITKReader.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;

/**
 @page itk2vol itk2vol
 @brief  Converts itk file into a volumetric file (.vol, .pgm3d).

@b Usage: itk2vol [input] [output]

@b Allowed @b options @b are:

@code
  -h [ --help ]           display this message
  -i [ --input ] arg      Any file format in the ITK library (mhd, mha, ...) 
  -o [ --output ] arg     volumetric file (.vol, .pgm3d) 
  --inputMin arg (=-1000) set minimum density threshold on Hounsfield scale
  --inputMax arg (=3000)  set maximum density threshold on Hounsfield scale
@endcode

@b Example:
@code
$itk2vol -i image.mhd --dicomMin -500 --dicomMax -100 -o sample.vol 
@endcode

@see itk2vol.cpp

*/


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3DChar;
  typedef ImageContainerBySTLVector < Z3i::Domain,  int > Image3D;
  
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "Any file format in the ITK library (mhd, mha, ...) " )
    ("output,o", po::value<std::string>(), "volumetric file (.vol, .pgm3d) " )
    ("inputMin", po::value<int>()->default_value(-1000), "set minimum density threshold on Hounsfield scale")
    ("inputMax", po::value<int>()->default_value(3000), "set maximum density threshold on Hounsfield scale");

  
  
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
		<< "Converts itk file into a volumetric file (.vol, .pgm3d). "
		<< general_opt << "\n";
      std::cout << "Example:\n"
		<< "itk2vol -i image.mhd --inputMin -500 --inputMax -100 -o sample.vol \n";
      return 0;
    }
  
  if(! vm.count("input") ||! vm.count("output"))
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;      
      return 0;
    }

  
  string inputFilename = vm["input"].as<std::string>();
  string outputFilename = vm["output"].as<std::string>();
  int inputMin = vm["inputMin"].as<int>();
  int inputMax = vm["inputMax"].as<int>();
  typedef DGtal::functors::Rescaling<int ,unsigned char > RescalFCT;
  
  trace.info() << "Reading input input file " << inputFilename ; 
  Image3D inputImage = ITKReader< Image3D  >::importITK(inputFilename);
  trace.info() << " [done] " << std::endl ; 
  trace.info() << " converting into vol file... " ; 
  RescalFCT rescaleCustom(inputMin, inputMax, 0, 255);

  DGtal::GenericWriter<Image3D, 3, unsigned char, RescalFCT>::exportFile(outputFilename, inputImage, "UInt8Array3D", rescaleCustom);
  

  trace.info() << " [done] " << std::endl ;   


  return 0;
  
}


