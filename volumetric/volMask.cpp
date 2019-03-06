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
 * @file
 * @ingroup volumetric
 * @author Bertrand Kerautret (\c bertrand.kerautret@univ-lyon2.fr )
 * 
 *
 * @date 2019/03/01
 *
 * Source file of the tool volMask
 *
 * This file is part of the DGtal library/DGtalTools Project.
 */

///////////////////////////////////////////////////////////////////////////////
#include "DGtal/base/Common.h"
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace DGtal;
///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;


/**
 @page volMask volMask
 
 @brief  Description of the tool...

 @b Usage:   volMask [input]

 @b Allowed @b options @b are :
 
 @code
  -h [ --help ]           display this message
  -i [ --input ] arg      an input file... 
  -p [ --parameter] arg   a double parameter...
 @endcode

 @b Example: 

 @code
   	volMask -i  $DGtal/examples/samples/....
 @endcode

 @image html resvolMask.png "Example of result. "

 @see
 @ref volMask.cpp

 */

typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;


int main( int argc, char** argv )
{
  // parse command line -------------------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string >(), "an input vol file. " )
    ("mask,a", po::value<std::string >(), "the mask image." )
    ("output,o", po::value<std::string >(), "the output masked image." )
  ("maskValue,m", po::value<Image3D::Value>()->default_value(1), "the masking value." );


  bool parseOK=true;
  po::variables_map vm;
  try
    {
      po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }
  catch(const std::exception& ex)
    {
      parseOK=false;
      trace.info()<< "Error checking program options: "<< ex.what()<< endl;
    }
  

  // check if min arguments are given and tools description ------------------
  po::notify(vm);
  if( !parseOK || vm.count("help")||argc<=1)
    {
      std::cout << "Usage: " << argv[0] << " [input]\n"
                << "The tools description... \n"
                << general_opt << "\n"
                << "Typical use example:\n \t volMask -i ... \n";
      return 0;
    }  
  if(! vm.count("input"))
    {
      trace.error() << " The file name was not defined" << endl;
      return 1;
    }



  //  recover the  args ----------------------------------------------------
  string inputFileName = vm["input"].as<string>();
  string maskFileName = vm["mask"].as<string>();
  string outputFileName = vm["output"].as<string>();

  Image3D::Value maskValue = vm["maskValue"].as<Image3D::Value>();

  trace.info() << "Reading input image...";
  Image3D inputImage = DGtal::GenericReader<Image3D>::import(inputFileName);
  trace.info() << "[done]"<< std::endl;
  trace.info() << "Reading mask image...";
  Image3D maskImage = DGtal::GenericReader<Image3D>::import(maskFileName);
  trace.info() << "[done]"<< std::endl;

  // Some nice processing  --------------------------------------------------
  Image3D::Domain d =  inputImage.domain();
  // First step getting the bounding box of the domain:

  Z3i::Point minP = inputImage.domain().upperBound();
  Z3i::Point maxP = inputImage.domain().lowerBound();
  
  Z3i::Point::Iterator minIt;
  Z3i::Point::Iterator maxIt;

  for(const auto &p: inputImage.domain())
  {
    minIt = minP.begin();
    maxIt = maxP.begin();
    maxIt = maxP.begin();
    if( maskImage(p) ) // no noise on mask image
	  {
      for(auto pIt=p.begin(); pIt!=p.end();pIt++ )
      {
        if( *pIt < *minIt ){*minIt = *pIt;}
        if( *pIt > *maxIt ){*maxIt = *pIt;}
        minIt++;
        maxIt++;
      }
	  }
  }
    
  // offset to avoid problems on borders
  Z3i::Point offset(5,5,5);
  minP -= offset;
  maxP += offset;
  
  trace.info() << "sub-domain:" << minP << " " << maxP << std::endl;

  Image3D outputImage( Image3D::Domain(minP,maxP) );
  
  // Second step: masking source image
  for (const auto &p: outputImage.domain())
  {
    if (maskImage(p) ==  maskValue)
    {
      outputImage.setValue(p, inputImage(p) );
    }
  }

  trace.info() << "writing output image...";
  GenericWriter<Image3D>::exportFile(outputFileName, outputImage);
  trace.info() << "[Done]" << std::endl;
  return 0;
}

