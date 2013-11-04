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
 * @file vol2sdp.cpp
 * @ingroup conerters
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
#include <fstream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/readers/VolReader.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input-file,i", po::value<std::string>(), "volumetric file (.vol) " )
    ("output-file,o", po::value<std::string>(), "sequence of discrete point file (.sdp) " )
    ("thresholdMin,m", po::value<int>(), "min threshold (default 128)" )
    ("thresholdMax,M", po::value<int>(), "max threshold (default 255)" );
  
  
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
      std::cout << "Usage: " << argv[0] << " [input-file] [output-file]\n"
		<< "Convert volumetric  file into a digital set of points from a given threshold."
		<< general_opt << "\n";
      std::cout << "Example:\n"
		<< "vol2sdp -i ${DGtal}/examples/samples/lobster.vol -o volumeList.p3d \n";
      return 0;
    }
  
  if(! vm.count("input-file") ||! vm.count("output-file"))
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;      
      return 0;
    }

  
  string inputFilename = vm["input-file"].as<std::string>();
  string outputFilename = vm["output-file"].as<std::string>();
  
  trace.info() << "Reading input file " << inputFilename ; 
  Image3D inputImage = DGtal::VolReader<Image3D>::importVol(inputFilename);
  trace.info() << " [done] " << std::endl ; 
  std::ofstream outStream;
  outStream.open(outputFilename.c_str());
  int minTh = vm["thresholdMin"].as<int>();
  int maxTh = vm["thresholdMax"].as<int>();

  trace.info() << "Processing image to output file " << outputFilename ; 
  //Processing all points
  for(Image3D::Domain::ConstIterator it=inputImage.domain().begin(); it != inputImage.domain().end(); ++it){
    if(inputImage(*it) >= minTh && inputImage(*it) <= maxTh ){
      outStream << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << std::endl;
    }
  }
 outStream.close();

  trace.info() << " [done] " << std::endl ;   


  return 0;
  
}




