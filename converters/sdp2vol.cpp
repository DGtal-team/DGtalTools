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
 * @date 2013/07/21
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
#include "DGtal/io/writers/GenericWriter.h"

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
    ("sdp,s", po::value<std::string>(), "Sequence of 3d Discrete points (.sdp) " )
    ("output-file,o", po::value<std::string>(), "Vol file  (.vol, .longvol, .pgm3d) " )
    ("foregroundVal,f", po::value<int>()->default_value(128), "value which will represent the foreground object in the resulting image (default 128)")
    ("invertY", "Invert the Y axis (image flip in the y direction)")
    ("backgroundVal,b", po::value<int>()->default_value(0), "value which will represent the background outside the  object in the resulting image (default 0)")
    ("domain,d",  po::value<std::vector <int> >()->multitoken(), "The domain of the resulting image xmin ymin zmin xmax ymax zmax ");
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);    
  if( !parseOK || vm.count("help"))
    {
      std::cout << "Usage: " << argv[0] << " [input-file] [output-file]\n"
		<< "Convert volumetric  file into a digital set of points from a given threshold."
		<< general_opt << "\n";
      std::cout << "Example:\n"
		<< "vol2sdp -i ${DGtal}/examples/samples/lobster.vol -o volumeList.p3d \n";
      return 0;
    }
  if(! vm.count("sdp") ||! vm.count("output-file") || !vm.count("domain") )
    {
      trace.error() << " Input/ output filename and domain are needed to be defined" << endl;      
      return 0;
    }
  
  std::vector<int> domainCoords= vm["domain"].as<std::vector <int> >();
  Z3i::Point ptLower(domainCoords[0],domainCoords[1], domainCoords[2]);
  Z3i::Point ptUpper(domainCoords[3],domainCoords[4], domainCoords[5]);
  Image3D::Domain imageDomain(ptLower, ptUpper);
  
  string inputSDP = vm["sdp"].as<std::string>();
  string outputFilename = vm["output-file"].as<std::string>();
  int foregroundVal = vm["foregroundVal"].as<int>();
  int backgroundVal = vm["backgroundVal"].as<int>();
  

  vector<unsigned int> vPos;
  vPos.push_back(0);
  vPos.push_back(1);
  vPos.push_back(2);
  trace.info() << "Reading input SDP file: " << inputSDP  ; 
  std::vector<Z3i::Point> vectPoints=  PointListReader<Z3i::Point>::getPointsFromFile(inputSDP, vPos); 
  trace.info() << " [done] " << std::endl ; 
  
  Image3D imageResult(imageDomain); 
  for(Image3D::Domain::ConstIterator iter = imageResult.domain().begin(); iter!= imageResult.domain().end();
      iter++){
    imageResult.setValue(*iter, backgroundVal);
  }    
  for(int i=0; i<vectPoints.size(); i++){
    if(vm.count("invertY")){
      vectPoints[i][1]=ptUpper[1]-vectPoints[i][1];
    }
    imageResult.setValue(vectPoints[i], foregroundVal);
  }
  
  GenericWriter<Image3D>::exportFile(outputFilename, imageResult);

  return 0;
  
}




