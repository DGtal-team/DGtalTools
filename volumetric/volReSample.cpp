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
 * @file volResSample.cpp
 * @ingroup volumetric
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2014/07/15
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
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"

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
  typedef DGtal::functors::BasicDomainSubSampler< HyperRectDomain<SpaceND<3, int> >,  
                                                  DGtal::int32_t, double >   ReSampler; 
  typedef DGtal::ConstImageAdapter<Image3D, Image3D::Domain, ReSampler,
				   Image3D::Value,  DGtal::functors::Identity >  SamplerImageAdapter;


  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string >(), "input volumetric file (.vol, .longvol, .pgm3d) " )
    ("output,o", po::value<std::string>(), "the new volumetric file (.vol, .longvol, .pgm3d) " )
    ("gridSize,g", po::value<std::vector<double> >()->multitoken(), "size_x size_y size_z : the grid size of the re sampling ");

  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);


  if(! vm.count("input")||! vm.count("output"))
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;
      parseOK=false;
    }

  if( !parseOK || vm.count("help") || !vm.count("gridSize") )
    {
      std::cout << "Usage: " << argv[0] << " [input-files] [output-file]\n"
		<< "Re sample a 3D volumetric image (.vol, .longvol, .pgm3d)  with a given grid size."
		<< general_opt << "\n";
      std::cout << "Example: to re sampling an image with scale x,y,z  0.98, 0.98, 5.0 you can obtain it by: \n"
		<< "volResSample -i image3d.vol -g 1.02 1.02 0.2  -o imageReSampled.vol \n" << endl;
      return 0;
    }




  std::vector<  double > aGridSizeReSample = vm["gridSize"].as<std::vector<double > >();
  if(aGridSizeReSample.size()!=3){
    trace.error() << "The grid size should contains 3 elements" << std::endl;
    return 0;
  }

  std::string inputFileName = vm["input"].as<std::string>();
  std::string outputFileName = vm["output"].as<std::string>();
 
 
  trace.info()<< "Importing volume file :  " << inputFileName<< " ... " ;
  Image3D input3dImage = GenericReader<Image3D>::import(inputFileName);
  trace.info()<< "[done]" << endl;
  
  
  PointVector<3,int> shiftVector3D(0 ,0, 0);      
  DGtal::functors::BasicDomainSubSampler< HyperRectDomain<SpaceND<3, int> >,  
                                          DGtal::int32_t, double > reSampler(input3dImage.domain(),
                                                                             aGridSizeReSample,  shiftVector3D);  
  SamplerImageAdapter sampledImage (input3dImage,reSampler.getSubSampledDomain(), reSampler, functors::Identity());
  GenericWriter<SamplerImageAdapter>::exportFile(outputFileName, sampledImage);
  

  return 0;
}
