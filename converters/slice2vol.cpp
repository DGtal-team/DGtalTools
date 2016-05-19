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
 * @file slice2vol.cpp
 * @ingroup surfaceTools
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
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>



using namespace std;
using namespace DGtal;


/**
 @page slice2vol slice2vol
 @brief  Converts set of 2D images into volumetric file  (pgm3d, vol, longvol).

@b Usage: slice2vol [input] [output]

@b Allowed @b options @b are:

@code
  -h [ --help ]                      display this message
  -s [ --sliceOrientation ] arg (=2) specify the slice orientation for which the slice are considered (by default =2 (Z direction))
  -i [ --input ] arg                 input 2D files (.pgm) 
  -o [ --output ] arg                volumetric file (.vol, .longvol .pgm3d) 
@endcode

@b Example:
@code
  $ slice2vol -i slice1.pgm slice2.pgm slice3.pgm  -o out.vol
@endcode

@see
slice2vol.cpp

*/


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char > Image2D;
  
  
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("sliceOrientation,s", po::value<unsigned int>()->default_value(2), "specify the slice orientation for which the slice are considered (by default =2 (Z direction))" )
    ("input,i", po::value<std::vector <std::string> >()->multitoken(), "input 2D files (.pgm) " )
    ("output,o", po::value<std::string>(), "volumetric file (.vol, .longvol .pgm3d) " );
  
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);    
  

  if( ! vm.count("input") || ! vm.count("output") || !parseOK || vm.count("help"))
    {
      std::cout << "Usage: " << argv[0] << " [input-files] [output]\n"
		<< "Converts set of 2D images into volumetric file  (pgm3d, vol, longvol). "
		<< general_opt << "\n";
      std::cout << "Example:\n"
                << "slice2vol -i slice1.pgm slice2.pgm slice3.pgm  -o out.vol \n"
                << "see vol2slice"<<endl; 
       
      return 0;
    }
  
  if(! (vm.count("input") && vm.count("output")) )
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;      
      return 0;
    }



  std::string outputFileName = vm["output"].as<std::string>();
  std::vector<string> vectImage2DNames = vm["input"].as<std::vector<std::string> >();
  unsigned int sliceOrientation = vm["sliceOrientation"].as<unsigned int>();
  std::vector<Image2D> vectImages2D; 
  // Reading all images
  for(unsigned int i=0; i< vectImage2DNames.size(); i++){
    trace.info() << "Reading image " << i ;
    Image2D image = GenericReader<Image2D>::import(vectImage2DNames.at(i));
    vectImages2D.push_back(image);
    trace.info() << " [done]" << std::endl;
  }
  
  Image2D::Domain domImage2D =  vectImages2D.at(0).domain();
  DGtal::functors::Projector<DGtal::Z3i::Space> projIn3Dlower(0); 
  DGtal::functors::Projector<DGtal::Z3i::Space> projIn3Dupper(vectImages2D.size()-1); 
  projIn3Dlower.initAddOneDim(sliceOrientation);
  projIn3Dupper.initAddOneDim(sliceOrientation);
  Image3D::Domain domImage3D (projIn3Dlower(vectImages2D.at(0).domain().lowerBound()),
  			      projIn3Dupper(vectImages2D.at(0).domain().upperBound()));

  Image3D imageResult (domImage3D);
  for( unsigned int i=0; i<vectImages2D.size();  i++){
    Image2D sliceImage = vectImages2D.at(i);
    DGtal::functors::Projector<DGtal::Z3i::Space> projIn3D(i); 
    projIn3D.initAddOneDim(sliceOrientation);
    for(Image2D::Domain::ConstIterator it = sliceImage.domain().begin();  
         it!= sliceImage.domain().end(); it++){
         Z3i::Point pt =projIn3D(*it);
         imageResult.setValue(pt, sliceImage(*it));
    }
  }
  trace.info() << "Exporting 3d image ... " << std::endl ;
  GenericWriter<Image3D>::exportFile(outputFileName, imageResult);
  trace.info()  << "[done]";
  return 0;  
}




