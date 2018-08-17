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
 * @file heightfield2vol.cpp
 * @ingroup converters
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2015/03/18
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
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

/**
 @page heightfield2vol heightfield2vol
 @brief  Converts a 2D heightfield image into a volumetric file.

@b Usage: heightfield2vol [input] [output]

@b Allowed @b options @b are:

@code
  -h [ --help ]                       display this message
  -i [ --input ] arg                  heightfield file.
  -o [ --output ] arg                 volumetric file 
  -s [ --scale ] arg (=1)             set the scale factor on height values. 
                                      (default 1.0)
  -z [ --volZ ] arg (=255)            set the Z max value of domain.
  -f [ --foregroundValue ] arg (=128) specify the foreground value of the 
                                      resulting voxel.
  -b [ --backgroundValue ] arg (=0)   specify the background value of the 
                                      resulting volumetric file.
@endcode

@b Example:
@code
  $ heightfield2vol -i ${DGtal}/examples/samples/church.pgm -o volResu.vol -s 0.3 -z 50  

@endcode
You will obtain such image:
@image html  resHeightfield2vol.png "Resulting image."
@see heightfield2vol.cpp

*/



// Defining a Helper to get the 3D point functor from an 2DImage 
template<typename TImage2D, typename TPoint3D >
struct Image3DPredicatFrom2DImage{
  typedef  TPoint3D Point3D;
  typedef HyperRectDomain<Z3i::Space>  Domain;
  typedef typename TImage2D::Value Value;
  /**
   *  Construct the predicat given a 2D Image
   **/
  Image3DPredicatFrom2DImage(DGtal::ConstAlias<TImage2D> anImage, double aScale, 
                             unsigned int maxHeight,
                             unsigned int fg, unsigned int bg
                             ):myImageRef(anImage), 
                               myScale(aScale),
                               myMaxHeight(maxHeight),
                               myFG(fg), myBG(bg) {
  }   
  inline
  unsigned int operator()(const Point3D &aPoint)  const {
    functors::Projector<SpaceND<2, typename TImage2D::Integer> > projXY;
    return  (*myImageRef)(projXY(aPoint))*myScale >= aPoint[2] ? myFG: myBG  ;    
  }
  
  inline
  Domain domain() const {
    return Domain(Z3i::Point(0,0,0), Z3i::Point(myImageRef->domain().upperBound()[0],
                                                myImageRef->domain().upperBound()[1],
                                                myMaxHeight) );
  }
  CountedConstPtrOrConstPtr<TImage2D> myImageRef;
  double myScale;  
  unsigned int myMaxHeight;
  unsigned int myFG;
  unsigned int myBG;
};



int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char> Image2D;

  
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "heightfield file." )
    ("output,o", po::value<std::string>(), "volumetric file ") 
    ("scale,s", po::value<double>()->default_value(1.0), "set the scale factor on height values. (default 1.0)")
    ("volZ,z", po::value<unsigned int>()->default_value(255), "set the Z max value of domain.")    
    ("foregroundValue,f", po::value<unsigned int>()->default_value(128), "specify the foreground value of the resulting voxel." )
    ("backgroundValue,b", po::value<unsigned int>()->default_value(0), "specify the background value of the resulting volumetric file.");
  
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
		<< "Convert a 2D heightfield image into a volumetric file. "
		<< general_opt << "\n";
      std::cout << "Example:\n"
		<< "heightfield2vol -i ${DGtal}/examples/samples/church.pgm -o volResu.vol -s 0.3 -z 50  \n";
      return 0;
    }
  
  if(! vm.count("input") ||! vm.count("output"))
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;      
      return 0;
    }
  
  string inputFilename = vm["input"].as<std::string>();
  string outputFilename = vm["output"].as<std::string>();

  trace.info() << "Reading input file " << inputFilename ; 
  Image2D inputImage = DGtal::GenericReader<Image2D>::import(inputFilename);  
  double scale = vm["scale"].as<double>(); 
  unsigned int  maxZ = vm["volZ"].as<unsigned int>();
  trace.info() << " [done] " << std::endl ; 
  
  unsigned int foregroundValue = vm["foregroundValue"].as<unsigned int>();
  unsigned int backgroundValue = vm["backgroundValue"].as<unsigned int>();
  
  typedef Image3DPredicatFrom2DImage<Image2D, Z3i::Point> HeightMapVol;
  Image3DPredicatFrom2DImage<Image2D, Z3i::Point> image3Dpredicate(inputImage, scale, maxZ, foregroundValue, backgroundValue);  
  trace.info() << "Processing image to output file " << outputFilename ; 
    
  VolWriter<HeightMapVol>::exportVol(outputFilename, image3Dpredicate);

  trace.info() << " [done] " << std::endl ;   
  return 0;  
}
