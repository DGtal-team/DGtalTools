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

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;



/**
 @page heightfield2vol heightfield2vol
 @brief  Converts a 2D heightfield image into a volumetric file.

@b Usage: heightfield2vol [OPTIONS] 1 [2]

@b Allowed @b options @b are:

@code
Positionals:
  1 TEXT REQUIRED                       input heightfield file (2D image).
  2 TEXT                                output volumetric file.

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT REQUIRED              input heightfield file (2D image).
  -o,--output TEXT                      output volumetric file.
  -s,--scale FLOAT                      set the scale factor on height values (default 1.0)
  -z,--volZ UINT                        set the Z max value of domain.
  -f,--foregroundValue UINT             specify the foreground value of the resulting voxel.
  -b,--backgroundValue UINT             specify the background value of the resulting volumetric file.
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

// parse command line using CLI ----------------------------------------------
   CLI::App app;
   std::string inputFileName;
   std::string outputFileName {"result.vol"};
    
   unsigned int foregroundValue = {128};
   unsigned int backgroundValue = {0};
   double scale {1.0};
   unsigned int  maxZ {255};
 

   app.description("Convert a 2D heightfield image into a volumetric file.\n Example: \n heightfield2vol -i ${DGtal}/examples/samples/church.pgm -o volResu.vol -s 0.3 -z 50  \n");
   app.add_option("-i,--input,1", inputFileName, "input heightfield file (2D image).")
     ->required();
   app.add_option("-o,--output,2", outputFileName,"output volumetric file.", true);
   app.add_option("-s,--scale", scale, "set the scale factor on height values (default 1.0)");
   app.add_option("-z,--volZ", maxZ, "set the Z max value of domain.");
   app.add_option("-f,--foregroundValue", foregroundValue, "specify the foreground value of the resulting voxel.");
   app.add_option("-b,--backgroundValue", backgroundValue, "specify the background value of the resulting volumetric file.");

   app.get_formatter()->column_width(40);
   CLI11_PARSE(app, argc, argv);
   // END parse command line using CLI ----------------------------------------------
  
   
  trace.info() << "Reading input file " << inputFileName ; 
  Image2D inputImage = DGtal::GenericReader<Image2D>::import(inputFileName);  

  trace.info() << " [done] " << std::endl ; 
  
  
  typedef Image3DPredicatFrom2DImage<Image2D, Z3i::Point> HeightMapVol;
  Image3DPredicatFrom2DImage<Image2D, Z3i::Point> image3Dpredicate(inputImage, scale, maxZ, foregroundValue, backgroundValue);  
  trace.info() << "Processing image to output file " << outputFileName ; 
    
  VolWriter<HeightMapVol>::exportVol(outputFileName, image3Dpredicate);
  trace.info() << " [done] " << std::endl ;   
  return EXIT_SUCCESS;  
}
