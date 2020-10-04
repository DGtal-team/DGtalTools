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

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;


/**
 @page slice2vol slice2vol
 @brief  Converts set of 2D images into volumetric file  (pgm3d, vol, longvol).

@b Usage: slice2vol [input] [output]

@b Allowed @b options @b are:

@code

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT ... REQUIRED          input 2D files (.pgm)
  -o,--output TEXT                      volumetric file (.vol, .longvol .pgm3d)
  -s,--sliceOrientation UINT:{0,1,2}=2  specify the slice orientation for which the slice are defined (by default =2 (Z direction))
@endcode

@b Example:
@code
  $ slice2vol -i slice1.pgm slice2.pgm slice3.pgm  -o out.vol
@endcode

@see
slice2vol.cpp

*/



int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char > Image2D;
  
  // parse command line using CLI ----------------------------------------------
   CLI::App app;
   std::vector<std::string> vectImage2DNames;
   std::string outputFileName {"result.vol"};
   unsigned int sliceOrientation {2};

   app.description("Converts set of 2D images into volumetric file  (pgm3d, vol, longvol).\nExample:\n slice2vol -i slice1.pgm slice2.pgm slice3.pgm  -o out.vol see vol2slice");
   app.add_option("-i,--input", vectImage2DNames, "input 2D files (.pgm)")
     -> required();
   app.add_option("-o,--output", outputFileName, "volumetric file (.vol, .longvol .pgm3d)");
   app.add_option("--sliceOrientation,-s", sliceOrientation, "specify the slice orientation for which the slice are defined (by default =2 (Z direction))", true)
     -> check(CLI::IsMember({0, 1, 2}));
   
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------
  


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
  return EXIT_SUCCESS;

}




