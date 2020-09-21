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
 * @file volFlip.cpp
 * @ingroup volumetric/voltools
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2013/08/06
 *
 * 
 *
 * This file is part of the DGtalTools.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ConstImageAdapter.h>

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
 @page volFlip volFlip
 
 @brief  Flips 2D slice image of an 3D vol image (mirror transformation).

 @b Usage:  volFlip --input \<volFileName\> --imagePlane 0 1 --flipDimension 0 --o \<volOutputFileName\> (vol, longvol, p3d format)



 @b Allowed @b options @b are : 
 @code
 Positionals:
   1 TEXT:FILE REQUIRED                  Input vol file.

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         Input vol file.
   --imagePlane UINT=[0,1] x 2           arg=  {0,1,2} x {0,1,2} defines the axis of the slice image which will be transformed (by default arg= 0 1  i.e. the slice image defined in the X,Y plane (Z=cst)
   --flipDimension UINT=0                specify which axis will be used to apply the flip.
   -o,--output TEXT=result.vol           Output filename.
   
 
 
 @endcode

 @b Example: 

 @code
 $ volFlip --imagePlane 0 1 --flipDimension 0 -i ${DGtal}/examples/samples/lobster.vol -o flippedXxyLobster.vol 
 @endcode


 You should obtain such a result:
 @image html resVolFlip.png "(a) source image (b) flipped version."
 
 @see
 @ref volFlip.cpp

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
  exit ( 1 );
}


int main(int argc, char**argv)
{
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.vol"};
  
  std::vector <unsigned int> vectImgPlane {0,1};
  unsigned int dimFlip {0};
  
  app.description("Flip 2D slice image of an 3D vol image (mirror transformation)\nBasic usage:\n \t volFlip --input <volFileName> --imagePlane 0 1 --flipDimension 0 --o <volOutputFileName> (vol, longvol, p3d format)\n Example:\n volFlip --imagePlane 0 1 --flipDimension 0 -i ${DGtal}/examples/samples/lobster.vol -o flippedXxyLobster.vol \n The resulting Z slice images (Z= cst) of flippedXxyLobster.p3d will appears flipped according the x axis. \n");
  app.add_option("-i,--input,1", inputFileName, "Input vol file." )
  ->required()
  ->check(CLI::ExistingFile);
  
  app.add_option("--imagePlane", vectImgPlane, "arg=  {0,1,2} x {0,1,2} defines the axis of the slice image which will be transformed (by default arg= 0 1  i.e. the slice image defined in the X,Y plane (Z=cst)", true)
   ->expected(2);
  app.add_option("--flipDimension",dimFlip,"specify which axis will be used to apply the flip.", true);
  app.add_option("-o,--output",outputFileName, "Output filename.", true );

 
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  unsigned int dimFirstImg = vectImgPlane.at(0);
  unsigned int dimSecondImg = vectImgPlane.at(1);
  
  unsigned int normalImgDim = (dimFirstImg!=2 && dimSecondImg!=2)? 2 :( (dimFirstImg!=1 && dimSecondImg!=1)? 1:  0    );  
    
    
  trace.beginBlock("Loading file");
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  Image3D;
  Image3D  imageSRC =  GenericReader<Image3D>::import ( inputFileName );
  trace.endBlock();
  Image3D  imageRes(imageSRC.domain());
  for(int i=0; i <= imageSRC.domain().upperBound()[normalImgDim]; i++){
    Point startPoint(0,0, 0);
    startPoint[normalImgDim]=i;
    for( Domain::ConstSubRange::ConstIterator 
	   it = imageSRC.domain().subRange(dimFirstImg, dimSecondImg, 
					   startPoint).begin(),
	   itend =  imageSRC.domain().subRange(dimFirstImg, dimSecondImg, startPoint).end();
	 it != itend; ++it){
      Point pt = *it;
      pt[dimFlip]= imageSRC.domain().upperBound()[dimFlip] - pt[dimFlip] ;
      imageRes.setValue(*it, imageSRC(pt)); 
    }
  }

  trace.beginBlock("Exporting...");
  bool res =  VolWriter< Image3D>::exportVol(outputFileName, imageRes);
  trace.endBlock();
  if (res) return 0; else return 1;
}
