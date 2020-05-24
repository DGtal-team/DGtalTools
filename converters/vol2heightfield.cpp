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
 * @file vol2heightfield.cpp
 * @ingroup Converters
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
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;



/**
 @page vol2heightfield vol2heightfield
 @brief  Converts volumetric  file into a projected 2D image given from a normal direction N and from a starting point P.

 The 3D volume is scanned in this normal direction N starting from P with a step 1. If the intensity of the 3d point is inside the given thresholds its 2D gray values are set to the current scan number.


@b Usage: vol2heightfield [input] [output]

@b Allowed @b options @b are:

@code

Positionals:
  1 TEXT:FILE REQUIRED                  vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
  -o,--output TEXT                      sequence of discrete point file (.sdp)
  -m,--thresholdMin INT=128             threshold min (excluded) to define binary shape.
  -M,--thresholdMax INT=255             threshold max (included) to define binary shape.
  --rescaleInputMin INT=0               min value used to rescale the input intensity (to avoid basic cast into 8  bits image).
  --rescaleInputMax INT=255             max value used to rescale the input intensity (to avoid basic cast into 8  bits image).
  --nx FLOAT=0                          set the x component of the projection direction.
  --ny FLOAT=0                          set the y component of the projection direction.
  --nz FLOAT=1                          set the z component of the projection direction.
  -x,--centerX UINT=0                   choose x center of the projected image.
  -y,--centerY UINT=0                   choose y center of the projected image.
  -z,--centerZ UINT=0                   choose z center of the projected image.
  -w,--width UINT=100                   set the width of the resulting height Field image.
  --height UINT=100                     set the height of the resulting height Field image.
  --heightFieldMaxScan UINT             set the maximal scan deep.
  --setBackgroundLastDepth              change the default background (black with the last filled intensity).

@endcode

@b Example:
@code 
$ vol2heightfield -i ${DGtal}/examples/samples/lobster.vol -m 60 -M 500  --nx 0 --ny 0.7 --nz -1 -x 150 -y 0 -z 150 --width 300 --height 300 --heightFieldMaxScan 350  -o resultingHeightMap.pgm 
@endcode

You should obtain such a resulting image:
@image html resVol2heightfield.png "resulting image."
@see
@ref vol2heightfield.cpp

*/


int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char> Image2D;
  typedef DGtal::ConstImageAdapter<Image3D, Z2i::Domain, DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>,
                                   Image3D::Value,  DGtal::functors::Identity >  ImageAdapterExtractor;

// parse command line using CLI ----------------------------------------------
   CLI::App app;
   std::string inputFileName;
   std::string outputFileName {"result.raw"};
   int thresholdMin {128};
   int thresholdMax {255};
   DGtal::int64_t rescaleInputMin {0};
   DGtal::int64_t rescaleInputMax {255};
   double nx {0};
   double ny {0};
   double nz {1};
   unsigned int centerX {0};
   unsigned int centerY {0};
   unsigned int centerZ {0};
   unsigned int heightImageScan {100};
   unsigned int widthImageScan {100};
   unsigned int heightFieldMaxScan {255};
   unsigned int maxScan;
   bool bgLastDepth = false;
   
   app.description("Convert volumetric  file into a projected 2D image given from a normal direction N and from a starting point P. The 3D volume is scanned in this normal direction N starting from P with a step 1. If the intensity of the 3d point is inside the given thresholds its 2D gray values are set to the current scan number.\n  Example:\n vol2heightfield -i ${DGtal}/examples/samples/lobster.vol -m 60 -M 500  --nx 0 --ny 0.7 --nz -1 -x 150 -y 0 -z 150 --width 300 --height 300 --heightFieldMaxScan 350  -o resultingHeightMap.pgm");
   app.add_option("-i,--input,1", inputFileName, "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("-o,--output", outputFileName, "sequence of discrete point file (.sdp)");
  app.add_option("--thresholdMin,-m", thresholdMin, "threshold min (excluded) to define binary shape.", true);
  app.add_option("--thresholdMax,-M", thresholdMax, "threshold max (included) to define binary shape.", true);
  app.add_option("--rescaleInputMin", rescaleInputMin, "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);
  app.add_option("--rescaleInputMax", rescaleInputMax, "max value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);
 
  app.add_option("--nx",nx, "set the x component of the projection direction.", true);
  app.add_option("--ny",ny, "set the y component of the projection direction.", true);
  app.add_option("--nz",nz, "set the z component of the projection direction.", true);
  app.add_option("--centerX,-x", centerX, "choose x center of the projected image.", true);
  app.add_option("--centerY,-y", centerY, "choose y center of the projected image.", true);
  app.add_option("--centerZ,-z", centerZ, "choose z center of the projected image.", true);
  app.add_option("--width,-w", widthImageScan, "set the width of the resulting height Field image.", true);
  app.add_option("--height", heightImageScan, "set the height of the resulting height Field image.", true);
  app.add_option("--heightFieldMaxScan",maxScan, "set the maximal scan deep.");
  app.add_flag("--setBackgroundLastDepth", bgLastDepth,"change the default background (black with the last filled intensity).");
  
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------


  trace.info() << "Reading input file " << inputFileName ; 

  typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
  Image3D inputImage =  GenericReader< Image3D >::importWithValueFunctor( inputFileName,RescalFCT(rescaleInputMin,
                                                                                                  rescaleInputMax,
                                                                                                  0, 255) );

  trace.info() << " [done] " << std::endl ;   
  std::ofstream outStream;
  outStream.open(outputFileName.c_str());
  
  trace.info() << "Processing image to output file " << outputFileName << std::endl;   

  if(maxScan > std::numeric_limits<Image2D::Value>::max()){
    maxScan = std::numeric_limits<Image2D::Value>::max();
    trace.warning()<< "value --setBackgroundLastDepth outside mox value of image. Set to max value:" << maxScan << std::endl; 
  }
  

  Image2D::Domain aDomain2D(DGtal::Z2i::Point(0,0), 
                          DGtal::Z2i::Point(widthImageScan, heightImageScan));
  Z3i::Point ptCenter (centerX, centerY, centerZ);
  Z3i::RealPoint normalDir (nx, ny, nz);
  Image2D resultingImage(aDomain2D);
  
  for(Image2D::Domain::ConstIterator it = resultingImage.domain().begin(); 
      it != resultingImage.domain().end(); it++){
    resultingImage.setValue(*it, 0);
  }
  DGtal::functors::Identity idV;
  
  unsigned int maxDepthFound = 0;
  for(unsigned int k=0; k < maxScan; k++){
    Z3i::Point c (ptCenter+normalDir*k, DGtal::functors::Round<>());
    DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedder(inputImage.domain(), 
                                                                        c,
                                                                        normalDir,
                                                                        widthImageScan);
    ImageAdapterExtractor extractedImage(inputImage, aDomain2D, embedder, idV);
    for(Image2D::Domain::ConstIterator it = extractedImage.domain().begin(); 
        it != extractedImage.domain().end(); it++){
      if( resultingImage(*it)== 0 &&  extractedImage(*it) < thresholdMax &&
          extractedImage(*it) > thresholdMin){
        maxDepthFound = k;
        resultingImage.setValue(*it, maxScan-k);
      }
    }    
  }
  if (bgLastDepth) {
    for(Image2D::Domain::ConstIterator it = resultingImage.domain().begin(); 
        it != resultingImage.domain().end(); it++){
      if( resultingImage(*it)== 0 ){
        resultingImage.setValue(*it, maxScan-maxDepthFound);
      }
    }
  }   
  resultingImage >> outputFileName;  
  trace.info() << " [done] " << std::endl ;   
  return EXIT_SUCCESS;
}




