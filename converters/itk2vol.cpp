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
 * @file itk2vol.cpp
 * @ingroup converters
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2015/06/06
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
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/readers/ITKReader.h"

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;

/**
 @page itk2vol itk2vol
 @brief  Converts itk file into a volumetric file (.vol, .pgm3d).

@b Usage: itk2vol [input] [output]

@b Allowed @b options @b are:

@code
Positionals:
  1 TEXT:FILE REQUIRED                  Any file format in the ITK library (mhd, mha, ...).
  2 TEXT                                volumetric file (.vol, .pgm3d).

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         Any file format in the ITK library (mhd, mha, ...).
  -o,--output TEXT                      volumetric file (.vol, .pgm3d).
  -m,--maskImage TEXT                   Use a mask image to remove image part (where mask is 0). The default mask value can be changed using mask default value.
  -r,--maskRemoveLabel INT              Change the label value that defines the part of input image to be removed by the option --maskImage.
  --inputMin INT                        set minimum density threshold on Hounsfield scale.
  --inputMax INT                        set maximum density threshold on Hounsfield scale.
  -t,--inputType TEXT:{int,double}      to sepcify the input image type (int or double).

@b Example:
@code
$itk2vol -i image.mhd --dicomMin -500 --dicomMax -100 -o sample.vol 
@endcode

@see itk2vol.cpp

*/



template<typename TImage, typename TImageMask>
void
applyMaskImage( TImage &imageInput,  const  TImageMask &maskImage,
                typename TImageMask::Value valRemove)
{
  for(const auto &p : imageInput.domain())
  {
    if (maskImage(p) == valRemove)
    {
      imageInput.setValue(p,0);
    }
  }
  
}




int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3DChar;
  typedef ImageContainerBySTLVector < Z3i::Domain,  double > Image3D_D;
  typedef ImageContainerBySTLVector < Z3i::Domain,  int > Image3D_I;


   // parse command line using CLI ----------------------------------------------
   CLI::App app;
   std::string inputFileName;
   std::string outputFileName {"result.vol"};
   DGtal::int64_t inputMin {-1000};
   DGtal::int64_t inputMax  {3000};
   int  maskRemoveLabel {0};  
   string inputMask {""};
   string inputType {""};
   
   app.description("Converts itk file into a volumetric file (.vol, .pgm3d). \n Example:\n itk2vol -i image.mhd --inputMin -500 --inputMax -100 -o sample.vol \n");
   app.add_option("-i,--input,1", inputFileName, "Any file format in the ITK library (mhd, mha, ...)." )
    ->required()
    ->check(CLI::ExistingFile);
   app.add_option("-o,--output,2", outputFileName, "volumetric file (.vol, .pgm3d).");
   app.add_option("-m,--maskImage", inputMask, "Use a mask image to remove image part (where mask is 0). The default mask value can be changed using mask default value.");
   app.add_option("-r,--maskRemoveLabel", maskRemoveLabel,"Change the label value that defines the part of input image to be removed by the option --maskImage." );
   app.add_option("--inputMin", inputMin, "set minimum density threshold on Hounsfield scale.");
   app.add_option("--inputMax", inputMax, "set maximum density threshold on Hounsfield scale.");
   app.add_option("-t,--inputType", inputType, "to sepcify the input image type (int or double).")
     -> check(CLI::IsMember({"int", "double"}));
   
   app.get_formatter()->column_width(40);
   CLI11_PARSE(app, argc, argv);
   // END parse command line using CLI ----------------------------------------------

   
  if (inputType == "double") {
      typedef DGtal::functors::Rescaling<double ,unsigned char > RescalFCT;
      trace.info() << "Reading input input file " << inputFileName ; 
      Image3D_D inputImage = ITKReader< Image3D_D  >::importITK(inputFileName);
      trace.info() << " [done] " << std::endl ; 
      trace.info() << " converting into vol file... " ; 
      if ( inputMask != "")
      {
        Image3D_I maskImage = ITKReader< Image3D_I  >::importITK(inputMask);
        applyMaskImage(inputImage, maskImage, maskRemoveLabel);
      }

      RescalFCT rescaleCustom(inputMin, inputMax, 0, 255);
      DGtal::GenericWriter<Image3D_D, 3, unsigned char, RescalFCT>::exportFile(outputFileName,
                                                                               inputImage,
                                                                               "UInt8Array3D",
                                                                               rescaleCustom);
  }else {
     typedef DGtal::functors::Rescaling<int ,unsigned char > RescalFCT;
      trace.info() << "Reading input input file " << inputFileName ; 
      Image3D_I inputImage = ITKReader< Image3D_I  >::importITK(inputFileName);
      trace.info() << " [done] " << std::endl ; 
      trace.info() << " converting into vol file... " ; 
      RescalFCT rescaleCustom(inputMin, inputMax, 0, 255);
      if (inputMask != "")
      {
        Image3D_I maskImage = ITKReader< Image3D_I  >::importITK(inputMask);
        applyMaskImage(inputImage, maskImage, maskRemoveLabel);
      }
      
      DGtal::GenericWriter<Image3D_I, 3, unsigned char, RescalFCT>::exportFile(outputFileName,
                                                                               inputImage,
                                                                               "UInt8Array3D",
                                                                               rescaleCustom);
  }
  trace.info() << " [done] " << std::endl ;   
  return EXIT_SUCCESS;  
}
