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
 * @file
 * @ingroup volumetric
 * @author Bertrand Kerautret (\c bertrand.kerautret@univ-lyon2.fr )
 * @author Jonas Lamy (\c jonas.lamy@univ-lyon2.fr )
 *
 *
 * @date 2019/03/01
 *
 * Source file of the tool volMask
 *
 * This file is part of the DGtal library/DGtalTools Project.
 */

///////////////////////////////////////////////////////////////////////////////
#include "DGtal/base/Common.h"
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/GenericWriter.h>

#include "CLI11.hpp"

#ifdef WITH_ITK
#include "DGtal/io/readers/ITKReader.h"
#include "DGtal/io/writers/ITKWriter.h"

#endif

///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace DGtal;
///////////////////////////////////////////////////////////////////////////////

/**
 @page volMask volMask
 
 @brief  Extracts a new image from the a mask image that represents the regions of the image which are selected and copied in the resulting image. Elements outside the regions defined by the mask are set to 0.
 
 @b Usage:   volMask [input]
 
 @b Allowed @b options @b are :
 
 @code
 Positionals:
   1 TEXT:FILE REQUIRED                  an input 3D image vol (or ITK: .nii, mha, ... ) file.
   2 TEXT=result.vol                     the output masked image.

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         an input 3D image vol (or ITK: .nii, mha, ... ) file.
   -t,--inputType TEXT                   to specify the input image type (int or double).
   -a,--mask TEXT:FILE                   the mask image that represents the elements that are copied as output in the resulting image (by default set to 1 you can change this value by using --maskValue).
   -o,--output TEXT=result.vol           the output masked image.
   -f,--offsetBorder UINT=0              add a border offset to the bounding box of the masked value domain.
   -m,--maskValue INT=1                  the masking value.
 @endcode
 
 @b Example:
 
 @code
 volMask -i ${DGtal}/examples/samples/lobster.vol  -a ${DGtal}/examples/samples/lobster.vol -o lobsMasked.vol -m 100
 @endcode
 
 @image html resvolMask.png "Example of result. "
 
 @see
 @ref volMask.cpp
 
 */




typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
#ifdef WITH_ITK
typedef ImageContainerBySTLVector < Z3i::Domain,  double > Image3D_D;
typedef ImageContainerBySTLVector < Z3i::Domain,  int > Image3D_I;
#endif

/**
 * Computes the minimal subdomain containing a masked value.
 * If no masked value appears it returns the image domain.
 **/
template<typename TImage, typename TImageMask>
typename TImageMask::Domain
subDomainMasked(const TImage &image, const  TImageMask &maskImage,
                typename TImageMask::Value maskValue, unsigned int domainOffset=1){
  typename TImageMask::Domain res;
  Z3i::Point minP = image.domain().upperBound();
  Z3i::Point maxP = image.domain().lowerBound();
  
  Z3i::Point::Iterator minIt;
  Z3i::Point::Iterator maxIt;
  bool foundMaskedVal = false;
  for(const auto &p: image.domain())
  {
    minIt = minP.begin();
    maxIt = maxP.begin();
    if( maskImage.domain().isInside(p) && maskImage(p) ==  maskValue ) // no noise on mask image
    {
      foundMaskedVal=true;
      for(auto pIt=p.begin(); pIt!=p.end();pIt++ )
      {
        if( *pIt < *minIt ){*minIt = *pIt;}
        if( *pIt > *maxIt ){*maxIt = *pIt;}
        minIt++;
        maxIt++;
      }
    }
  }
  if (!foundMaskedVal)
  {
    trace.info() << "No masked value found resulting image will be empty." << std::endl;
    return image.domain();

  }
  // offset to avoid problems on borders
  Z3i::Point offset(domainOffset,domainOffset,domainOffset);
  minP -= offset;
  maxP += offset;
  trace.info() << "sub-domain:" << minP << " " << maxP << std::endl;
  return typename TImageMask::Domain(minP, maxP);
}


template<typename TImage, typename TImageMask>
void
applyMask(const TImage &inputImage,TImage &outputImage,
          const TImageMask &maskImage, typename TImageMask::Value maskValue)
{
  for (const auto &p: outputImage.domain())
  {
    if (inputImage.domain().isInside(p) &&  maskImage(p) ==  maskValue)
    {
      outputImage.setValue(p, inputImage(p) );
    }
  }
  
}

template<typename TImage, typename TImageMask>
void
processImage(const TImage &inputImage, const TImageMask &maskImage,
             typename TImageMask::Value maskValue, std::string outputFileName, unsigned int offsetBorder ){
  // First step getting the bounding box of the domain:
  auto subDm = subDomainMasked(inputImage, maskImage, maskValue, offsetBorder);
  TImage outputImage( subDm );
  // Second step: masking source image
  applyMask(inputImage, outputImage, maskImage,  maskValue);
  trace.info() << "writing output image...";
  GenericWriter<TImage>::exportFile(outputFileName, outputImage);
}


int main( int argc, char** argv )
{

  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.vol"};
  std::string inputType;
  std::string maskFileName;
  unsigned int offsetBorder {0};
  int maskValue {1};
  
  app.description("Outputs a new image from two input images, one representing the data, one representing the selection mask. The size of output image is the size of the bounding box of selected values, plus the chosen border offset. \n Typical use example:\n \t volMask -i ${DGtal}/examples/samples/lobster.vol  -a ${DGtal}/examples/samples/lobster.vol -o lobsMasked.vol -m 100  \n");
  
#ifdef WITH_ITK
  app.add_option("-i,--input,1", inputFileName, "an input 3D image vol (or ITK: .nii, mha, ... ) file." )
      ->required()
      ->check(CLI::ExistingFile);
  app.add_option("--inputType,-t",inputType, "to specify the input image type (int or double).");

  #else
  app.add_option("-i,--input,1", inputFileName, "an input vol file." )
      ->required()
      ->check(CLI::ExistingFile);
#endif
  
  app.add_option("--mask,-a",maskFileName, "the mask image that represents the elements that are copied as output in the resulting image (by default set to 1 you can change this value by using --maskValue). ")
  ->check(CLI::ExistingFile);
  app.add_option("-o,--output,2", outputFileName, "the output masked image.", true );
  app.add_option("--offsetBorder,-f", offsetBorder, "add a border offset to the bounding box of the masked value domain.", true);
  app.add_option("--maskValue,-m", maskValue, "the masking value.", true);
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  
  trace.info() << "Reading mask image...";
  Image3D maskImage = DGtal::GenericReader<Image3D>::import(maskFileName);
  trace.info() << "[done]"<< std::endl;
  trace.info() << "Reading input image...";
#ifdef WITH_ITK
  if (inputType=="double")
  {
    Image3D_D inputImage = DGtal::GenericReader<Image3D_D>::import(inputFileName);
    trace.info() << "[done]"<< std::endl;
    processImage(inputImage, maskImage, maskValue, outputFileName, offsetBorder);
  }
  else if (inputType=="int")
  {
    Image3D_I inputImage = DGtal::GenericReader<Image3D_I>::import(inputFileName);
    trace.info() << "[done]"<< std::endl;
    processImage(inputImage, maskImage, maskValue, outputFileName, offsetBorder);
  }
  else
  {
    Image3D inputImage = DGtal::GenericReader<Image3D>::import(inputFileName);
    trace.info() << "[done]"<< std::endl;
    processImage(inputImage, maskImage, maskValue, outputFileName, offsetBorder);
  }
#else
  Image3D inputImage = DGtal::GenericReader<Image3D>::import(inputFileName);
  processImage(inputImage, maskImage, maskValue, outputFileName, offsetBorder);
#endif
  
  trace.info() << "[Done]" << std::endl;
  return 0;
}

