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
 * @file volImageMetrics.cpp
 *
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2013/07/20
 *
 * This file is part of the DGtal library.
 */


/**
 @page volImageMetrics volImageMetrics
 
 @brief  Applies basic image measures (RMSE, PSNR) between two volumetric images A and B.

 @b Usage:  volImageMetrics --volA \<volAFilename\> --volB \<volBFilename\> 


 @b Allowed @b options @b are : 
 @code
 
  Positionals:
    1 TEXT:FILE REQUIRED                  Input filename of volume A (vol format, and other pgm3d can also be used).
    2 TEXT:FILE REQUIRED                  Input filename of volume B (vol format, and other pgm3d can also be used).

  Options:
    -h,--help                             Print this help message and exit
    -a,--volA TEXT:FILE REQUIRED          Input filename of volume A (vol format, and other pgm3d can also be used).
    -b,--volB TEXT:FILE REQUIRED          Input filename of volume B (vol format, and other pgm3d can also be used).

 @endcode

 @b Example: 

 @code
 # generating another input vol file using tutorial example (eroded.vol):
 $DGtal/build/examples/tutorial-examples/FMMErosion
 # compare the two images:
 $ volImageMetrics -a eroded.vol -b $DGtal/examples/samples/cat10.vol 
 @endcode


 You should obtain such an output:
@verbatim
# Image based measures (generated with volImageMetrics) given with the image A: eroded.voland the image B: /Users/kerautre/EnCours/DGtal/examples/samples/cat10.vol
#  RMSE PSNR 
 33.9411 171.331
@endverbatim
 
 @see  \ref volImageMetrics.cpp

 */

#include <iostream>
#include <limits>

#include <DGtal/base/Common.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/math/Statistic.h>

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;
using namespace Z3i;

typedef ImageContainerBySTLVector < Z3i::Domain,  int > Image3D;
typedef ImageContainerBySTLVector < Z2i::Domain,  int > Image2D;

double
getRMSE(const Image3D & imageA, const Image3D &imageB){
  double sumDiff=0;
  for(Image3D::Domain::ConstIterator it = imageA.domain().begin(); it!=imageA.domain().end(); it++){
    sumDiff+=(imageA(*it)-imageB(*it))*(imageA(*it)-imageB(*it));
  }
  return sqrt(sumDiff/imageA.domain().size());
}

double
getPSNR(const Image3D & imageA, const Image3D &imageB, double rmsd){
  unsigned long long int d =  std::numeric_limits<Image3D::Value>::max();
  return 10.0*log10(d*d/rmsd);
}


int main(int argc, char**argv)
{
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileNameVolA;
  std::string inputFileNameVolB;

  app.description("Apply basic image measures (RMSE, PSNR) between two volumetric images A and B. \n Basic usage:\n \t volImageMetrics --volA <volAFilename> --volB <volBFilename> \n Typical use :\n  volImageMetrics -a imageA.vol  -b imageB.vol \n");
  app.add_option("-a,--volA,1", inputFileNameVolA, "Input filename of volume A (vol format, and other pgm3d can also be used)." )
  ->required()
  ->check(CLI::ExistingFile);
  app.add_option("-b,--volB,2", inputFileNameVolB, "Input filename of volume B (vol format, and other pgm3d can also be used)." )
  ->required()
  ->check(CLI::ExistingFile);
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  Image3D imageA = GenericReader<Image3D>::import(inputFileNameVolA);
  Image3D imageB = GenericReader<Image3D>::import(inputFileNameVolB);
   
  std::cout << "# Image based measures (generated with volImageMetrics) given with the image A: "<< inputFileNameVolA<< " and the image B: "<< inputFileNameVolB << endl;
  std::cout << "#  RMSE PSNR "<< endl;    
  
  double rmse= getRMSE(imageA, imageB);
  double psnr= getPSNR(imageA, imageB, rmse);
    
  std::cout << " " << rmse << " " << psnr << endl;
  return 1;
}
