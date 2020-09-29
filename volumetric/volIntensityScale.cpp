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
 * @file volIntensityScale.cpp
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2015/01/12
 *
 * 
 *
 * This file is part of the DGtalTools.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/VolReader.h>
#include <DGtal/io/writers/VolWriter.h>
#include "DGtal/io/writers/GenericWriter.h"
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ConstImageAdapter.h>

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
 @page volIntensityScale volIntensityScale
 
 @brief Applies a linear rescaling of the image intensity from an input intensity interval [InMin, InMax] into an output interval [OutMin, OutMax].

 @b Usage: 	  ./volumetric/volIntensityScale [OPTIONS] 1 (image files can be independently in vol, pgm3D, p3d format)


 @b Allowed @b options @b are : 
 @code

 Positionals:
   1 TEXT:FILE REQUIRED                  Input vol file.

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         Input vol file.
   -o,--output TEXT=result.vol           volumetric output file (.vol, .pgm, .pgm3d, .longvol)
   -m,--inMin INT=0                      the min value of the input image.
   -M,--inMax INT=255                    the max value of the input image.
   --outMin INT=0                        the min value of the output image.
   --outMax INT=255                      the max value of the output image.
   
 @endcode

 @b Example: 

 In this example, we select all intensities included in [0,100] and scale them into [0, 255]: 

 @code
  $ volIntensityScale -i  ${DGtal}/examples/samples/lobster.vol  --inMin 0 --inMax 100 -o lobster0-100.vol 

 @endcode


 You should obtain such a result:
 @image html resVolIntensityScale.png "(a) source image, (b) intensity scaled image."
 
 @see
 @ref volIntensityScale.cpp

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
}


int main(int argc, char**argv)
{

  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.vol"};
  int inMin {0};
  int inMax {255};
  int outMin {0};
  int outMax {255};
  
  
  app.description("Apply a linear rescaling of the image intensity from an input intensity interval [InMin, InMax] into an output interval [OutMin, OutMax].\n Basic usage:\n volIntensityScale --input <volFileName> --output <volOutputFileName>  (both files can be independently in vol, pgm3D, p3d format)\n Example: \n volIntensityScale -i  ${DGtal}/examples/samples/lobster.vol  --inMin 0 --inMax 100 -o lobster0-100.vol");
  
  app.add_option("-i,--input,1", inputFileName, "Input vol file." )
  ->required()
  ->check(CLI::ExistingFile);
  
  app.add_option("-o,--output",outputFileName, "volumetric output file (.vol, .pgm, .pgm3d, .longvol) ", true);
  app.add_option("-m,--inMin", inMin,  "the min value of the input image.", true);
  app.add_option("-M,--inMax", inMax,  "the max value of the input image.", true);
  app.add_option("--outMin", outMin,  "the min value of the output image.", true);
  app.add_option("--outMax", outMax,  "the max value of the output image.", true);
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  trace.beginBlock("Loading file");

  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;
  typedef DGtal::functors::Rescaling<int ,unsigned char > RescalFCT;
  MyImageC image = VolReader< MyImageC,  RescalFCT  >::importVol( inputFileName,
                                                                  RescalFCT(inMin,
                                                                            inMax,
                                                                            outMin, outMax));  
  trace.endBlock();

  trace.beginBlock("Exporting...");
  bool res =  GenericWriter<MyImageC>::exportFile(outputFileName, image);
  trace.endBlock();
  if (res) return 0; else return 1;   
}
