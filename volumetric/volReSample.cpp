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
 * @file volReSample.cpp
 * @ingroup volumetric
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2014/07/15
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
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;


/**
 @page volReSample volReSample
 
 @brief Re samples a 3D volumetric image (.vol, .longvol, .pgm3d)  with a given grid size.
 @ingroup volumetrictools
 @b Usage: v ./volumetric/volReSample [OPTIONS] 1 [2]


 @b Allowed @b options @b are : 
 @code
  Positionals:
    1 TEXT:FILE REQUIRED                  input volumetric file (.vol, .longvol, .pgm3d).
    2 TEXT                                the new volumetric file (.vol, .longvol, .pgm3d).

  Options:
    -h,--help                             Print this help message and exit
    -i,--input TEXT:FILE REQUIRED         input volumetric file (.vol, .longvol, .pgm3d).
    -o,--output TEXT                      the new volumetric file (.vol, .longvol, .pgm3d).
    -g,--gridSize FLOAT x 3               size_x size_y size_z : the grid size of the re sampling
 @endcode

 @b Example: 

Here is an example of re sampling with different grid sizes 2, 4 and 8:

 @code
 $ volReSample  $DGtal/examples/samples/Al.100.vol AlRS2.vol -g 2 2 2
 $ volReSample  $DGtal/examples/samples/Al.100.vol AlRS4.vol -g 4 4 4
 $ volReSample  $DGtal/examples/samples/Al.100.vol AlRS8.vol -g 8 8 8
 @endcode

We can convert the resulting volumetric files into a sequence of discrete points with the tool  @ref vol2sdp :
@code 
$ vol2sdp $DGtal/examples/samples/Al.100.vol -m 1 -o AlRS1.sdp
$ vol2sdp AlRS2.vol -m 1 -o AlRS2.sdp
$ vol2sdp AlRS4.vol -m 1 -o AlRS4.sdp
$ vol2sdp AlRS8.vol -m 1 -o AlRS8.sdp
$ cat AlRS{1,2,4,8}.sdp >> AlRS1_2_4_8.sdp
# display the resulting file:
3dSDPViewer AlRS1_2_4_8.sdp    
@endcode

 Note that if DGtal is compiled with the  option DGTAL_WITH_ITK set to ON, you can export the image in format ITK format and integrating image spacing.

 You should obtain such a result:
 @image html resVolReSample.png "Resulting of re sampling with grid size = 2, 4 and 8."
 
 @see
 @ref volReSample.cpp

 */

int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  typedef DGtal::functors::BasicDomainSubSampler< HyperRectDomain<SpaceND<3, int> >,  
                                                  DGtal::int32_t, double >   ReSampler; 
  typedef DGtal::ConstImageAdapter<Image3D, Image3D::Domain, ReSampler,
				   Image3D::Value,  DGtal::functors::Identity >  SamplerImageAdapter;


  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.vol"};
  std::vector<double> aGridSizeReSample;
  
  app.description("Re sample a 3D volumetric image (.vol, .longvol, .pgm3d)  with a given grid size. \n Example:\n to re sample an image with scale x,y,z  = 0.98, 0.98, 5.0,  you can do:\n volResSample -i image3d.vol -g 1 1 2  -o imageReSampled.vol \n " 
      "Note that if DGtal is compiled with the  option DGTAL_WITH_ITK set to ON, you can export the image in format ITK format and integrating image spacing.");
  app.add_option("-i,--input,1", inputFileName, "input volumetric file (.vol, .longvol, .pgm3d)." )
  ->required()
  ->check(CLI::ExistingFile);
  app.add_option("-o,--output,2", outputFileName, "the new volumetric file (.vol, .longvol, .pgm3d)." );
  app.add_option("-g,--gridSize", aGridSizeReSample, "size_x size_y size_z : the grid size of the re sampling ")
   ->expected(3);
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  trace.info()<< "Importing volume file :  " << inputFileName<< " ... " ;
  Image3D input3dImage = GenericReader<Image3D>::import(inputFileName);
  trace.info()<< "[done]" << endl;
  
  PointVector<3,int> shiftVector3D(0 ,0, 0);      
  DGtal::functors::BasicDomainSubSampler< HyperRectDomain<SpaceND<3, int> >,  
                                          DGtal::int32_t, double > reSampler(input3dImage.domain(),
                                                                             aGridSizeReSample,    shiftVector3D);
  const functors::Identity aFunctor{};
  SamplerImageAdapter sampledImage ( input3dImage, reSampler.getSubSampledDomain(), reSampler, aFunctor );
#ifdef DGTAL_WITH_ITK
    const std::string ext = outputFileName.substr( outputFileName.find_last_of(".") + 1 );
    if (std::find(ITK_IO_IMAGE_EXT.begin(),
                  ITK_IO_IMAGE_EXT.end(), ext) != ITK_IO_IMAGE_EXT.end() )
    {
        ITKWriter<SamplerImageAdapter>::exportITK(outputFileName, sampledImage,
                                                  Z3i::RealPoint(aGridSizeReSample[0],
                                                                 aGridSizeReSample[1],
                                                                 aGridSizeReSample[2]));
    }
    else
    {
        GenericWriter<SamplerImageAdapter>::exportFile(outputFileName, sampledImage);
    }

#else
  GenericWriter<SamplerImageAdapter>::exportFile(outputFileName, sampledImage);
#endif
  
  return 0;
}
