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
 * @file vol2slice.cpp
 * @ingroup Converters
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
#include <boost/format.hpp>

#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////

/**
 @page vol2slice vol2slice
 @brief  Convert a volumetric file (.vol, .longvol, .pgm3d) into a set of 2D slice  images.

@b Usage: vol2slice [input] [output]

@b Allowed @b options @b are:

@code
 Typical use: to extract all slices defined in Y plane (y=cst): 
 vol2slice -i image3d.vol -s 1  -o slice.pgm 

Usage: ./converters/vol2slice [OPTIONS] 1 [2]

Positionals:
  1 TEXT:FILE REQUIRED                  vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha                                         or mhd formats, the input values are linearly scaled between 0 and 255.
  2 TEXT=result.pgm                     base_name.extension:  extracted 2D slice volumetric files (will result n files base_name_xxx.extension)

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha                                         or mhd formats, the input values are linearly scaled between 0 and 255.
  -o,--output TEXT=result.pgm           base_name.extension:  extracted 2D slice volumetric files (will result n files base_name_xxx.extension)
  -f,--setFirstSlice INT:NUMBER=0       Set the first slice index to be extracted.
  -l,--setLastSlice INT:NUMBER          Set the last slice index to be extracted (by default set to maximal value according to the given volume).
  -s,--sliceOrientation UINT:{0,1,2}=2  specify the slice orientation for which the slice are defined (by default =2 (Z direction))
  --rescaleInputMin INT=0               min value used to rescale the input intensity (to avoid basic cast into 8  bits image).
  --rescaleInputMax INT=255             max value used to rescale the input intensity (to avoid basic cast into 8  bits image).



@endcode

@b Example:
@code 
   # Export Z slice images (-s 2): 
   $ vol2slice -i ${DGtal}/examples/samples/lobster.vol -o slice.pgm  -f 10 -l 15 -s 2
@endcode

You should obtain such a visualization:
@image html resVol2slice.png "resulting visualisation."

@see
@ref vol2slice.cpp

*/


int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char > Image2D;
  typedef DGtal::ConstImageAdapter<Image3D, Image2D::Domain, DGtal::functors::Projector< DGtal::Z3i::Space>,
				   Image3D::Value,  DGtal::functors::Identity >  SliceImageAdapter;

  CLI::App app;


  std::string inputFileName;
  std::string outputFileName = "result.pgm";
  DGtal::int64_t rescaleInputMin {0};
  DGtal::int64_t rescaleInputMax {255};
  int userStartSlice {0};
  int userEndSlice;
  unsigned int sliceOrientation {2};
  

  app.description("Convert a volumetric file (.vol, .longvol, .pgm3d) into a set of 2D slice  images. \n Typical use: to extract all slices defined in Y plane (y=cst): \n vol2slice -i image3d.vol -s 1  -o slice.pgm \n");

  app.add_option("-i,--input,1", inputFileName, "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.")
    ->required()
    ->check(CLI::ExistingFile);
  
  app.add_option("--output,-o,2",outputFileName ,"base_name.extension:  extracted 2D slice volumetric files (will result n files base_name_xxx.extension)", true);
  app.add_option("--setFirstSlice,-f", userStartSlice, "Set the first slice index to be extracted.", true)
    -> check(CLI::Number);
  app.add_option("--setLastSlice,-l", userEndSlice, "Set the last slice index to be extracted (by default set to maximal value according to the given volume).")
    -> check(CLI::Number);
  app.add_option("--sliceOrientation,-s", sliceOrientation, "specify the slice orientation for which the slice are defined (by default =2 (Z direction))", true)
    -> check(CLI::IsMember({0, 1, 2}));
  app.add_option("--rescaleInputMin", rescaleInputMin, "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);
  app.add_option("--rescaleInputMax", rescaleInputMax, "max value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);
  app.get_formatter()->column_width(40);

  CLI11_PARSE(app, argc, argv);
  
  std::string outputExt = outputFileName.substr(outputFileName.find_last_of(".")+1);
  std::string outputBasename = outputFileName.substr(0, outputFileName.find_last_of("."));
     

  trace.info()<< "Importing volume file base name:  " << outputBasename << " extension: " << outputExt << " ..." ;
  typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
  Image3D input3dImage =  GenericReader< Image3D >::importWithValueFunctor( inputFileName,RescalFCT(rescaleInputMin,
                                                                                                    rescaleInputMax,
                                                                                                    0, 255) );
  
  trace.info()<< "[done]" << endl;
  
  unsigned int startSlice=0;
  unsigned int endSlice=input3dImage.domain().upperBound()[sliceOrientation];

  if(userStartSlice !=0){
    startSlice = userStartSlice;
  }
  if(userEndSlice != 0){
    endSlice = userEndSlice;
  }

  //Processing each slice
#pragma omp parallel for schedule(dynamic)
   for( unsigned int i=startSlice; i <= endSlice; i++){
    trace.info() << "Exporting slice image "<< i ;
    DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(sliceOrientation);
    DGtal::Z2i::Domain domain2D(invFunctor(input3dImage.domain().lowerBound()),
  				invFunctor(input3dImage.domain().upperBound()));
    DGtal::functors::Projector<DGtal::Z3i::Space> aSliceFunctor(i); aSliceFunctor.initAddOneDim(sliceOrientation);
    const DGtal::functors::Identity identityFunctor{};
    SliceImageAdapter sliceImage( input3dImage, domain2D, aSliceFunctor, identityFunctor );
    stringstream outName; outName << outputBasename << "_" <<  boost::format("%|05|")% i <<"."<< outputExt ;
    trace.info() << ": "<< outName.str() ;
    GenericWriter<SliceImageAdapter>::exportFile(outName.str(), sliceImage);
    trace.info() << " [done]"<< endl;
  }


  return 0;
}
