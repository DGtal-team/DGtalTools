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
 * @file vol2obj.cpp
 * @ingroup Converters
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr)
 *
 * @date 2013/10/13
 *
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/readers/PointListReader.h"

#include "DGtal/images/ImageSelector.h"

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;
using namespace Z3i;

/**
 @page vol2obj vol2obj
 @brief Converts any volumetric file to an OBJ one. Each grid point with value
between
 [@a thresholdMin,@a thresholdMax] is exported as a unit cube.


@b Usage: ./converters/vol2obj [OPTIONS] 1 [2]

@b Allowed @b options @b are:

@code

Positionals:
  1 TEXT:FILE REQUIRED                  vol file (.vol, .longvol .p3d, .pgm3d or .sdp and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
  2 TEXT                                output file (.obj or .off).

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         vol file (.vol, .longvol .p3d, .pgm3d or .sdp and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
  -o,--output TEXT                      output file (.obj or .off).
  -m,--thresholdMin INT=128             threshold min (excluded) to define binary shape.
  -M,--thresholdMax INT=255             threshold max (included) to define binary shape.
  --rescaleInputMin INT=0               min value used to rescale the input intensity (to avoid basic cast into 8  bits image).
  --rescaleInputMax INT=255             max value used to rescale the input intensity (to avoid basic cast into 8  bits image).

@endcode

@see
vol2obj.cpp

*/
///////////////////////////////////////////////////////////////////////////////


int main( int argc, char** argv )
{

  // parse command line using CLI ----------------------------------------------
   CLI::App app;
   std::string inputFileName;
   std::string outputFileName {"result.obj"};
   int thresholdMin {128};
   int thresholdMax {255};  
   DGtal::int64_t rescaleInputMin {0};
   DGtal::int64_t rescaleInputMax {255};

   
  app.description(" Converts any volumetric file (or .sdp file) to an OBJ one. Each grid point with value between [thresholdMin, thresholdMax] is exported as a unit cube.");
  app.add_option("-i,--input,1", inputFileName, "vol file (.vol, .longvol .p3d, .pgm3d or .sdp and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("--output,-o,2",outputFileName ,"output file (.obj or .off).");
  app.add_option("--thresholdMin,-m", thresholdMin, "threshold min (excluded) to define binary shape.", true);
  app.add_option("--thresholdMax,-M", thresholdMax, "threshold max (included) to define binary shape.", true);
  app.add_option("--rescaleInputMin", rescaleInputMin, "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);
  app.add_option("--rescaleInputMax", rescaleInputMax, "max value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);


  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------


  Board3D<> board;

  typedef ImageSelector<Domain, unsigned char>::Type Image;
  string extension = inputFileName.substr(inputFileName.find_last_of(".") + 1);

  if(extension!="sdp")
    {
      typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
      Image image =  GenericReader< Image >::importWithValueFunctor( inputFileName,RescalFCT(rescaleInputMin,
                                                                                             rescaleInputMax,
                                                                                             0, 255) );
      trace.info() << "Image loaded: "<<image<< std::endl;
      Domain domain = image.domain();
      for(Domain::ConstIterator it = domain.begin(), itend=domain.end(); it!=itend; ++it){
        unsigned char  val= image( (*it) );
        if(val<=thresholdMax && val >=thresholdMin){
          board << *it;
        }
      }
    }
  else
    if(extension=="sdp")
      {
        vector<Z3i::Point> vectVoxels = PointListReader<Z3i::Point>::getPointsFromFile(inputFileName);
        for(unsigned int i=0;i< vectVoxels.size(); i++){
          board << vectVoxels.at(i);
        }
      }

  board.saveOBJ(outputFileName);
  return EXIT_SUCCESS;
}
