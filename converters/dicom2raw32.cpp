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
 * @file dicom2raw32.cpp
 * @ingroup Converters
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2015/5/29
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
#include "DGtal/io/writers/RawWriter.h"
#include "DGtal/io/readers/DicomReader.h"

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;

/**
   @ingroup convertertools
   @page dicom2raw32 dicom2raw32
   @section dicom2raw32_sec dicom2raw32
   @brief Convert dicom file into a volumetric raw file (.raw).

   @b Usage: dicom2raw32 [input] [output]

   @b Allowed @b options @b are:

@code
 -h [ --help ]           display this message
  -i [ --input ] arg      dicom image  (.dcm) 
  -o [ --output ] arg     volumetric file (.vol, .longvol .pgm3d) 
  --dicomMin arg (=-1000) set minimum density threshold on Hounsfield scale
  --dicomMax arg (=3000)  set maximum density threshold on Hounsfield scale
@endcode

@b Example:
@code  $dicom2raw32 -i ${DGtal}/tests/samples/dicomSample/1629.dcm --dicomMin -500 --dicomMax -100 -o sample.vol
@endcode

@see dicom2raw32.cpp

*/

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  
  typedef ImageContainerBySTLVector < Z3i::Domain, DGtal::uint32_t > Image3D;

  // parse command line using CLI ----------------------------------------------
   CLI::App app;
   std::string inputFileName;
   std::string outputFileName {"result.raw"};
   DGtal::int64_t dicomMin {-1000};
   DGtal::int64_t dicomMax {3000};

   app.description("Convert dicom file into a volumetric raw file (.raw). Example:\n dicom2raw32 -i ${DGtal}/tests/samples/dicomSample/1629.dcm --dicomMin -500 --dicomMax -100 -o sample.raw \n ");
   app.add_option("-i,--input,1", inputFileName, "dicom image  (.dcm)." )
    ->required()
    ->check(CLI::ExistingFile);
   app.add_option("-o,--output", outputFileName, "volumetric raw file (.raw)");
   app.add_option("--dicomMin",dicomMin,"set minimum density threshold on Hounsfield scale" );
   app.add_option("--dicomMax",dicomMax,"set maximum density threshold on Hounsfield scale" );


  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------
  
  
  typedef DGtal::functors::Rescaling<int ,DGtal::uint32_t > RescalFCT;

  trace.info() << "Reading input dicom file " << inputFileName ; 
  Image3D inputImage = DicomReader< Image3D,  RescalFCT  >::importDicom(inputFileName, 
									RescalFCT(dicomMin,dicomMax, 0, 
                                                                                  std::numeric_limits<DGtal::uint32_t>::max()) );
  trace.info() << " [done] " << std::endl ; 
  trace.info() << " converting into longvol file... " ; 
  DGtal::RawWriter<Image3D>::exportRaw32(outputFileName, inputImage) ; 
  trace.info() << " [done] " << std::endl ;   
  
  return EXIT_SUCCESS;
  
}




