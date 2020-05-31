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
 * @file dicom2vol.cpp
 * @ingroup conerters
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2013/10/30
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
#include "DGtal/io/readers/DicomReader.h"


#include "CLI11.hpp"


using namespace std;
using namespace DGtal;


/**
 @page dicom2vol dicom2vol
 @brief Converts dicom file into a volumetric file (.vol, .longvol .pgm3d).

@b Usage: dicom2vol [input] [output]

@b Allowed @b options @b are:

@code
Positionals:
  1 TEXT:FILE REQUIRED                  dicom image  (.dcm).

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         dicom image  (.dcm).
  -o,--output TEXT                      volumetric file (.vol, .longvol .pgm3d, .raw)
  --dicomMin INT                        set minimum density threshold on Hounsfield scale
  --dicomMax INT                        set maximum density threshold on Hounsfield scale
@endcode

@b Example:
@code  
$  dicom2vol -i ${DGtal}/tests/samples/dicomSample/1629.dcm --dicomMin 0 --dicomMax 300 -o sample.vol
@endcode

@see dicom2vol.cpp

*/



int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  
   // parse command line using CLI ----------------------------------------------
   CLI::App app;
   std::string inputFileName;
   std::string outputFileName {"result.raw"};
   DGtal::int64_t dicomMin {-1000};
   DGtal::int64_t dicomMax {3000};

   app.description("Convert dicom file into a volumetric file (.vol, .longvol .pgm3d).\n Example:\n dicom2vol -i ${DGtal}/tests/samples/dicomSample/1629.dcm --dicomMin 0 --dicomMax 300 -o sample.vol.");
   app.add_option("-i,--input,1", inputFileName, "dicom image  (.dcm)." )
    ->required()
    ->check(CLI::ExistingFile);
   app.add_option("-o,--output,2", outputFileName, "volumetric file (.vol, .longvol .pgm3d, .raw)", true);
   app.add_option("--dicomMin",dicomMin,"set minimum density threshold on Hounsfield scale" );
   app.add_option("--dicomMax",dicomMax,"set maximum density threshold on Hounsfield scale" );

   app.get_formatter()->column_width(40);
   CLI11_PARSE(app, argc, argv);
   // END parse command line using CLI ----------------------------------------------
   
   typedef DGtal::functors::Rescaling<int ,unsigned char > RescalFCT;
  
   trace.info() << "Reading input dicom file " << inputFileName ; 
  Image3D inputImage = DicomReader< Image3D,  RescalFCT  >::importDicom(inputFileName, 
									RescalFCT(dicomMin,dicomMax, 0, 255) );
  trace.info() << " [done] " << std::endl ; 
  trace.info() << " converting into vol file... " ; 
  inputImage >> outputFileName; 
  trace.info() << " [done] " << std::endl ;   

  return EXIT_SUCCESS;
}




