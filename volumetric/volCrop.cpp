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
 * @file volCrop.cpp
 * @ingroup volumetric/voltools
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2013/05/05
 *
 * 
 *
 * This file is part of the DGtalTools.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/VolReader.h>
#include <DGtal/io/writers/VolWriter.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ConstImageAdapter.h>

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;
using namespace Z3i;

/**
 @page volCrop volCrop
 
 @brief  Crops a 3D vol image from domain coordinates.

 @b Usage: 	 volCrop --input \<volFileName\> --o \<volOutputFileName\> (both files can be independently in vol, pgm3D, p3d format)


 @b Allowed @b options @b are : 
 @code

 Positionals:
   1 TEXT:FILE REQUIRED                  Input vol file.

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         Input vol file.
   -o,--output TEXT=output.vol           Output filename.
   --xMin INT=0                          x coordinate of lower point.
   --yMin INT=0                          y coordinate of lower point.
   --zMin INT=0                          z coordinate of lower point.
   --xMax INT REQUIRED                   x coordinate of upper point.
   --yMax INT REQUIRED                   y coordinate of upper point.
   --zMax INT REQUIRED                   z coordinate of upper point.
 @endcode

 @b Example: 

 @code
$ volCrop --xMin 50 --yMin 50 --zMin 10 --xMax 150 --yMax 150 --zMax 50 -i ${DGtal}/examples/samples/lobster.vol -o croppedLobster.vol 
$ 3dImageViewer -i croppedLobster.vol
 @endcode


 You should obtain such a visualization:
 @image html resVolCrop.png "Resulting visualization."
 
 @see
 @ref volCrop.cpp

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
  exit ( 1 );
}

int main(int argc, char**argv)
{
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"output.vol"};
  int xMin {0};
  int yMin {0};
  int zMin {0};
  int xMax;
  int yMax;
  int zMax;
  
  app.description("Crops a 3D vol image from domain coordinates.\n Basic usage: \n \t volCrop --input <volFileName> --o <volOutputFileName> (both files can be independently in vol, pgm3D, p3d format)\nExample:\n volCrop --xMin 50 --yMin 50 --zMin 10 --xMax 150 --yMax 150 --zMax 50 -i ${DGtal}/examples/samples/lobster.vol -o croppedLobster.vol \n");
  
  app.add_option("-i,--input,1", inputFileName, "Input vol file." )
  ->required()
  ->check(CLI::ExistingFile);
  app.add_option("--output,-o", outputFileName, "Output filename.", true);

  app.add_option("--xMin",xMin, "x coordinate of lower point.", true);
  app.add_option("--yMin",yMin, "y coordinate of lower point.", true);
  app.add_option("--zMin",zMin, "z coordinate of lower point.", true);

  app.add_option("--xMax",xMax, "x coordinate of upper point.")
   ->required();
  app.add_option("--yMax",yMax, "y coordinate of upper point.")
   ->required();
  app.add_option("--zMax",zMax, "z coordinate of upper point.")
   ->required();
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------
  
  Z3i::Point ptLow( xMin, yMin, zMin);
  Z3i::Point ptMax(xMax , yMax, zMax);
  
  trace.beginBlock("Loading file");
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;
  MyImageC  imageC = VolReader< MyImageC >::importVol ( inputFileName );
  functors::Identity df;  
  
  typedef ConstImageAdapter<MyImageC, Domain, functors::Identity, MyImageC::Value, functors::Identity > ConstImageAdapterForSubImage;
  Domain subDomain(ptLow, ptMax);
  ConstImageAdapterForSubImage subImage(imageC, subDomain, df, df);
  trace.endBlock();

  trace.beginBlock("Exporting...");
  bool res =  VolWriter< ConstImageAdapterForSubImage>::exportVol(outputFileName, subImage);
  trace.endBlock();
  if (res) return 0; else return 1;
}
