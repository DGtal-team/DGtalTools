
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
 * @file volTrValues.cpp
 * @ingroup volumetric/
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2014/08/07
 *
 * 
 *
 * This file is part of the DGtalTools.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ConstImageAdapter.h>

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;
using namespace Z3i;

 
/**
 @page volTrValues volTrValues
 
 @brief  Applies basic vol image transform from the input values to output values.

 @b Usage:  	 volTrValues --input <volFileName> --o <volOutputFileName> -s 1 99 -r 100 200  

=> all voxels of values 1 (resp. 99) will be 100 (resp. 200) in the resulting image.   


 @b Allowed @b options @b are : 
 @code

 Positionals:
   1 TEXT:FILE REQUIRED                  Input vol file.

 Options:

 Positionals:
   1 TEXT:FILE REQUIRED                  Input vol file.
   2 TEXT=result.vol                     Output filename.

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         Input vol file.
   -o,--output TEXT=result.vol           Output filename.
   -s,--inputVals UINT ... REQUIRED      specify the values which will be transformed with the output values (given with --outputVals).
   -r,--outputVals UINT ... REQUIRED     specify the values which will be transformed with the output values (given with --outputVals).
   
 @endcode

 @b Example: 

 This tool can be useful to apply simple intensity transforms. For
 instance if you want to transform all intensities starting from 0 to 50 into interval 200 250 you can do as follows:

 @code
 $ volTrValues -i $DGtal/examples/samples/lobster.vol -s {0..50} -r {200..250} -o lobsterTr.vol
 $ 3dImageViewer -i  lobsterTr.vol
 @endcode

 By using  @ref Doc3dImageViewer ou should obtain such a result:
 @image html resVolTrValues.png "Result visualization."
 
 @see
 @ref volTrValues.cpp

 */


int main(int argc, char**argv)
{

  // parse command line ----------------------------------------------
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.vol"};
  std::vector<unsigned int> inputVals;
  std::vector<unsigned int>  outputVals;

  app.description("Apply basic vol image transform from the input values to output values.\n Basic usage:\n \t volTrValues --input <volFileName> --o <volOutputFileName> -s 1 99 -r 100 200 \n\t => all voxel of values 1 (resp. 99) will be 100 (resp. 200) in the resulting image.");
  
  app.add_option("-i,--input,1", inputFileName, "Input vol file." )
  ->required()
  ->check(CLI::ExistingFile);
  
  app.add_option("--output,-o,2",outputFileName, "Output filename.", true);
  app.add_option("--inputVals,-s", inputVals, "specify the values which will be transformed with the output values (given with --outputVals).") ->required();
  app.add_option("--outputVals,-r", outputVals, "specify the values which will be transformed with the output values (given with --outputVals).") ->required();
  
 
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------


  if(inputVals.size()!=outputVals.size()){
    trace.error()<< "Transformation not possible the two sets of input/output values should have the same size." << std::endl;
    exit(1);
  }
  
  trace.beginBlock("Loading file");
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;
  MyImageC  image = GenericReader< MyImageC >::import( inputFileName );
  trace.endBlock();  
  unsigned int val;
  for(MyImageC::Domain::ConstIterator it = image.domain().begin(),
        itend = image.domain().end(); it != itend; ++it)
    {
      val = image(*it);
      for(unsigned int i = 0; i< inputVals.size(); i++){
        if(inputVals.at(i)==val){
          image.setValue( *it , outputVals.at(i));
        }
      }
    } 

  trace.beginBlock("Exporting...");
  bool res =  GenericWriter<MyImageC>::exportFile(outputFileName, image);
  trace.endBlock();

  if (res) return 0; else return 1;
}
