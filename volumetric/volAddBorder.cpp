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
 * @file voAddBorder.cpp
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/05/01
 *
 * This file is part of the DGtal library.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/VolReader.h>
#include <DGtal/io/writers/VolWriter.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;
using namespace Z3i;

/**
 @page volAddBorder volAddBorder
 
 @brief Adds a border of one voxel with value 0 around a vol file.

 @b Usage: 	volAddBorder --input \<volFileName\> --o \<volOutputFileName\> 


 @b Allowed @b options @b are : 
 @code

 Positionals:
   1 TEXT:FILE REQUIRED                  Input vol file.

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         Input vol file.
   --inside                              Sets zero value to domain boundary voxels without changing the domain extent.
   -o,--output TEXT                      Output filename.

 @endcode

 @b Example: 

 @code
 $  volAddBorder -i $DGtal/examples/samples/Al.100.vol -o Al.100border.vol
 @endcode
 
 You should obtain an resulting image with domain:
-1, -1, -1 x 100, 100, 100 instead  0, 0, 0 x 99, 99, 99

 @see
 @ref volAddBorder.cpp

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
  std::string outputFileName {"result.vol"};
  bool addInside {false};

  app.description("Add a border of one voxel with value 0 around a vol file.\n \tvolAddBorder --input <volFileName> --o <volOutputFileName>\n Example: \n \t volAddBorder -i $DGtal/examples/samples/Al.100.vol -o Al.100border.vol");
  
  app.add_option("-i,--input,1", inputFileName, "Input vol file." )
  ->required()
  ->check(CLI::ExistingFile);
  app.add_flag("--inside", addInside , "Sets zero value to domain boundary voxels without changing the domain extent.");
  app.add_option("--output,-o,2", outputFileName, "Output filename.");
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;
  
  MyImageC  imageC = VolReader< MyImageC >::importVol ( inputFileName );
  Z3i::Domain rDom  ( imageC.domain().lowerBound() -  Vector().diagonal(addInside ? 0: 1),
                      imageC.domain().upperBound() + Vector().diagonal(addInside ? 0: 1));
  MyImageC  outputImage(rDom);
  Z3i::Domain iDom  ( imageC.domain().lowerBound() +  Vector().diagonal( 1),
                      imageC.domain().upperBound() - Vector().diagonal( 1));

  //Fast Copy
  for(MyImageC::Domain::ConstIterator it = imageC.domain().begin(),
          itend = imageC.domain().end(); it != itend; ++it){
     if(!addInside){
          outputImage.setValue( *it , imageC(*it));
     }
     else
     {
       if (!iDom.isInside(*it)){
         imageC.setValue( *it , 0); 
       }
     }      
     
  }  
  bool res = true;
  if (!addInside) {
    res=  VolWriter< MyImageC>::exportVol(outputFileName, outputImage);
  }
  else{
    res=  VolWriter< MyImageC>::exportVol(outputFileName, imageC);
  }
  if (res) return 0; else return 1;
}
