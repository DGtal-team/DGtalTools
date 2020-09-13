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
 * @file volFillInterior.cpp
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2018/09/01
 *
 *
 * This file is part of the DGtal library.
 */


/**
 @page volFillInterior volFillInterior
 
 @brief  Fills the interior of a voxel shape. The process can be sketched as follows: First the volume is filled in a
 breath-first approach from the point (0,0,0) (supposed to be exterior) using 6-adjacency. Then the complement is returned.
 
 @b Usage:  volFillInterior --input \<volFileName\> -o \<volOutputFileName\> (vol, longvol, p3d format)
 
 
 
 @b Allowed @b options @b are :
 @code
 Positionals:
   1 TEXT:FILE REQUIRED                  Input vol file.

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         Input vol file.
   -o,--output TEXT=result.vol           Output filename.
 @endcode
 
 @b Example:
 
 @code
 $ volFlip -i ${DGtal}/examples/samples/lobster.vol -o filled.vol
 @endcode
 
 @see
 @ref volFillInterior
 
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
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam ( const std::string &param )
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

  app.description("Fill the interior of a voxel set by filling the exterior using the 6-adjacency.\nThe exterior is the set of voxels with value zero and the interior voxels have value 128\n Basic usage:\n\tvolFillInterior --input <volFileName> --o <volOutputFileName> ");
  
  app.add_option("-i,--input,1", inputFileName, "Input vol file." )
  ->required()
  ->check(CLI::ExistingFile);
  app.add_option("-o,--output",outputFileName, "Output filename.", true);
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  trace.beginBlock("Loading");
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;
  MyImageC  image = VolReader< MyImageC >::importVol ( inputFileName );
  trace.info() << image << std::endl;
  trace.endBlock();
  
  //Flag image
  ImageContainerBySTLVector<Z3i::Domain, bool>  imageFlag(image.domain());
  for(auto &p: imageFlag.domain())
  {
    if (image(p) != 0)
      imageFlag.setValue(p, true);
  }
  
  std::stack<Z3i::Point> pstack;
  pstack.push(*(image.domain().begin()));
  FATAL_ERROR_MSG(image(pstack.top())==0, "Starting point of the domain must be equal to zero.");
  
  //6-Pencil for the exterior propagation
  std::vector<Z3i::Point> pencil6= { {1,0,0}, {0,1,0}, {0,0,1},{-1,0,0}, {0,-1,0}, {0,0,-1}  };
  
  trace.beginBlock("Filling");
  while (!pstack.empty())
  {
    Z3i::Point p = pstack.top();
    pstack.pop();
    imageFlag.setValue(p, true);

    for(auto & delta: pencil6)
      if ((image.domain().isInside(p + delta)) &&
          (imageFlag( p + delta) == false))
        pstack.push( p + delta);
  }
  trace.endBlock();
  
  trace.beginBlock("Complement");
  for(auto &p :  image.domain())
    if ((image(p) == 0) && (!imageFlag(p)))
      image.setValue(p,128);
  trace.endBlock();
  
  trace.beginBlock("Saving");
  bool res =  VolWriter< MyImageC >::exportVol(outputFileName, image);
  trace.endBlock();
  
  if (res) return 0; else return 1;
}
