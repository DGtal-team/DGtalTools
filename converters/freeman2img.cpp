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
 * @file freeman2img.cpp
 * @ingroup Tools
 * @author Bertrand Kerautret (\c kerautre@loria.fr) and Jacques-Olivier Lachaud 
 * LORIA (CNRS, UMR 7503), University of Nancy, France 
 * (backport from ImaGene)
 * @date 2012/19/05
 *
 * convert freeman chain to a Sequence of Discrete Points.  
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

//image
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/writers/GenericWriter.h"

//contour
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/SurfelSetPredicate.h"

//STL
#include <vector>
#include <string>

#include "CLI11.hpp"


using namespace DGtal;


/**
 @page freeman2img freeman2img
 @brief Transforms one or several freeman chains into a pgm file by filling their interior areas.

 The transformation can fill shapes with hole by using the freemanchain orientation. The interior is considered on the left according to a freeman chain move, i.e. a clockwise oriented contour represents a hole in the shape.

@b Usage: freeman2img [input] [output]

@b Allowed @b options @b are:

@code

Positionals:
  1 TEXT:FILE REQUIRED                  Input freeman chain file name.

Options:

Positionals:
  1 TEXT:FILE REQUIRED                  Input freeman chain file name.

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         Input freeman chain file name.
  -b,--border UINT                      add a border in the resulting image (used only in the automatic mode i.e when --space is not used.
  -o,--output TEXT=result.pgm           the output fileName
  -s,--space INT x 4                    Define the space from its bounding box (lower and upper coordinates) else the space is automatically defined from the freemanchain bounding boxes.
@endcode

@b Example:
@code
  $freeman2img -i ${DGtal}/tests/samples/contourS.fc -o sample.pgm
@endcode
You will obtain such image:
@image html  resFreeman2img.png "Resulting image"


@b Example @b with @b several @b contours:

The file located in $DGtal/examples/samples/contourS2.fc contains different contours with some ones corresponds to hole.  We can apply the same conversion as the previous example:

@code
$ freeman2img -i  $DGtal/examples/samples/contourS2.fc  -o sample2.pgm
@endcode

You will obtain such image:
@image html  resFreeman2img2.png "Resulting image"



@see @ref img2freeman
@ref freeman2img.cpp

*/


int main( int argc, char** argv )
{

  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName="result.pgm";
  unsigned int border {0};
  std::vector<int> space;

  app.description("Transform one or several freeman chains into an grayscale image file by filling its interior areas.\nThe transformation can fill shapes with hole by using the freemanchain orientation.  The interior is considered on the left according to a freeman chain move, i.e. a clockwise oriented contour represents a hole in the shape. Basic example:\n \t freeman2img  -i  inputChain.fc -o contourDisplay.pgm -b 5");    
  app.add_option("-i,--input,1", inputFileName, "Input freeman chain file name." )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("-b,--border",border, "add a border in the resulting image (used only in the automatic mode i.e when --space is not used.");
  app.add_option("-o,--output", outputFileName, "the output fileName", true);
  app.add_option("-s,--space", space, "Define the space from its bounding box (lower and upper coordinates) else the space is automatically defined from the freemanchain bounding boxes." )
    ->expected(4);

  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

    
  typedef KhalimskySpaceND<2, int>::Cell Cell;
  typedef KhalimskySpaceND<2, int>::SCell SCell;
  typedef FreemanChain<Z2i::Integer> FreemanChain;
  typedef DGtal::KhalimskySpaceND< 2, int > KSpace;
  typedef DGtal::ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image2D ; 
      
    
  std::vector< FreemanChain > vectFcs =  PointListReader< Z2i::Point >::getFreemanChainsFromFile<Z2i::Integer> (inputFileName);    
  int minx=std::numeric_limits<int>::max();
  int miny=std::numeric_limits<int>::max();
  int maxx=std::numeric_limits<int>::min();
  int maxy=std::numeric_limits<int>::min();

  
  if(space.size()!=4){
      for(std::vector< FreemanChain >::const_iterator it = vectFcs.begin(); it!= vectFcs.end(); it++){
        FreemanChain fc = *it;
        int t_minx=std::numeric_limits<int>::max();
        int t_miny=std::numeric_limits<int>::max();
        int t_maxx=std::numeric_limits<int>::min();
        int t_maxy=std::numeric_limits<int>::min();
        fc.computeBoundingBox(t_minx, t_miny, t_maxx, t_maxy);
        minx = t_minx > minx? minx: t_minx;
        miny = t_miny > miny? miny: t_miny;
        maxx = t_maxx < maxx? maxx: t_maxx;
        maxy = t_maxy < maxy? maxy: t_maxy;
      }
      minx-=border; miny-=border; maxx+=border;   maxy+=border;
    }else{
      minx = space[0];
      miny = space[1];
      maxx = space[2];
      maxy = space[3];      
    }
    KSpace aKSpace;
    aKSpace.init(Z2i::Point(minx, miny), Z2i::Point(maxx, maxy), true);
    std::set<SCell> boundarySCell;
    std::set<Cell> interiorCell;
    for(std::vector< FreemanChain >::const_iterator it = vectFcs.begin(); it!= vectFcs.end(); it++){
      FreemanChain fc = *it;
      FreemanChain::getInterPixelLinels(aKSpace, fc, boundarySCell, true);
    }
    
    Image2D imageResult (Z2i::Domain(Z2i::Point(minx, miny), Z2i::Point(maxx, maxy))); 
    Surfaces<KSpace>::uFillInterior(aKSpace, functors::SurfelSetPredicate<std::set<SCell>,SCell>(boundarySCell), imageResult, 255, false, false );  
    imageResult >> outputFileName;
  

}

