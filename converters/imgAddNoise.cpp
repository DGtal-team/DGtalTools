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
 * @file imgAddNoise
 * @ingroup converters
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr)
 *
 * @date 2015/03/24
 *
 *
 * This file is part of the DGtal library.
 */


#include <DGtal/base/Common.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageSelector.h"
#include <DGtal/geometry/volumes/KanungoNoise.h>
#include <DGtal/images/IntervalForegroundPredicate.h>

#include "CLI11.hpp"

#include <vector>
#include <string>
#include <climits>



using namespace DGtal;


/**
 @page imgAddNoise imgAddNoise
 @brief  Adds noise (Kanungo's) to a binary 2D object.

@b Usage: imgAddNoise [input] [output]

@b Allowed @b options @b are:

@code

Positionals:
  1 TEXT:FILE REQUIRED                  input image file name (any 2D image format accepted by DGtal::GenericReader).
  2 TEXT=result.png                     output image file name (any 2D image format accepted by DGtal::GenericWriter)

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         input image file name (any 2D image format accepted by DGtal::GenericReader).
  -o,--output TEXT=result.png           output image file name (any 2D image format accepted by DGtal::GenericWriter)
  -n,--noise FLOAT=0.5                  Kanungo noise level in ]0,1[ (default 0.5)
@endcode

@b Example:
@code
  $ imgAddNoise -i ${DGtal}/examples/samples/klokan.pgm -o noise.pgm 

@endcode
You will obtain such image:
@image html  resImgAddNoise.png "Resulting image."
@see imgAddNoise.cpp

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

typedef ImageSelector < Z2i::Domain, unsigned char>::Type MyImage;


int main( int argc, char** argv )
{
// parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.png"};
  double noise {0.5};

  app.description("Add Kanungo noise to a binary object with 0 values as background points and values >0 for the foreground ones.\n Example: \n imgAddNoise -i ${DGtal}/examples/samples/klokan.pgm -o noise.pgm ");
  app.add_option("-i,--input,1", inputFileName, "input image file name (any 2D image format accepted by DGtal::GenericReader)." )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("-o,--output,2", outputFileName, "output image file name (any 2D image format accepted by DGtal::GenericWriter)", true);
  app.add_option("-n,--noise", noise, "Kanungo noise level in ]0,1[ (default 0.5)", true);

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

    
  typedef functors::IntervalForegroundPredicate<MyImage> Binarizer;
  MyImage image = GenericReader<MyImage>::import( inputFileName );
  trace.info() <<"Input image: "<< image<<std::endl;
  Binarizer predicate(image, 0,255);
  
  KanungoNoise<Binarizer, Z2i::Domain> kanungo(predicate, image.domain(), noise);
  
  MyImage result(image.domain());
  for(Z2i::Domain::ConstIterator it = image.domain().begin(), itend = image.domain().end(); it!= itend; ++it)
  {
    if (kanungo(*it))
      result.setValue(*it, 255);
    else
      result  .setValue(*it, 0);
  }

  result >> outputFileName;
  
  return EXIT_SUCCESS;
}



