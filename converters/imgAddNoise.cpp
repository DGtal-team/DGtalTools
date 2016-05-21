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
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <vector>
#include <string>
#include <climits>

 
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;


/**
 @page imgAddNoise imgAddNoise
 @brief  Adds noise (Kanungo's) to a binary 2D object.

@b Usage: imgAddNoise [input] [output]

@b Allowed @b options @b are:

@code
  -h [ --help ]             display this message
  -i [ --input ] arg        input image file name (any 2D image format accepted
                            by DGtal::GenericReader)
  -o [ --output ] arg       output image file name (any 2D image format 
                            accepted by DGtal::GenericWriter)
  -n [ --noise ] arg (=0.5) Kanungo noise level in ]0,1[ (default 0.5)
                                      resulting volumetric file.
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
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "input image file name (any 2D image format accepted by DGtal::GenericReader)")
    ("output,o", po::value<std::string>(), "output image file name (any 2D image format accepted by DGtal::GenericWriter)")
    ("noise,n", po::value<double>()->default_value(0.5), "Kanungo noise level in ]0,1[ (default 0.5)")  ;
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
    parseOK=false;
  }
  po::notify(vm);    
  if(vm.count("help")||argc<=1|| !parseOK)
    {
      trace.info()<< "Add Kanungo noise to a binary object with 0 values "
                  << "as background points and values >0 for the foreground ones."
                  <<std::endl << "Basic usage: "<<std::endl
                  << "\t imgAddNoi0se [options] --input <imageName> --output <outputImage>"
                  << "-noise 0.3" <<std::endl
                  << general_opt << "\n"
                  << "Example: \n"
                  << "imgAddNoise -i ${DGtal}/examples/samples/klokan.pgm -o noise.pgm "
                  << std::endl;

      return 0;
    }
  
  //Parameters
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  const std::string input = vm["input"].as<std::string>();
  if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
  const std::string  output = vm["output"].as<std::string>();
  const double noise = vm["noise"].as<double>();
  
  typedef functors::IntervalForegroundPredicate<MyImage> Binarizer;
  MyImage image = GenericReader<MyImage>::import( input );
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

  result >> output;
  
  return 0;
}



