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

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;
using namespace Z3i;

/**
 @page vol2obj vol2obj
 @brief Converts any volumetric file to an OBJ one. Each grid point with value between
 [@a thresholdMin,@a thresholdMax] is exported as a unit cube.


@b Usage: vol2obj [input-file]

@b Allowed @b options @b are:

@code
   -h [ --help ]                    display this message
   -i [ --input ] arg               vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm
                                    (with 3 dims)) file or sdp (sequence of
                                    discrete points)
   -o [ --output ] arg              Output OBJ filename
   -m [ --thresholdMin ] arg (=0)   threshold min to define binary shape
   -M [ --thresholdMax ] arg (=255) threshold max to define binary shape
@endcode

@see
vol2obj.cpp

*/
///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
    ("output,o", po::value<std::string>(), "Output OBJ filename" )
    ("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min to define binary shape" )
    ("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max to define binary shape" );

  bool parseOK=true;
  po::variables_map vm;
  try
    {
      po::store(po::parse_command_line(argc, argv, general_opt), vm);
    } catch(const std::exception& ex)
    {
      parseOK=false;
      trace.info()<< "Error checking program options: "<< ex.what()<< endl;
    }
  po::notify(vm);
  if( !parseOK || vm.count("help")||argc<=1)
    {
      std::cout << "Usage: " << argv[0] << " [input-file]\n"
                << "Convert a  volume file into OBJ format\n"
                << general_opt << "\n";
      return 0;
    }

  if(! vm.count("input"))
    {
      trace.error() << " The input filename was defined" << endl;
      return 0;
    }
  if(! vm.count("output"))
    {
      trace.error() << " The output filename was defined" << endl;
      return 0;
    }

  string inputFilename = vm["input"].as<std::string>();
  string outputFilename = vm["output"].as<std::string>();
  int thresholdMin = vm["thresholdMin"].as<int>();
  int thresholdMax = vm["thresholdMax"].as<int>();

  Board3D<> board;

  typedef ImageSelector<Domain, unsigned char>::Type Image;
  string extension = inputFilename.substr(inputFilename.find_last_of(".") + 1);
  if(extension!="vol" && extension != "p3d" && extension != "pgm3D" &&
     extension != "pgm3d" && extension != "sdp" && extension != "pgm")
    {
      trace.info() << "File extension not recognized: "<< extension << std::endl;
      return 0;
    }

  if(extension=="vol" || extension=="pgm3d" || extension=="pgm3D")
    {
      Image image = GenericReader<Image>::import (inputFilename );
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
        vector<Z3i::Point> vectVoxels = PointListReader<Z3i::Point>::getPointsFromFile(inputFilename);
        for(unsigned int i=0;i< vectVoxels.size(); i++){
          board << vectVoxels.at(i);
        }
      }


  board.saveOBJ(outputFilename);
  return 0;
}
