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
 * @file freeman2sdp.cpp
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


//contour
#include "DGtal/geometry/curves/FreemanChain.h"


//boost
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

//STL
#include <vector>
#include <string>

using namespace DGtal;

/**
 @page freeman2sdp freeman2sdp
 @brief Transform freeman chain into a Sequence of Discrete Points. Result is given to std output.

 @b Usage: freeman2sdp [input] > output.sdp

@b Allowed @b options @b are:

@code
  -h [ --help ]         display this message
  -i [ --input ] arg    Input freeman chain file name.
  --info                adds some info as comments at the beginning of the 
                        file.
  -o [ --oneLine ]       output the digital contour in one line like: X0 Y0 X1 
                        Y1 ... XN YN
@endcode

@b Example:
@code
  freeman2sdp -i ${DGtal}/tests/samples/contourS.fc > contourS.sdp
@endcode
You will obtain such result:
@verbatim
$ more contourS.sdp  
# grid curve 1/1 closed
13 60
14 60
14 59
14 58
15 58
15 57
16 57
16 56
17 56
17 55
17 54
18 54
...
@endverbatim
@see freeman2sdp.cpp

*/



///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "Input freeman chain file name.")
    ("info", "adds some info as comments at the beginning of the file.") 
    ("oneLine,o", " output the digital contour in one line like: X0 Y0 X1 Y1 ... XN YN");

  
  
  typedef FreemanChain<Z2i::Integer> FreemanChain;
  
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }

  po::notify(vm);    
  if(!parseOK||vm.count("help")||argc<=1 || (!(vm.count("input")) ) )
    {
      trace.info()<< "Transform freeman chain into a Sequence of Discrete Points. Result is given to std output. " <<std::endl << "Basic usage: "<<std::endl
		  << "\t freeman2sdp [input] > out.sdp  "<<std::endl
		  << general_opt << "\n";
      trace.info() << "Example:\n"
                   << "freeman2sdp -i ${DGtal}/tests/samples/contourS.fc > contourS.sdp \n";

      return 0;
    }
  
  bool oneline = vm.count("oneLine");
  if( vm.count("input") ){
    std::string fileName = vm["input"].as<std::string>();
    std::vector< FreemanChain > vectFcs =  PointListReader< Z2i::Point >:: getFreemanChainsFromFile<Z2i::Integer> (fileName); 
    
    for(unsigned int i=0; i< vectFcs.size(); i++){
      bool isClosed = vectFcs.at(i).isClosed(); 
      std::cout << "# grid curve " << i+1 << "/" << vectFcs.size()
                << ( (isClosed)?" closed":" open" ) << std::endl;
      if ( vm.count("info"))
        std::cout << "# SDP contour" << i+1<< "/" << vectFcs.size() << " "
                  << "# size=" << vectFcs.at(i).size() << std::endl;
      std::vector<Z2i::Point> vectPts; 
      FreemanChain::getContourPoints( vectFcs.at(i), vectPts ); 
      for(unsigned int k=0; k < vectPts.size(); k++){
        std::cout << vectPts.at(k)[0] << " "<< vectPts.at(k)[1] ;
	if(!oneline){
          std::cout << std::endl;
	}else
          {
            std::cout << " "; 
          }
          
      }
      std::cout << std::endl;
    }
    
  }
  

}

