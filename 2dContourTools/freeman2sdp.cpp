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
 * @author Bertrand Kerautret (\c kerautre@loria.fr)
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
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
#include "DGtal/geometry/curves/representation/FreemanChain.h"


//boost
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

//STL
#include <vector>
#include <string>

using namespace DGtal;




///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("info,i", "adds some info as comments at the beginning of the file.") 
    ("oneLine,o", " output the digital contour in one line like: X0 Y0 X1 Y1 ... XN YN") 
    ("FreemanChain,f", po::value<std::string>(), "FreemanChain file name");
  
  
  typedef FreemanChain<Z2i::Integer> FreemanChain;
  
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }

  po::notify(vm);    
  if(!parseOK||vm.count("help")||argc<=1 || (!(vm.count("FreemanChain")) ) )
    {
      trace.info()<< "Transform freeman chain into a Sequence of Discrete Points. Result is given to std output " <<std::endl << "Basic usage: "<<std::endl
		  << "\t freeman2sdp [options] --FreemanChain  <fileName>  "<<std::endl
		  << general_opt << "\n";
      return 0;
    }
  
  bool oneline = vm.count("oneLine");
  if( vm.count("FreemanChain") ){
    string fileName = vm["FreemanChain"].as<string>();
    vector< FreemanChain > vectFcs =  PointListReader< Z2i::Point >:: getFreemanChainsFromFile<Z2i::Integer> (fileName); 
    
    for(int i=0; i< vectFcs.size(); i++){
      bool isClosed = vectFcs.at(i).isClosed(); 
      cout << "# grid curve " << i+1 << "/" << vectFcs.size()
	   << ( (isClosed)?" closed":" open" ) << endl;
      if ( vm.count("info"))
	cout << "# SDP contour" << i+1<< "/" << vectFcs.size() << " "
	     << "# size=" << vectFcs.at(i).size() << endl;
      vector<Z2i::Point> vectPts; 
      FreemanChain::getContourPoints( vectFcs.at(i), vectPts ); 
      for(int k=0; k < vectPts.size(); k++){
	cout << vectPts.at(k)[0] << " "<< vectPts.at(k)[1];
	if(!oneline){
	  cout << endl;
	}
      }
      cout << endl;
    }
    
  }
  

}

