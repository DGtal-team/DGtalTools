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
 * @file eulerCharacteristic.cpp
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2014/06/24
 *
 *
 * This file is part of the DGtal library.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/io/readers/VolReader.h>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/IntervalForegroundPredicate.h>

using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;




/**
 @page eulerCharacteristic eulerCharacteristic
 
 @brief Computes the Euleur Characteristic of  a vol to a 8-bit raw file.

 The vol file is first binarized using interval [m,M[ thresholds and
 the Eucler characteristic is given from the cubical complex.

 @b Usage: 	 eulerCharacteristic --input <volFileName> -m <minlevel> -M <maxlevel> 


 @b Allowed @b options @b are : 
 @code
  -h [ --help ]                    display this message.
  -i [ --input ] arg               Input vol file.
  -m [ --thresholdMin ] arg (=0)   threshold min (excluded) to define binary 
                                   shape
  -M [ --thresholdMax ] arg (=255) threshold max (included) to define binary 
 @endcode

 @b Example: 

 @code
eulerCharacteristic -i $DGtal/examples/samples/cat10.vol -m 0
 @endcode


 You should obtain such a result:
@verbatim
Got 72479 cells
Got 10128 pointels 28196 linels  26112 surfels and 8043  bells
Volumetric Euler Characteristic = 1
@endverbatim
 
 @see
 @ref eulerCharacteristic.cpp

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

  // parse command line ----------------------------------------------
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>(), "Input vol file." )
    ("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min (excluded) to define binary shape" )
    ("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max (included) to define binary shape" );
  bool parseOK=true;
  po::variables_map vm;


  
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify ( vm );
  if (!parseOK || vm.count ( "help" ) ||argc<=1 )
    {
      trace.info() << "Compute the Euleur Characteristic of  a vol to a 8-bit raw file. The vol file is first binarized using interval [m,M[ thresholds and the Eucler characteristic is given from the cubical complex."<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\t eulerCharacteristic --input <volFileName> -m <minlevel> -M <maxlevel> "<<std::endl
                   << general_opt << "\n";
      return 0;
    }

  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  std::string filename = vm["input"].as<std::string>();
  int thresholdMin = vm["thresholdMin"].as<int>();
  int thresholdMax = vm["thresholdMax"].as<int>();
 
  //Importing the Vol
  trace.beginBlock("Loading the vol file");
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;
  MyImageC  imageC = VolReader< MyImageC >::importVol ( filename );
  trace.info()<<imageC<<std::endl;
  trace.endBlock();  

  //Constructing the cubical complex
  trace.beginBlock("Construting the cubical complex");
  KSpace::CellSet myCellSet;
  KSpace  ks;
  bool space_ok = ks.init( imageC.domain().lowerBound(), imageC.domain().upperBound(), true );
  if (!space_ok)
    {
      trace.error() << "Error in the Khamisky space construction."<<std::endl;
      return 2;
    }
  functors::IntervalForegroundPredicate<MyImageC> interval(imageC, thresholdMin,thresholdMax);  
  for(MyImageC::Domain::ConstIterator it =imageC.domain().begin(), itend= imageC.domain().end();
      it != itend; ++it)
    {
      if (interval( *it ))
        {
          Domain dom( 2*(*it), 2*(*it) + Point::diagonal(2));
          for(Domain::ConstIterator itdom = dom.begin(), itdomend = dom.end(); itdom != itdomend; ++itdom)
            myCellSet.insert( ks.uCell( *itdom) );
        }
    }
  trace.info() << "Got "<< myCellSet.size()<< " cells"<<std::endl;
  trace.endBlock();

  trace.beginBlock("Computing the characteristics");
  std::vector<int> cells(4,0);
  
  for(KSpace::CellSet::const_iterator it = myCellSet.begin(), itend = myCellSet.end(); it !=itend; ++it)
    cells[ ks.uDim(*it) ] ++; 
  
  trace.info() << "Got "<< cells[0]<< " pointels "<<cells[1]<<" linels  "<< cells[2]<<" surfels and "<<cells[3]<<"  bells"<<std::endl;
  trace.endBlock();

  trace.info() << "Volumetric Euler Characteristic = "<<cells[0] - cells[1] + cells[2] - cells[3]<<std::endl;

  return 0;
}
