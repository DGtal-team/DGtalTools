
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
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;


/**
 @page volTrValues volTrValues
 
 @brief  Applies basic vol image transform from the input values to output values.

 @b Usage:  	 volTrValues --input <volFileName> --o <volOutputFileName> -s 1 99 -r 100 200  

=> all voxels of values 1 (resp. 99) will be 100 (resp. 200) in the resulting image.   


 @b Allowed @b options @b are : 
 @code
  -h [ --help ]                     display this message.
  -i [ --input ] arg                Input vol file.
  -o [ --output ] arg (=output.vol) Output filename.
  -s [ --inputVals ] arg            specify the values which will be 
                                    transformed with the output values (given 
                                    with --outputVals).
  -r [ --outputVals ] arg           specify the output values to transformed 
                                    according to the input values (given with 
                                    --inputVals).
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
    ( "output,o", po::value<std::string>()->default_value("output.vol"),"Output filename." )
    ("inputVals,s", po::value<std::vector<unsigned int > >()->multitoken(), "specify the values which will be transformed with the output values (given with --outputVals)." ) 
    ("outputVals,r", po::value<std::vector<unsigned int > >()->multitoken(), "specify the output values to transformed according to the input values (given with --inputVals)." ) ;
  bool parseOK=true;
  
  po::variables_map vm;
  
  
  try{
    po::store ( po::parse_command_line ( argc, argv, general_opt ), vm );
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify ( vm );
  if (!parseOK || !vm.count("inputVals")|| !vm.count("outputVals") || !vm.count("input") || !vm.count("output")   || vm.count ( "help" ))
    {
      trace.info() << "Apply basic vol image transform from the input values to output values."<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\t volTrValues --input <volFileName> --o <volOutputFileName> -s 1 99 -r 100 200  "<<std::endl
                   << "\t => all voxel of values 1 (resp. 99) will be 100 (resp. 200) in the resulting image.   "<<std::endl
                   << general_opt << "\n";
      if( !vm.count("inputVals")){
        missingParam("inputVals");
      }
      if( !vm.count("outputVals")){
        missingParam("outputVals");
      }
      if( !vm.count("input")){
        missingParam("input");
      }
      if( !vm.count("output")){
        missingParam("output");
      }
    return 0;
    }


  //Parse options

  
  std::string filename = vm["input"].as<std::string>();
  if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
  std::string outputFileName = vm["output"].as<std::string>();


  std::vector<unsigned int> inputVals = vm["inputVals"].as<std::vector<unsigned int > >();
  std::vector<unsigned int>  outputVals = vm["outputVals"].as<std::vector<unsigned int > >();

  if(inputVals.size()!=outputVals.size()){
    trace.error()<< "Transformation not possible the two sets of input/output values should have the same size." << std::endl;
    exit(1);
  }
  
  trace.beginBlock("Loading file");
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;
  MyImageC  image = GenericReader< MyImageC >::import( filename );
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




