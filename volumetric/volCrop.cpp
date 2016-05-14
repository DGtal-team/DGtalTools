
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
 * @file volCrop.cpp
 * @ingroup volumetric/voltools
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2013/05/05
 *
 * 
 *
 * This file is part of the DGtalTools.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/VolReader.h>
#include <DGtal/io/writers/VolWriter.h>
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
 @page volCrop volCrop
 
 @brief  Crops a 3D vol image from domain coordinates.

 @b Usage: 	 volCrop --input \<volFileName\> --o \<volOutputFileName\> (both files can be independently in vol, pgm3D, p3d format)



 @b Allowed @b options @b are : 
 @code
  -h [ --help ]                     display this message.
  -i [ --input ] arg                Input vol file.
  --xMin arg (=0)                   x coordinate of lower point.
  --yMin arg (=0)                   y coordinate of lower point.
  --zMin arg (=0)                   z coordinate of lower point.
  --xMax arg                        x coordinate of upper point.
  --yMax arg                        y coordinate of upper point.
  --zMax arg                        z coordinate of upper point.
  -o [ --output ] arg (=output.vol) Output filename.
 @endcode

 @b Example: 

 @code
$ volCrop --xMin 50 --yMin 50 --zMin 10 --xMax 150 --yMax 150 --zMax 50 -i ${DGtal}/examples/samples/lobster.vol -o croppedLobster.vol 
$ 3dImageViewer -i croppedLobster.vol
 @endcode


 You should obtain such a visualization:
 @image html resVolCrop.png "Resulting visualization."
 
 @see
 @ref volCrop.cpp

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
    ( "xMin", po::value<unsigned int>()->default_value(0), "x coordinate of lower point." )
    ( "yMin", po::value<unsigned int >()->default_value(0), "y coordinate of lower point." )
    ( "zMin", po::value<unsigned int >()->default_value(0), "z coordinate of lower point." )
    ( "xMax", po::value<unsigned int >(), "x coordinate of upper point." )
    ( "yMax", po::value<unsigned int>(), "y coordinate of upper point." )
    ( "zMax", po::value<unsigned int>(), "z coordinate of upper point." )
    ( "output,o", po::value<string>()->default_value("output.vol"),"Output filename." );
  bool parseOK=true;
  
  po::variables_map vm;
  
  
  try{
    po::store ( po::parse_command_line ( argc, argv, general_opt ), vm );
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify ( vm );
  if (!parseOK || vm.count ( "help" ) || !vm.count("input"))
    {
      trace.info() << " Crops a 3D vol image from domain coordinates."<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\t volCrop --input <volFileName> --o <volOutputFileName> (both files can be independently in vol, pgm3D, p3d format)"<<std::endl
                   << general_opt << "\n";
      std::cout << "Example:\n"
		<< "volCrop --xMin 50 --yMin 50 --zMin 10 --xMax 150 --yMax 150 --zMax 50 -i ${DGtal}/examples/samples/lobster.vol -o croppedLobster.vol \n";
      return 0;
    }


  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  if ( ! ( vm.count ( "xMax" ) ) ) missingParam ( "--xMax" );
  if ( ! ( vm.count ( "yMax" ) ) ) missingParam ( "--yMax" );
  if ( ! ( vm.count ( "zMax" ) ) ) missingParam ( "--zMax" );

  std::string filename = vm["input"].as<std::string>();
  if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
  std::string outputFileName = vm["output"].as<std::string>();

  Z3i::Point ptLow( vm["xMin"].as<unsigned int>(), vm["yMin"].as<unsigned int>(),vm["zMin"].as<unsigned int>());
  Z3i::Point ptMax( vm["xMax"].as<unsigned int>(), vm["yMax"].as<unsigned int>(),vm["zMax"].as<unsigned int>());
    
  trace.beginBlock("Loading file");
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;
  MyImageC  imageC = VolReader< MyImageC >::importVol ( filename );
  functors::Identity df;  
  
  typedef ConstImageAdapter<MyImageC, Domain, functors::Identity, MyImageC::Value, functors::Identity > ConstImageAdapterForSubImage;
  Domain subDomain(ptLow, ptMax);
  ConstImageAdapterForSubImage subImage(imageC, subDomain, df, df);
  trace.endBlock();

  trace.beginBlock("Exporting...");
  bool res =  VolWriter< ConstImageAdapterForSubImage>::exportVol(outputFileName, subImage);
  trace.endBlock();
  if (res) return 0; else return 1;
}




