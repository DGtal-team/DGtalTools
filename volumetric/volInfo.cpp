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
 * @file volInfo.cpp
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 * @ingroup Volumetric 
 * @date 2022/02/11
 *
 *
 * This file is part of the DGtal library.
 */

#include <iostream>
#include <map>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;
using namespace Z3i;

/**
 @page volInfo volInfo
 
 @brief Get information from a vol file (size and values).
 @ingroup volumetrictools

 @b Usage: 	./volumetric/volInfo [OPTIONS] 1 [2]


 @b Allowed @b options @b are : 
 @code
 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         Input vol file.
  
 Positionals:
   1 TEXT:FILE REQUIRED                  Input vol file.
 @endcode

 @b Example: 
You can retrieve information from the vol file as:
 @code
 $ volInfo $DGtal/examples/samples/lobster.vol
 @endcode

 @see
 @ref volInfo.cpp

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
    // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
    
  app.description("Retreive information from vol file\n Basic usage: \n \tvolInfo  <volFileName>  ");
  app.add_option("-i,--input,1", inputFileName, "Input vol file." )->required()->check(CLI::ExistingFile);
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  trace.beginBlock("Loading file");
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;

  MyImageC  imageC = GenericReader<MyImageC>::import(inputFileName);
  
  trace.info()<<imageC<<std::endl;
  trace.info()<<"Scanning the values:"<<std::endl;

  std::map<unsigned char, size_t> values;
  for(auto v: imageC.range())
    values[ v ] ++;
  
  for(auto val: values)
    std::cout<<"\tvalue "<< (int)val.first<<": "<< val.second<<" voxels."<<std::endl;
  
  trace.endBlock();
  return 0;
}
