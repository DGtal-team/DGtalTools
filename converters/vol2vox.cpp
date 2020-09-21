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
 * @file vol2vox.cpp
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2018/01/25
 *
 *
 * This file is part of the DGtal library.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/VolReader.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
 @page vol2vox
 @brief  Converts a vol file to a MagicaVoxel VOX file (https://ephtracy.github.io).


 @b Usage: vo2vox -i [input] -o [output]

 @b Allowed @b options @b are:

 @code

Positionals:
  2 TEXT=result.vox                     Output filename.

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         Input vol file.
  -o,--output TEXT=result.vox           Output filename.
 @endcode

 @b Example:
 @code
 $ vol2vox -i ${DGtal}/examples/samples/Al.100.vol -o Al.100.vox

 @endcode
 */



/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam ( const std::string &param )
{
  trace.error() <<" Parameter: "<<param<<" is required..";
  trace.info() <<std::endl;
  exit ( 1 );
}

template <typename Word>
static
std::ostream& write_word( std::ostream& outs, Word value )
{
  for (unsigned size = sizeof( Word ); size; --size, value >>= 8)
    outs.put( static_cast <char> (value & 0xFF) );
  return outs;
}


int main(int argc, char**argv)
{
  

  
   // parse command line using CLI ----------------------------------------------
   CLI::App app;
   std::string inputFileName;
   std::string outputFileName{"result.vox"};
  app.description("Convert a vol (up to [0,126]^3) to a vox.\n The current exporter does not support custom colormaps, only the voxel values are exported \n Example: vol2vox  -i ${DGtal}/examples/samples/Al.100.vol -o Al.100.vox");
  app.add_option("-i,--input", inputFileName, "Input vol file.")
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("-o,--output,2", outputFileName, "Output filename.", true);
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);

  
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;
  trace.beginBlock("Loading..");
  MyImageC  imageL = VolReader< MyImageC >::importVol ( inputFileName );
  trace.endBlock();
 
   if ( (imageL.domain().upperBound() - imageL.domain().lowerBound()).max() > 126 )
   {
     trace.error() << "Vol file too large (width > 126)."<<std::endl;
     trace.info() << std::endl;
     exit(2);
   }
 
   trace.beginBlock("Exporting...");
   ofstream myfile;
   myfile.open (outputFileName, ios::out | ios::binary);

   DGtal::uint32_t cpt=0;
   for(auto it = imageL.range().begin(), itend = imageL.range().end();
       it!=itend; ++it)
     if (*it != 0)
       cpt++;

   Point size = imageL.domain().upperBound() - imageL.domain().lowerBound();

   /*
    * VOX file format:
    *
    4b VOX' '
    4b version (150)

    4b MAIN (chunckid)
    4b size chucnk content (n)
    4b size chunck children (m)

    4b SIZE (chunckid)
    4b chunck content
    4b chunck children
    4bx3  x,y,z

    4b VOXEL (chunckid)
    4b chunck content
    4b chunck children
    4b number of voxels
    1b x 4 (x,y,z,idcol)  x numVoxels

    */

   trace.info()<<size<<std::endl;

   //HEADER
   myfile <<'V'<<'O'<<'X'<<' ';
   //  version 150
   write_word(myfile, (DGtal::uint32_t)150);

   //Chunck MAIN
   myfile <<'M'<<'A'<<'I'<<'N';
   write_word(myfile,DGtal::uint32_t(0)); //size content
   write_word(myfile,DGtal::uint32_t( (4+4+4 + 12 ) + ((4+ 4 + 4) + (4+ cpt*4)))); //size children

   //Chunck SIZE
   myfile <<'S'<<'I'<<'Z'<<'E';
   write_word(myfile,DGtal::uint32_t(12)); //3x4
   write_word(myfile,DGtal::uint32_t(0));  //0 children
   write_word(myfile,DGtal::uint32_t(size[0]+1));
   write_word(myfile,DGtal::uint32_t(size[1]+1));
   write_word(myfile,DGtal::uint32_t(size[2]+1));

   //Chunck VOXEL
   myfile << 'X'<<'Y'<<'Z'<<'I';
   write_word(myfile, (DGtal::uint32_t)(4+cpt*4));  // 4 + numvoxel * 4
   write_word(myfile,DGtal::uint32_t(0));  //0 children
   write_word(myfile, DGtal::uint32_t(cpt)); //numvoxels

   trace.info() << "Number of voxels= "<<cpt<<std::endl;

   //Data
   for(auto it = imageL.domain().begin(), itend = imageL.domain().end();
       it!=itend; ++it)
     if (imageL(*it) != 0)
     {
       Point p = (*it) - imageL.domain().lowerBound();
       myfile.put((DGtal::uint8_t)p[0]);
       myfile.put( (DGtal::uint8_t)p[1]);
       myfile.put( (DGtal::uint8_t)p[2]);
       myfile.put(imageL(*it));
     }

   myfile.close();
   trace.endBlock();

  return 0;
}
