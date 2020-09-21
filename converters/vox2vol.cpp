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
 * @file vox2vol.cpp
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
#include <DGtal/io/writers/GenericWriter.h>

#include "CLI11.hpp"



using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
   @page vox2vol
   @brief  Converts a MagicaVoxel VOX file (https://ephtracy.github.io) to a vol file.


   @b Usage: vox2vol -i [input] -o [output]

   @b Allowed @b options @b are:

   @code

   Positionals:
   1 TEXT:FILE REQUIRED                  
   2 TEXT=result.vol                     
   
   Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         
   -o,--ouput TEXT=result.vol            
   @endcode

   @b Example:
   @code
   $ vox2vol -i Al.100.vox -o Al.100.vol

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

/**
 * Generic read word (binary mode) in little-endian mode.
 *
 * @param fin input stream.
 * @param aValue value to write.
 *
 * @return modified stream.
 */
template <typename Word>
static
std::istream& read_word( std::istream&  fin, Word& aValue )
{
  aValue = 0;
  char c;
  for (auto size = 0; size < sizeof( Word ); ++size)
  {
    fin.get( c ) ;
    unsigned char cc=static_cast<unsigned char>(c);
    aValue |= (cc << (8 * size));
  }
  return fin;
}
template <typename Word>
static
std::ostream& write_word( std::ostream& outs, Word value )
{
  for (unsigned size = sizeof( Word ); size; --size, value >>= 8)
    outs.put( static_cast <char> (value & 0xFF) );
  return outs;
}

DGtal::uint32_t toInt(const char a,
                      const char b,
                      const char c,
                      const char d)
{
  return (static_cast<DGtal::uint32_t>((unsigned char)a) +
          ((static_cast<DGtal::uint32_t>((unsigned char) b)) <<8) +
          ((static_cast<DGtal::uint32_t>((unsigned char) c)) <<16) +
          ((static_cast<DGtal::uint32_t>((unsigned char) d)) <<24));
}


int main(int argc, char**argv)
{

  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.vol"};


  app.description("Convert a vox file to a vol. Basic usage:\n vox2vol --input <volFileName> --o <volOutputFileName> ");
  app.add_option("-i,--input,1", inputFileName, "" )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("-o,--ouput,2", outputFileName, "", true );   

   
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  trace.beginBlock("Loading...");
  ifstream myfile;
  myfile.open (inputFileName, ios::in | ios::binary);
   
  /* VOX file format:
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

  //HEADER
  char a,b,c,d;
  myfile.get(a);
  myfile.get(b);
  myfile.get(c);
  myfile.get(d);


  if ( ( a != 'V' || (b != 'O') || (c!='X') || (d != ' ')))
  {
    trace.error() << "Magic number error"<<std::endl;
    trace.error() << (int)a<<" "<<(int)b<<" "<<(int) c<<" "<<(int)d<<std::endl;
    trace.error() << a<<" "<<b<<" "<< c<<" "<<d<<std::endl;

    exit(2);
  }

  myfile.get(a);
  myfile.get(b);
  myfile.get(c);
  myfile.get(d);
  DGtal::uint32_t version = toInt(a,b,c,d);
  trace.info()<<"Version = "<<version<<std::endl;
  //read_word(myfile, version);
  if (version != 150)
  {
    trace.error() << "Version error   "<<version<<std::endl;
    trace.error() << (unsigned int)a<<" "<<(int)b<<" "<<(int) c<<" "<<(int)d<<std::endl;
    trace.error() << a<<" "<<b<<" "<< c<<" "<<d<<std::endl;
    exit(2);
  }

  read_word(myfile, version);
  DGtal::uint32_t main = toInt('M','A','I','N');
  trace.info()<< main << std::endl;
  trace.info() <<version <<std::endl;
  if ( version != main)
  {
    trace.error() << "MAIN number error"<<std::endl;
    trace.error() << (int)a<<" "<<(int)b<<" "<<(int) c<<" "<<(int)d<<std::endl;
    trace.error() << a<<" "<<b<<" "<< c<<" "<<d<<std::endl;
    exit(2);
  }

  DGtal::uint32_t XYZI= toInt('X','Y','Z','I');
  read_word(myfile,version);
  while ( version != XYZI)
    read_word(myfile,version);

  //trash two ints
  read_word(myfile,version);
  read_word(myfile,version);
  DGtal::uint32_t cpt;
  read_word(myfile,cpt);

  Z3i::Domain domain(Z3i::Point(0,0,0), Z3i::Point(126,126,126));
  ImageContainerBySTLVector<Z3i::Domain, unsigned char>  image(domain);
  trace.info()<< "Number of voxels in this chunk = "<<version<<std::endl;
  for(auto i=0 ; i<cpt; ++i)
  {
    myfile.get(a);
    myfile.get(b);
    myfile.get(c);
    myfile.get(d);
    image.setValue(Z3i::Point((unsigned int)(unsigned char)a,
                              (unsigned int)(unsigned char)b,
                              (unsigned int)(unsigned char)c),
                   (unsigned char) d);
  }

  image >> outputFileName;

  return EXIT_SUCCESS;

}
