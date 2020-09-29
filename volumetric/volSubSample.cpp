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
 * @file voAddBorder.cpp
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/05/01
 *
 *
 * This file is part of the DGtal library.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/VolReader.h>
#include <DGtal/io/writers/VolWriter.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
 @page volSubSample volSubSample
 
 @brief Brutally sub samples a vol file (division by 2 in each direction).


 @b Usage: 	./volumetric/volSubSample [OPTIONS] 1 [2]


 @b Allowed @b options @b are : 
 @code
 Positionals:
   1 TEXT:FILE REQUIRED                  Input vol file.
   2 TEXT=result.vol                     Output filename.

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         Input vol file.
   -o,--output TEXT=result.vol           Output filename.
   -f,--function TEXT:{mean,none,max,min,mean}=mean
                                         Function used to the down-sampling: {none,max, min, mean}

 Positionals:
   1 TEXT:FILE REQUIRED                  Input vol file.
   2 TEXT=result.vol                     Output filename.

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         Input vol file.
   -o,--output TEXT=result.vol           Output filename.
   -f,--function TEXT:{mean,none,max,min,mean}=mean
                                         Function used to the down-sampling: {none,max, min, mean}

 @endcode

 @b Example: 
You can apply several sub sampling:
 @code
 $ volSubSample -i $DGtal/examples/samples/lobster.vol -o lobster2.vol -f mean 
 $ volSubSample -i lobster2.vol -o lobster4.vol -f mean 
 $ volSubSample -i lobster4.vol -o lobster8.vol -f mean 
 @endcode

You can  display the result by extracting the surface using \ref 3dVolMarchingCubes:

@code 
$ 3dVolMarchingCubes -i $DGtal/examples/samples/lobster.vol -o lobster.off -t 30
$ 3dVolMarchingCubes -i lobster2.vol -t 30 -o lobster2.off
$ 3dVolMarchingCubes -i lobster4.vol -t 30 -o lobster4.off
$ 3dVolMarchingCubes -i lobster8.vol -t 30 -o lobster8.off
$ meshViewer -i lobster.off lobster2.off  lobster4.off  lobster8.off -n
@endcode

 You should obtain such a result:
 @image html  resVolSubSample.png "Resulting visualization."
 
 @see
 @ref volSubSample.cpp

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

template<typename Val, typename Image, typename Point, typename Domain>
Val maxVal(Image const& image, Point const& p, Domain const& domain)
{
  typename Image::Domain dom( p*2, p*2 + Point::diagonal(1));
  Val v=image(p*2);
  for(typename Image::Domain::ConstIterator it=dom.begin(), itend=dom.end();
      it != itend;
      ++it)
    if ( domain.isInside(*it) &&  image( *it) > v) v=image(*it);
  
  return v;    
} 
template<typename Val, typename Image, typename Point, typename Domain>
Val minVal(Image const& image, Point const& p, Domain const& domain)
{
  typename Image::Domain dom( p*2, p*2 + Point::diagonal(1));
  Val v=image(p*2);
  for(typename Image::Domain::ConstIterator it=dom.begin(), itend=dom.end();
      it != itend;
      ++it)
    if (  domain.isInside(*it)&&  image( *it) < v) v=image(*it);
  
  return v;    
} 
template<typename Val, typename Image, typename Point, typename Domain>
Val meanVal(Image const& image, Point const& p, Domain const& domain)
{
  typename Image::Domain dom( p*2, p*2 + Point::diagonal(1));
  int v=0;
  int nb=0;
  for(typename Image::Domain::ConstIterator it=dom.begin(), itend=dom.end();
      it != itend;
      ++it)
    if ( domain.isInside(*it) )
      {
	nb++;
	v+=image(*it);
      }
  return static_cast<unsigned char>( v/nb );
} 


int main(int argc, char**argv)
{
  
  // parse command line ----------------------------------------------
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.vol"};
  std::string function {"mean"};
  
  app.description("Brutally sub sample a vol file (division by 2 in each direction).\n Basic usage: \n \tvolSubSample --input <volFileName> --o <volOutputFileName> ");

  
  app.add_option("-i,--input,1", inputFileName, "Input vol file." )
  ->required()
  ->check(CLI::ExistingFile);
  app.add_option("-o,--output,2", outputFileName, "Output filename.",true);
  app.add_option("-f,--function", function, "Function used to the down-sampling: {none,max, min, mean}", true)
   -> check(CLI::IsMember({"mean", "none", "max", "min", "mean"}));
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  trace.beginBlock("Loading file");
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;

  MyImageC  imageC = VolReader< MyImageC >::importVol ( inputFileName );
  MyImageC  outputImage( Z3i::Domain( imageC.domain().lowerBound(),
                                      (imageC.domain().upperBound()-imageC.domain().lowerBound())/Vector().diagonal(2)));

  trace.endBlock();
  Point subvector = Vector().diagonal(2);
  Point p;
  unsigned char val;
  trace.beginBlock("Down-scaling the volume...");
  trace.info()<<"Function= "<<function<<std::endl;
  trace.info() << outputImage.domain() << std::endl;
  //Fast Copy
  if (function == "none")
    for(MyImageC::Domain::ConstIterator it = outputImage.domain().begin(),
	  itend = outputImage.domain().end(); it != itend; ++it)
      {
	p = (*it) * 2;
	outputImage.setValue( *it , imageC( p ));
      }
  else
    if (function == "max")
      for(MyImageC::Domain::ConstIterator it = outputImage.domain().begin(),
	    itend = outputImage.domain().end(); it != itend; ++it)
	{
	  val = maxVal<unsigned char, MyImageC, Point>(imageC, *it, imageC.domain());
	  outputImage.setValue( *it , val );
	}
    else
      if (function == "min")
	for(MyImageC::Domain::ConstIterator it = outputImage.domain().begin(),
	      itend = outputImage.domain().end(); it != itend; ++it)
	  {
	    val = minVal<unsigned char, MyImageC, Point>(imageC, *it, imageC.domain());
	    outputImage.setValue( *it , val);
	  }  
      else
	
      if (function == "mean")
	for(MyImageC::Domain::ConstIterator it = outputImage.domain().begin(),
	      itend = outputImage.domain().end(); it != itend; ++it)
	  {
	    val = meanVal<unsigned char, MyImageC, Point>(imageC, *it, imageC.domain());
	    outputImage.setValue( *it , val );
	  }
      else
	trace.error() << "Bad function !"<<std::endl;


  trace.endBlock();
  
  trace.beginBlock("Exporting...");
  bool res =  VolWriter< MyImageC>::exportVol(outputFileName, outputImage);
  trace.endBlock();
  if (res) return 0; else return 1;
}
