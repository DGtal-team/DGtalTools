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

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;



/**
 @page volSubSample volSubSample
 
 @brief Brutally sub samples a vol file (division by 2 in each direction).


 @b Usage: 	volSubSample --input <volFileName> --o <volOutputFileName>


 @b Allowed @b options @b are : 
 @code
  -h [ --help ]                 display this message.
  -i [ --input ] arg            Input vol file.
  -o [ --output ] arg           Output filename.
  -f [ --function ] arg (=mean) Function used to the down-sampling: {none,max, min, mean}
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
Val maxVal(Image &image, Point &p, Domain& domain)
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
Val minVal(Image &image, Point &p, Domain& domain)
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
Val meanVal(Image &image, Point &p, Domain& domain)
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
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>(), "Input vol file." )
    ( "output,o", po::value<string>(),"Output filename." )
    ("function,f",   po::value<string>()->default_value("mean"), "Function used to the down-sampling: {none,max, min, mean}" );
  

  po::variables_map vm;
  po::store ( po::parse_command_line ( argc, argv, general_opt ), vm );
  po::notify ( vm );
  if ( vm.count ( "help" ) ||argc<=1 )
    {
      trace.info() << "Brutally sub sample a vol file (division by 2 in each direction)."<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\tvolSubSample --input <volFileName> --o <volOutputFileName> "<<std::endl
                   << general_opt << "\n";
      return 0;
    }

  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  std::string filename = vm["input"].as<std::string>();
  std::string function = vm["function"].as<std::string>();
  if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
  std::string outputFileName = vm["output"].as<std::string>();
  
  
  trace.beginBlock("Loading file");
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;

  MyImageC  imageC = VolReader< MyImageC >::importVol ( filename );
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
