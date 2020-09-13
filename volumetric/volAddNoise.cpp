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
 * @file volAddNoise
 * @ingroup converters
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr)
 *
 * @date 2015/03/24
 *
 *
 * This file is part of the DGtal library.
 */

#include <DGtal/base/Common.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageSelector.h"
#include <DGtal/geometry/volumes/KanungoNoise.h>
#include <DGtal/images/IntervalForegroundPredicate.h>

#include "CLI11.hpp"

// Max component option
#include <boost/pending/disjoint_sets.hpp>

#include <vector>
#include <string>
#include <climits>

using namespace DGtal;


/**
 @page volAddNoise volAddNoise
 @brief  Adds Kanungo noise to a binary object with 0 values as background
 points and values >0 for the foreground ones.

 @b Usage:  volumetric/volAddNoise [OPTIONS] 1 [2]

 @b Allowed @b options @b are:

 @code

 Positionals:
   1 TEXT:FILE REQUIRED                  input image file name (any 3D image format accepted by DGtal::GenericReader)
   2 TEXT=result.vol                     output image file name (any 3D image format accepted by DGtal::GenericWriter)

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         input image file name (any 3D image format accepted by DGtal::GenericReader)
   -o,--output TEXT=result.vol           output image file name (any 3D image format accepted by DGtal::GenericWriter)
   -n,--noise FLOAT=0.5                  Kanungo noise level in ]0,1[ (default 0.5)
                                         
   -m,--max                              Extract only the largest 6-connected component.

  @endcode

 @b Example:
 @code
 $ volAddNoise -i $DGtal/examples/samples/Al.100.vol -o AlNoisy0.4.vol  -n 0.4
 # Converting in sdp to display:
 $ vol2sdp -i AlNoisy0.4.vol -o AlNoisy0.4.sdp
 # displaying sequence of points:
 $ 3dSDPViewer -i tmp.sdp
 @endcode

 You should obtain such a visualization:

 @image html resVolAddnoise.png "Resuling visualisation with 3dSDPViewer."

 @see volAddNoise.cpp

 */

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam( const std::string &param )
{
  trace.error() << " Parameter: " << param << " is required..";
  trace.info() << std::endl;
  exit( 1 );
}

typedef ImageSelector<Z3i::Domain, unsigned char>::Type MyImage;

int main( int argc, char ** argv )
{
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.vol"};
  double noise {0.5};
  bool MaxFlag {false};
  
  app.description("Adds Kanungo noise to a binary object with 0 values as "
  "background points and values >0 for the foreground ones.\n Basic usage:\n \t volAddNoi0se [options] --input <imageName> --output <outputImage> -noise 0.3 \n");
  app.add_option("-i,--input,1", inputFileName, "input image file name (any 3D image format accepted by DGtal::GenericReader)" )
  ->required()
  ->check(CLI::ExistingFile);
  app.add_option("-o,--output,2", outputFileName, "output image file name (any 3D image format accepted by DGtal::GenericWriter)", true);
  
  app.add_option("--noise,-n", noise, "Kanungo noise level in ]0,1[ (default 0.5)\n", true);
  app.add_flag("--max,-m", MaxFlag, "Extract only the largest 6-connected component.");
    
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

 
  typedef functors::IntervalForegroundPredicate<MyImage> Binarizer;
  MyImage image = GenericReader<MyImage>::import( inputFileName );
  trace.info() << "Input image: " << image << std::endl;
  Binarizer predicate( image, 0, 255 );

  KanungoNoise<Binarizer, Z3i::Domain> kanungo( predicate, image.domain(),
                                                noise );

  MyImage result( image.domain() );
  for ( Z3i::Domain::ConstIterator it    = image.domain().begin(),
                                   itend = image.domain().end();
        it != itend; ++it )
  {
    if ( kanungo( *it ) )
      result.setValue( *it, 255 );
    else
      result.setValue( *it, 0 );
  }

  // Exporting
  if ( !MaxFlag )
  {
    result >> outputFileName;
  }
  else
  {
    trace.beginBlock( "Extracting the largest 6-connected component" );
    typedef std::map<Z3i::Point, std::size_t> Rank; // => order on Element
    typedef std::map<Z3i::Point, Z3i::Point> Parent;
    Rank rankMap;
    Parent parentMap;
    boost::associative_property_map<Rank> rankPMap( rankMap );
    boost::associative_property_map<Parent> parentPMap( parentMap );
    boost::disjoint_sets<boost::associative_property_map<Rank>,
                         boost::associative_property_map<Parent>>
    dsets( rankPMap, parentPMap );

    trace.beginBlock( "Initial disjoint sets construction" );
    for ( auto e : result.domain() )
    {
      if ( result( e ) != 0 )
        dsets.make_set( e );
    }
    trace.endBlock();

    trace.beginBlock( "Merging neighboring sets" );
    typename Z3i::Point decx( 1, 0, 0 );
    typename Z3i::Point decy( 0, 1, 0 );
    typename Z3i::Point decz( 0, 0, 1 );

    // Merging process
    for ( auto e : result.domain() )
    {
      if ( result( e ) != 0 )
      {
        if ( result.domain().isInside( e + decx ) &&
             ( result( e ) == result( e + decx ) ) )
          dsets.union_set( e, e + decx );

        if ( result.domain().isInside( e + decy ) &&
             ( result( e ) == result( e + decy ) ) )
          dsets.union_set( e, e + decy );

        if ( result.domain().isInside( e + decz ) &&
             ( result( e ) == result( e + decz ) ) )
          dsets.union_set( e, e + decz );
      }
    }
    trace.endBlock();

    trace.beginBlock( "Counting component sizes" );
    // counting
    std::map<Z3i::Point, unsigned int> sizes;
    for ( auto p : result.domain() )
    {
      if ( result( p ) != 0 )
      {
        Z3i::Point ref = dsets.find_set( p );
        sizes[ ref ]++;
      }
    }
    unsigned int maxElement = 0;
    Z3i::Point maxP;
    for ( auto i : sizes )
    {
      if ( maxElement < i.second )
      {
        maxElement = i.second;
        maxP       = i.first;
      }
    }
    trace.info() << "Largest component has " << maxElement << " voxels."
                 << std::endl;
    trace.endBlock();

    trace.beginBlock( "Cleaning up" );
    // Cleaning
    // Merging process
    Z3i::Point largest = dsets.find_set( maxP );
    trace.info() << "Largest ref point: " << largest << std::endl;
    for ( auto e : result.domain() )
    {
      if ( result( e ) != 0 )
      {
        if ( dsets.find_set( e ) != largest )
          result.setValue( e, 0 );
      }
    }
    trace.endBlock();
    trace.endBlock();
    result >> outputFileName;
  }
  return 0;
}
