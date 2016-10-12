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
 * @file at-u0-v1.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2016/10/12
 *
 * A tool file named at-u0-v1.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include <sstream>
#include <string>
#include <boost/format.hpp>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

#include "ATu0v1.h"

using namespace std;
using namespace DGtal;

int main( int argc, char* argv[] )
{
  using namespace Z2i;

  // parse command line ----------------------------------------------
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<string>(), "the input image PPM filename." )
    ("output,o", po::value<string>()->default_value( "AT" ), "the output image basename." )
    ("lambda,l", po::value<double>(), "the parameter lambda." )
    ("lambda-1,1", po::value<double>()->default_value( 0.3125 ), "the initial parameter lambda (l1)." ) // 0.3125
    ("lambda-2,2", po::value<double>()->default_value( 0.00005 ), "the final parameter lambda (l2)." )
    ("lambda-ratio,q", po::value<double>()->default_value( sqrt(2) ), "the division ratio for lambda from l1 to l2." )
    ("alpha,a", po::value<double>()->default_value( 1.0 ), "the parameter alpha." )
    ("epsilon,e", po::value<double>()->default_value( 1.0 ), "the initial and final parameter epsilon of AT functional at the same time." )
    ("epsilon-1", po::value<double>(), "the initial parameter epsilon." )
    ("epsilon-2", po::value<double>(), "the final parameter epsilon." )
    ("epsilon-r", po::value<double>()->default_value( 2.0 ), "sets the ratio between two consecutive epsilon values of AT functional." )
    ("nbiter,n", po::value<int>()->default_value( 10 ), "the maximum number of iterations." )
    ("snr", "force computation of SNR." )
    ("image-snr", po::value<string>(), "the input image without deterioration." )
    ;

  bool parseOK=true;
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  } catch ( const exception& ex ) {
    parseOK = false;
    cerr << "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);
  if ( ! parseOK || vm.count("help")
                 || !vm.count("input")
                 || (vm.count("snr") && !vm.count("image-snr"))
     )
    {
      cerr << "Usage: " << argv[0] << " -i toto.pgm\n"
       << "Computes the Ambrosio-Tortorelli reconstruction/segmentation of an input image."
       << endl << endl
       << " / "
       << endl
       << " | a.(u-g)^2 + v^2 |grad u|^2 + le.|grad v|^2 + (l/4e).(1-v)^2 "
       << endl
       << " / "
       << endl
       << "Discretized as (u 0-form, v 1-form, A vertex-edge bdry, B edge-face bdy)" << endl
       << "E(u,v) = a(u-g)^t (u-g) +  u^t A^t diag(v)^2 A^t u + l e v^t (A A^t + B^t B) v + l/(4e) (1-v)^t (1-v)" << endl
       << endl
       << general_opt << "\n"
       << "Example: ./at-u0-v1 -i ../Images/cerclesTriangle64b02.pgm -o tmp -a 0.05 -e 1 --lambda-1 0.1 --lambda-2 0.00001 -g"
       << endl;
      return 1;
    }
  string f1 = vm[ "input" ].as<string>();
  string f2 = vm[ "output" ].as<string>();
  double l1  = vm[ "lambda-1" ].as<double>();
  double l2  = vm[ "lambda-2" ].as<double>();
  double lr  = vm[ "lambda-ratio" ].as<double>();
  if ( vm.count( "lambda" ) ) l1 = l2 = vm[ "lambda" ].as<double>();
  if ( l2 > l1 ) l2 = l1;
  if ( lr <= 1.0 ) lr = sqrt(2);
  double a  = vm[ "alpha" ].as<double>();
  double e  = vm[ "epsilon" ].as<double>();
  double e1 = vm.count( "epsilon-1" ) ? vm[ "epsilon-1" ].as<double>() : e;
  double e2 = vm.count( "epsilon-2" ) ? vm[ "epsilon-2" ].as<double>() : e;
  double er = vm[ "epsilon-r" ].as<double>();

  bool color_image = f1.size() > 4 && f1.compare( f1.size() - 4, 4, ".ppm" ) == 0;
  bool grey_image  = f1.size() > 4 && f1.compare( f1.size() - 4, 4, ".pgm" ) == 0;
  if ( ! color_image && ! grey_image ) 
    {
      trace.error() << "Input image file must be either a PGM (grey-level) or a PPM (color) image with these extensions."
                    << endl;
      return 2;
    }

  KSpace K;
  ATu0v1< KSpace > AT;

  typedef ImageContainerBySTLVector<Domain, Color> ColorImage;
  typedef ImageContainerBySTLVector<Domain, unsigned char> GreyLevelImage;
  //---------------------------------------------------------------------------
  if ( color_image ) 
    {
      trace.beginBlock("Reading PPM image");
      ColorImage image = PPMReader<ColorImage>::importPPM( f1 );
      K.init( image.domain().lowerBound(), image.domain().upperBound(), true );
      AT.init( K );
      AT.addInput( image, [] ( Color c ) { return (double) c.red()   / 255.0; } );
      AT.addInput( image, [] ( Color c ) { return (double) c.green() / 255.0; } );
      AT.addInput( image, [] ( Color c ) { return (double) c.blue()  / 255.0; } );
      trace.endBlock();
    }
  else if ( grey_image ) 
    {
      trace.beginBlock("Reading PGM image");
      GreyLevelImage image = PGMReader<GreyLevelImage>::importPGM( f1 );
      K.init( image.domain().lowerBound(), image.domain().upperBound(), true );
      AT.init( K );
      AT.addInput( image, [] (unsigned char c ) { return (double) c / 255.0; } );
      trace.endBlock();
    }
  //---------------------------------------------------------------------------
  AT.setAlpha( a );
  AT.setUFromInput();
  trace.info() << AT << std::endl;
  while ( l1 >= l2 )
    {
      trace.info() << "************ lambda = " << l1 << " **************" << endl;
      AT.setLambda( l1 );
      l1 /= lr;
    }
  return 0;
}
