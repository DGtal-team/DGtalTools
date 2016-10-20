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
 * @file at-u2-v0.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 * @author Marion Foare (\c marion.foare@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2016/10/12
 *
 * A tool file named at-u2-v0.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include <sstream>
#include <string>
#include <functional>
#include <boost/format.hpp>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

#include "ATu2v0.h"

/**
@page DocATu2v0 imageProcessing/at-u2-v0 

@brief Computes a piecewise smooth approximation of an image, by optimizing the Ambrosio-Tortorelli functional (with u a 2-form and v a 0-form).

@writers Marion Foare, Jacques-Olivier Lachaud

@b Usage: at-u2-v0 -i [input.pgm]

@b Usage: at-u2-v0 -i [input.ppm]

Computes the Ambrosio-Tortorelli reconstruction/segmentation of an input image, either grey-level (.pgm) or color image (.ppm).

\f$ AT_e = \int a.(u-g)^2 + v^2 | \nabla u|^2 + le.| \nabla v|^2 + (l/4e).(1-v)^2 \f$
 
Discretized as (u 0-form, v 1-form, A vertex-edge bdry, B edge-face bdy)

\f$ E(u,v) = a(u-g)^t (u-g) +  u^t A^t diag(v)^2 A^t u + l e v^t (A A^t + B^t B) v + l/(4e) (1-v)^t (1-v) \f$

\b Allowed \b options \b are:

\code
  -h [ --help ]                         display this message
  -i [ --input ] arg                    the input image PPM filename.
  -m [ --inpainting-mask ] arg          the input inpainting mask filename.
  -o [ --output ] arg (=AT)             the output image basename.
  -l [ --lambda ] arg                   the parameter lambda.
  -1 [ --lambda-1 ] arg (=0.3125)       the initial parameter lambda (l1).
  -2 [ --lambda-2 ] arg (=0.00050000000000000001)
                                        the final parameter lambda (l2).
  -q [ --lambda-ratio ] arg (=1.4142135623730951)
                                        the division ratio for lambda from l1 
                                        to l2.
  -a [ --alpha ] arg (=1)               the parameter alpha.
  -e [ --epsilon ] arg                  the initial and final parameter epsilon
                                        of AT functional at the same time.
  --epsilon-1 arg (=2)                  the initial parameter epsilon.
  --epsilon-2 arg (=0.25)               the final parameter epsilon.
  --epsilon-r arg (=2)                  sets the ratio between two consecutive 
                                        epsilon values of AT functional.
  -n [ --nbiter ] arg (=10)             the maximum number of iterations.
  --image-snr arg                       the input image without deterioration 
                                        if you wish to compute the SNR.
  -p [ --pixel-size ] arg (=1)          the pixel size for outputing images 
                                        (useful when one wants to see the 
                                        discontinuities v on top of u).
  -c [ --color-v ] arg (=0xff0000)      the color chosen for displaying the 
                                        singularities v (e.g. red is 0xff0000).
  -v [ --verbose ] arg (=0)             the verbose level (0: silent, 1: less 
                                        silent, etc).
\endcode

@image html resATu0v1-cb2-a1_00000-l1_0000000-u.png "AT alpha=1 lambda=1 on carre noise=0.2"

@b example:

\code
./at-u2-v0 -i ../Images/cerclesTriangle64b02.pgm -o AT -a 0.05 -e 1 --lambda-1 0.1 --lambda-2 0.00001
\endcode


<center>
<table>
<tr>
<td colspan=5>
epsilon scale space
</td>
</tr>
<tr>
  <td> <img height=100px src="resATu0v1-cb2-e2_0-a0_10000-l0_0062092-u0-v1.png"/> </td>
  <td> <img height=100px src="resATu0v1-cb2-e1_0-a0_10000-l0_0062092-u0-v1.png"/> </td>
  <td> <img height=100px src="resATu0v1-cb2-e0_5-a0_10000-l0_0062092-u0-v1.png"/> </td>
  <td> <img height=100px src="resATu0v1-cb2-e0_25-a0_10000-l0_0062092-u0-v1.png"/> </td>
</tr>
<tr>
      <td align = center rowspan="4"> ./build/at-u2-v0 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.1 --lambda 0.006 --epsilon-1 2.0 --epsilon-2 2.0</td>
</tr>
<tr>
      <td align = center rowspan="4"> ./build/at-u2-v0 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.1 --lambda 0.006 --epsilon-1 2.0 --epsilon-2 1.0</td>
</tr>
<tr>
      <td align = center rowspan="4"> ./build/at-u2-v0 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.1 --lambda 0.006 --epsilon-1 2.0 --epsilon-2 0.5</td>
</tr>
<tr>
      <td align = center rowspan="4"> ./build/at-u2-v0 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.1 --lambda 0.006 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
<tr>
<td colspan=5>
alpha scale space
</td>
</tr>
<tr>
    <td> <img height=200px src="resATu0v1-cb2-a1_00000-l1_0000000-u.png"/> </td>
    <td> <img height=200px src="resATu0v1-cb2-a0_50000-l1_0000000-u.png"/> </td>
    <td> <img height=200px src="resATu0v1-cb2-a0_10000-l1_0000000-u.png"/> </td>
    <td> <img height=200px src="resATu0v1-cb2-a0_05000-l1_0000000-u.png"/> </td>
    <td> <img height=200px src="resATu0v1-cb2-a0_01000-l1_0000000-u.png"/> </td>
</tr>
<tr>
        <td align = center rowspan="5"> ./build/at-u2-v0 -i Images/carre2Degradesb02.pgm -o cb2 -a 1.0 --lambda 1.0 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
<tr>
        <td align = center rowspan="5"> ./build/at-u2-v0 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.5 --lambda 1.0 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
<tr>
        <td align = center rowspan="5"> ./build/at-u2-v0 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.1 --lambda 1.0 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
<tr>
        <td align = center rowspan="5"> ./build/at-u2-v0 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.05 --lambda 1.0 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
<tr>
        <td align = center rowspan="5"> ./build/at-u2-v0 -i Images/carre2Degradesb02.pgm -o cb2 -a 0.01 --lambda 1.0 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
<tr>
<td colspan=5>
lambda scale space (lena)
</td>
</tr>
<tr>
<td> <img height=200px src="resATu0v1-lena-370-b02-a0_48000-l0_2000000-u0-v1.png"/> </td>
<td> <img height=200px src="resATu0v1-lena-370-b02-a0_48000-l0_1000000-u0-v1.png"/> </td>
<td> <img height=200px src="resATu0v1-lena-370-b02-a0_48000-l0_0500000-u0-v1.png"/> </td>
<td> <img height=200px src="resATu0v1-lena-370-b02-a0_48000-l0_0250000-u0-v1.png"/> </td>
<td> <img height=200px src="resATu0v1-lena-370-b02-a0_48000-l0_0125000-u0-v1.png"/> </td>
</tr>
<tr>
        <td align = center rowspan="5"> ./build/at-u2-v0 -i Images/lena-370-b02.ppm -o lena -a 0.48 --lambda-1 0.15 --lambda-2 0.0125 -- lambda-ratio 2.0 --epsilon-1 2.0 --epsilon-2 0.25</td>
</tr>
</table>
</center>

*/


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
    ("inpainting-mask,m", po::value<string>(), "the input inpainting mask filename." )
    ("output,o", po::value<string>()->default_value( "AT" ), "the output image basename." )
    ("metric-average,M", "use metric average to smooth L1-metric." )
    ("lambda,l", po::value<double>(), "the parameter lambda." )
    ("lambda-1,1", po::value<double>()->default_value( 0.3125 ), "the initial parameter lambda (l1)." )
    ("lambda-2,2", po::value<double>()->default_value( 0.0005 ), "the final parameter lambda (l2)." )
    ("lambda-ratio,q", po::value<double>()->default_value( sqrt(2) ), "the division ratio for lambda from l1 to l2." )
    ("alpha,a", po::value<double>()->default_value( 1.0 ), "the parameter alpha." )
    ("epsilon,e", po::value<double>(), "the initial and final parameter epsilon of AT functional at the same time." )
    ("epsilon-1", po::value<double>()->default_value( 2.0 ), "the initial parameter epsilon." )
    ("epsilon-2", po::value<double>()->default_value( 0.25 ), "the final parameter epsilon." )
    ("epsilon-r", po::value<double>()->default_value( 2.0 ), "sets the ratio between two consecutive epsilon values of AT functional." )
    ("nbiter,n", po::value<int>()->default_value( 10 ), "the maximum number of iterations." )
    ("image-snr", po::value<string>(), "the input image without deterioration if you wish to compute the SNR." )
    ("pixel-size,p", po::value<int>()->default_value( 1 ), "the pixel size for outputing images (useful when one wants to see the discontinuities v on top of u)." )
    ("color-v,c", po::value<string>()->default_value( "0xff0000" ), "the color chosen for displaying the singularities v (e.g. red is 0xff0000)." )
    ("verbose,v", po::value<int>()->default_value( 0 ), "the verbose level (0: silent, 1: less silent, etc)." )
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
  if ( ! parseOK || vm.count("help") || !vm.count("input") )
    {
      cerr << "Usage: " << argv[0] << " -i toto.pgm\n"
           << "Computes the Ambrosio-Tortorelli reconstruction/segmentation of an input image."
           << "It outputs 2 or 3 images (of basename given by option --output) giving the"
           << " reconstructed image u, and other images superposing u and the discontinuities v."
           << endl << endl
           << " / "
           << endl
           << " | a.(u-g)^2 + v^2 |grad u|^2 + le.|grad v|^2 + (l/4e).(1-v)^2 "
           << endl
           << " / "
           << endl
           << "Discretized as (u 2-form, v 0-form, A vertex-edge bdry, B edge-face bdy, M vertex-edge average)" << endl
           << "E(u,v) = a(u-g)^t (u-g) +  u^t B diag(M v)^2 B^t u + l e v^t A^t A v + l/(4e) (1-v)^t (1-v)" << endl
           << endl
           << general_opt << "\n"
           << "Example: ./at-u2-v0 -i ../Images/cerclesTriangle64b02.pgm -o tmp -a 0.05 -e 1 --lambda-1 0.1 --lambda-2 0.00001 -g"
           << endl;
      return 1;
    }
  string f1  = vm[ "input" ].as<string>();
  string f2  = vm[ "output" ].as<string>();
  bool metric= vm.count( "metric-average" );
  double l1  = vm[ "lambda-1" ].as<double>();
  double l2  = vm[ "lambda-2" ].as<double>();
  double lr  = vm[ "lambda-ratio" ].as<double>();
  if ( vm.count( "lambda" ) ) l1 = l2 = vm[ "lambda" ].as<double>();
  if ( l2 > l1 ) l2 = l1;
  if ( lr <= 1.0 ) lr = sqrt(2);
  double a   = vm[ "alpha" ].as<double>();
  double e1  = vm[ "epsilon-1" ].as<double>();
  double e2  = vm[ "epsilon-2" ].as<double>();
  if ( vm.count( "epsilon" ) )
    e1 = e2 =  vm[ "epsilon" ].as<double>();
  double er  = vm[ "epsilon-r" ].as<double>();
  int  verb  = vm[ "verbose" ].as<int>();
  int nbiter = vm[ "nbiter" ].as<int>();
  int pix_sz = vm[ "pixel-size" ].as<int>();
  string scv = vm[ "color-v" ].as<string>();
  bool snr   = vm.count( "image-snr" );
  string isnr= snr ? vm[ "image-snr" ].as<string>() : "";
  Color color_v( (unsigned int) std::stoul( scv, nullptr, 16 ), 255 );

  bool color_image = f1.size() > 4 && f1.compare( f1.size() - 4, 4, ".ppm" ) == 0;
  bool grey_image  = f1.size() > 4 && f1.compare( f1.size() - 4, 4, ".pgm" ) == 0;
  if ( ! color_image && ! grey_image ) 
    {
      trace.error() << "Input image file must be either a PGM (grey-level) or a PPM (color) image with these extensions."
                    << endl;
      return 2;
    }

  KSpace K;
  ATu2v0< KSpace > AT( verb );
  Domain domain;
  AT.setMetricAverage( metric );
  
  typedef ATu2v0<KSpace>::Calculus Calculus;
  typedef ImageContainerBySTLVector<Domain, Color> ColorImage;
  typedef ImageContainerBySTLVector<Domain, unsigned char> GreyLevelImage;
  //---------------------------------------------------------------------------
  if ( color_image ) 
    {
      trace.beginBlock("Reading PPM image");
      ColorImage image = PPMReader<ColorImage>::importPPM( f1 );
      trace.endBlock();
      trace.beginBlock("Building AT");
      domain = image.domain();
      K.init( domain.lowerBound(), domain.upperBound(), true );
      AT.init( K );
      AT.addInput( image, [] ( Color c ) -> double { return ((double) c.red())   / 255.0; } );
      AT.addInput( image, [] ( Color c ) -> double { return ((double) c.green()) / 255.0; } );
      AT.addInput( image, [] ( Color c ) -> double { return ((double) c.blue())  / 255.0; } );
      trace.endBlock();
    }
  else if ( grey_image ) 
    {
      trace.beginBlock("Reading PGM image");
      GreyLevelImage image = PGMReader<GreyLevelImage>::importPGM( f1 );
      trace.endBlock();
      trace.beginBlock("Building AT");
      domain = image.domain();
      K.init( domain.lowerBound(), domain.upperBound(), true );
      AT.init( K );
      AT.addInput( image, [] (unsigned char c ) { return ((double) c) / 255.0; } );
      trace.endBlock();
    }

  //---------------------------------------------------------------------------
  if ( snr && color_image ) 
    {
      trace.beginBlock("Reading ideal PPM image");
      ColorImage image = PPMReader<ColorImage>::importPPM( isnr );
      trace.endBlock();
      AT.addInput( image, [] ( Color c ) -> double { return ((double) c.red())   / 255.0; }, true );
      AT.addInput( image, [] ( Color c ) -> double { return ((double) c.green()) / 255.0; }, true );
      AT.addInput( image, [] ( Color c ) -> double { return ((double) c.blue())  / 255.0; }, true );
    }
  else if ( snr && grey_image ) 
    {
      trace.beginBlock("Reading ideal PGM image");
      GreyLevelImage image = PGMReader<GreyLevelImage>::importPGM( isnr );
      trace.endBlock();
      AT.addInput( image, [] (unsigned char c ) { return ((double) c) / 255.0; }, true );
    }
  
  //---------------------------------------------------------------------------
  // Prepare zoomed output domain
  Domain out_domain( pix_sz * domain.lowerBound(), 
                     pix_sz * domain.upperBound() + Point::diagonal( pix_sz ) );
  //---------------------------------------------------------------------------
  AT.setUFromInput();
  double g_snr = snr ? AT.computeSNR() : 0.0;

  if ( vm.count( "inpainting-mask" ) )
    {
      string fm  = vm[ "inpainting-mask" ].as<string>();
      trace.beginBlock("Reading inpainting mask");
      GreyLevelImage mask = GenericReader<GreyLevelImage>::import( fm );
      trace.endBlock();
      Calculus::PrimalForm2 m( AT.calculus );
      for ( Calculus::Index index = 0; index < m.myContainer.rows(); index++)
        {
          auto cell = m.getSCell( index );
          double col = ((double) mask( K.sCoords( cell ) )) / 255.0;
          m.myContainer( index ) = col > 0.0 ? 1.0 : 0.0;
        }
      AT.setAlpha( a, m );
      AT.setUFromInputAndMask();
      if ( grey_image )
        {
          ostringstream ossGM;
          ossGM << boost::format("%s-g-mask.pgm") %f2;
          GreyLevelImage image_mg( domain );
          const Calculus::PrimalForm2 mg = functions::dec::diagonal( m ) * AT.getG( 0 );
          functions::dec::form2ToGreyLevelImage
            ( AT.calculus, mg, image_mg, 0.0, 1.0, 1 ); 
          PGMWriter<GreyLevelImage>::exportPGM( ossGM.str(), image_mg );
        }
      else if ( color_image )
        {
          ostringstream ossGM;
          ossGM << boost::format("%s-g-mask.ppm") %f2;
          ColorImage image_mg( domain );
          const Calculus::PrimalForm2 mg0 = functions::dec::diagonal( m ) * AT.getG( 0 );
          const Calculus::PrimalForm2 mg1 = functions::dec::diagonal( m ) * AT.getG( 1 );
          const Calculus::PrimalForm2 mg2 = functions::dec::diagonal( m ) * AT.getG( 2 );
          functions::dec::threeForms2ToRGBColorImage
            ( AT.calculus, mg0, mg1, mg2, image_mg, 0.0, 1.0, 1 ); 
          PPMWriter<ColorImage, functors::Identity >::exportPPM( ossGM.str(), image_mg );
        }
    }
  else 
    AT.setAlpha( a );
  
  trace.info() << AT << std::endl;
  double n_v = 0.0;
  double eps = 0.0;
  while ( l1 >= l2 )
    {
      trace.info() << "************ lambda = " << l1 << " **************" << endl;
      AT.setLambda( l1 );
      for ( eps = e1; eps >= e2; eps /= er )
        {
          trace.info() << "  ======= epsilon = " << eps << " ========" << endl;
          AT.setEpsilon( eps );
          int n = 0;
          do {
            trace.progressBar( n, nbiter );
            AT.solveU();
            AT.solveV();
            AT.checkV();
            n_v = AT.computeVariation();
          } while ( ( n_v > 0.0001 ) && ( ++n < nbiter ) );
          trace.progressBar( n, nbiter );
          trace.info() << "[#### last variation = " << n_v << " " << endl;
        }
      if ( grey_image ) 
        {
          if ( verb > 0 ) trace.beginBlock("Writing u[0] as PGM image");
          ostringstream ossU, ossV, ossW;
          ossU << boost::format("%s-a%.5f-l%.7f-u.pgm") % f2 % a % l1;
          ossV << boost::format("%s-a%.5f-l%.7f-u-v.pgm") % f2 % a % l1;
          ossW << boost::format("%s-a%.5f-l%.7f-u-v.ppm") % f2 % a % l1;
          const Calculus::PrimalForm2 u = AT.getU( 0 );
          const Calculus::PrimalForm1 v = AT.M01 * AT.getV();
          // Restored image
          GreyLevelImage image_u( domain );
          functions::dec::form2ToGreyLevelImage
            ( AT.calculus, u, image_u, 0.0, 1.0, 1 ); 
          PGMWriter<GreyLevelImage>::exportPGM( ossU.str(), image_u );
          // Zoomed restored image with discontinuities (in black).
          GreyLevelImage image_uv( out_domain );
          functions::dec::form2ToGreyLevelImage
            ( AT.calculus, u, image_uv, 0.0, 1.0, pix_sz ); 
          functions::dec::primalForm1ToGreyLevelImage
            ( AT.calculus, v, image_uv, 0.0, 1.0, pix_sz ); 
          PGMWriter<GreyLevelImage>::exportPGM( ossV.str(), image_uv );
          // Zoomed restored image with discontinuities (in specified color).
          ColorImage cimage( out_domain );
          functions::dec::threeForms2ToRGBColorImage
            ( AT.calculus, u, u, u, cimage, 0.0, 1.0, pix_sz ); 
          functions::dec::primalForm1ToRGBColorImage
            ( AT.calculus, v, cimage, color_v, 0.0, 1.0, pix_sz ); 
          PPMWriter<ColorImage, functors::Identity >::exportPPM( ossW.str(), cimage );
          if ( verb > 0 ) trace.endBlock();
        }
      else if ( color_image )
        {
          if ( verb > 0 ) trace.beginBlock("Writing u[0,1,2] as PGM image");
          ostringstream ossU, ossV;
          ossU << boost::format("%s-a%.5f-l%.7f-u.ppm") % f2 % a % l1;
          ossV << boost::format("%s-a%.5f-l%.7f-u-v.ppm") % f2 % a % l1;
          const Calculus::PrimalForm2 u0 = AT.getU( 0 );
          const Calculus::PrimalForm2 u1 = AT.getU( 1 );
          const Calculus::PrimalForm2 u2 = AT.getU( 2 );
          const Calculus::PrimalForm1 v  = AT.M01 * AT.getV();
          // Restored image
          ColorImage image_u( domain );
          functions::dec::threeForms2ToRGBColorImage
            ( AT.calculus, u0, u1, u2, image_u, 0.0, 1.0, 1 ); 
          PPMWriter<ColorImage, functors::Identity >::exportPPM( ossU.str(), image_u );
          ColorImage image_uv( out_domain );
          functions::dec::threeForms2ToRGBColorImage
            ( AT.calculus, u0, u1, u2, image_uv, 0.0, 1.0, pix_sz ); 
          functions::dec::primalForm1ToRGBColorImage
            ( AT.calculus, v, image_uv, color_v, 0.0, 1.0, pix_sz ); 
          PPMWriter<ColorImage, functors::Identity >::exportPPM( ossV.str(), image_uv );
          if ( verb > 0 ) trace.endBlock();
        }
      // Compute SNR if possible
      if ( snr )
        {
          double u_snr = AT.computeSNR();
          trace.info() << "- SNR of u = " << u_snr << "   SNR of g = " << g_snr << endl;
        }
      l1 /= lr;
    }
  return 0;
}
