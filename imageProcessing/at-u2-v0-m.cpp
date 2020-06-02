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
 * @file at-u2-v0-m.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 * @author Noemie Tasca (\c noemie.tasca@etu.univ-smb.fr )
 * Master informatique, University of Savoie, France
 *
 * @date 2018/09/13
 *
 * A tool file named at-u2-v0-m.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include <sstream>
#include <string>
#include <functional>
#include <boost/format.hpp>
#include <ctime>
#include <fstream>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

#include "ATu2v0-m.h"

/**
@page DocATu2v0-m imageProcessing/at-u2-v0-m

@brief Computes a piecewise smooth approximation of a grey-level or color image, by optimizing the Ambrosio-Tortorelli functional (with u a 2-form and v a 0-form) including euclidean metric.

@writers Noemie Tasca

@b Usage: at-u2-v0-m -i [input.pgm]

(for grey-level image restoration)

@b Usage: at-u2-v0-m -i [input.ppm]

(for color image restoration)

The Ambrosio-Tortorelli functional is a classical relaxation of the
Mumford-Shah functional.

Given an input grayscale image, defined in an open bounded domain
\f$ \Omega \f$, we represent its gray levels by a function \f$ g
\in L^{\infty}(\Omega) \f$. <br>
In the Ambrosio-Tortorelli functional,
one wants to find a function \f$ u \in SBV(\Omega) \f$ which is a
smooth approximation of the input image \f$ g \f$. <br>
The Ambrosio-Tortorelli functional is defined by
\f[
  \displaystyle
  AT_{\varepsilon}(u,v)	= \int_\Omega \alpha |u-g|^2 + v^2 |\nabla u|^2
  + \lambda \varepsilon |\nabla v|^2 + \frac{\lambda}{4 \varepsilon} |1-v|^2 dx,
\f]
for functions \f$ u,v \in W^{1,2}(\Omega)\f$ with \f$ 0 \leq v \leq 1 \f$.


In AT functional, function \f$ v \f$ is a smooth approximation
of the set of discontinuities, and takes value close to 0 in this set,
while being close to 1 outside discontinuities. <br>
A remarkable property of this functional is that it \f$ \Gamma \f$-converges to (a
relaxation of) MS functional as \f$ \varepsilon \f$ tends to 0. <br>
The intuition is that a large \f$ \varepsilon \f$ induces a solution
with a fuzzy set of discontinuities, which is then progressively
narrowed to the crisp 1-dimensional set of discontinuites as
\f$ \varepsilon \f$ goes to 0.

We discretize AT with discrete calculus and we set \f$ u \f$ and \f$ g
\f$ to live on the faces and \f$ v \f$ to live on the vertices and
edges. <br>
We call this formulation AT20-m :

\f[
  \displaystyle
  \begin{array}{rcl}
        AT_\varepsilon^{2,0} (u,v)
        &=&
        \alpha \, \sum_{i=1}^n \left( u_i - g_i \right)^T \, G_2 \, \left( u_i - g_i \right) \\[0.2cm]
        &+&
        \sum_{i=1}^n  u_i^T \, G_2^T \, B'^T \, diag\left( M_{01} \, v \right)^2 \, B' \, G_2 \, u_i \\[0.2cm]
        &+&
        \lambda \, \varepsilon \, v^T \, A^T \, G_1 \, A \, v + \dfrac{\lambda}{4 \, \varepsilon} \, \left( 1 - v \right)^T \, G_0 \, \left( 1 - v \right)
    \end{array}
\f]

For more details, see \ref moduleAT

\b Allowed \b options \b are:

\code
  -h [ --help ]                         display this message
  -i [ --input ] arg                    the input image PPM filename.
  -m [ --inpainting-mask ] arg          the input inpainting mask filename.
  -o [ --output ] arg (=AT)             the output image basename.
  -l [ --lambda ] arg                   the parameter lambda.
  -1 [ --lambda-1 ] arg (=0.3125)       the initial parameter lambda (l1).
  -2 [ --lambda-2 ] arg (=0.0005)       the final parameter lambda (l2).
  -q [ --lambda-ratio ] arg (=1.414213) the division ratio for lambda from l1 to l2.
  -a [ --alpha ] arg (=1)               the parameter alpha.
  -e [ --epsilon ] arg                  the initial and final parameter epsilon of AT functional at the same time.
  --epsilon-1 arg (=2)                  the initial parameter epsilon.
  --epsilon-2 arg (=0.25)               the final parameter epsilon.
  --epsilon-r arg (=2)                  sets the ratio between two consecutive epsilon values of AT functional.
  -n [ --nbiter ] arg (=10)             the maximum number of iterations.
  --image-snr arg                       the input image without deterioration if you wish to compute the SNR.
  -p [ --pixel-size ] arg (=1)          the pixel size for outputing images (useful when one wants to see the discontinuities v on top of u).
  -c [ --color-v ] arg (=0xff0000)      the color chosen for displaying the singularities v (e.g. red is 0xff0000).
  -v [ --verbose ] arg (=0)             the verbose level (0: silent, 1: less silent, etc).
  -s [ --step ]                         The step size for defined metric tensors (default value is 1/N where N is the size of image
  --multiresolution                     The option for a multiresolution option for solving
  -t [ --images-size ]                  A string that contains various sizes for multiresolution option


\endcode

@b example without multireoslution:
\code
./imageProcessing/at-u2-v0-m -i ./../imageProcessing/Images/CarreSimple/Carre256.pgm --epsilon-1 2 --epsilon-2 1 -v 0 -o ./imageProcessing/simple -a 1 -l 0.1
\endcode

@b example with multiresolution:
\code
./imageProcessing/at-u2-v0-m -i ./../imageProcessing/Images/CarreSimple/Carre.pgm --epsilon-1 2 --epsilon-2 1 -v 0 -o ./imageProcessing/multi -a 1 -l 0.1 --multiresolution true -t "32 64 128"
\endcode

<center>
<table>
<tr>
<td> Input image \a g </td>
<td> Reconstructed image \a u </td>
<td> Perfect image </td>
</tr>
<tr>
<td>@image html DegradeBruit128.png "Input image (noise = 0.4)"</td>
<td>@image html DegradeBruit128Restaure.png "alpha=0.5 lambda=0.1 "</td>
<td>@image html Degrade128.png "Perfect image"</td>
</tr>
</table>
</center>

@note The last step of multiresolution for a noisy image can be problematic because of the law followed by \f$ \alpha \f$ parameter. <br>
@b example
<center>
<table>
<tr>
<td> Input image \a g </td>
<td> Reconstructed image \a u </td>
<td> Perfect image </td>
</tr>
<tr>
<td>@image html degrade-b04.png "Input image (noise = 0.4)"</td>
<td>@image html DegradeBruit256Restaure.png "alpha=0.5 lambda=0.1 "</td>
<td>@image html degrade.png "Perfect image"</td>
</tr>
</table>
</center>



*/

using namespace std;
using namespace DGtal;
using namespace Z2i;

typedef ATu2v0<KSpace>::Calculus Calculus;
typedef ImageContainerBySTLVector<Domain, Color> ColorImage;
typedef ImageContainerBySTLVector<Domain, unsigned char> GreyLevelImage;



int main( int argc, char* argv[] )
{

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
            ("alpha,a", po::value<double>()->default_value( 1.0 ), "the parameter alpha (should be multiplied by N^2 where NxN is the size of the finest image)." )
            ("epsilon,e", po::value<double>(), "the initial and final parameter epsilon of AT functional at the same time." )
            ("epsilon-1", po::value<double>()->default_value( 2.0 ), "the initial parameter epsilon." )
            ("epsilon-2", po::value<double>()->default_value( 0.25 ), "the final parameter epsilon." )
            ("epsilon-r", po::value<double>()->default_value( 2.0 ), "sets the ratio between two consecutive epsilon values of AT functional." )
            ("nbiter,n", po::value<int>()->default_value( 10 ), "the maximum number of iterations." )
            ("image-snr", po::value<string>(), "the input image without deterioration if you wish to compute the SNR." )
            ("pixel-size,p", po::value<int>()->default_value( 1 ), "the pixel size for outputing images (useful when one wants to see the discontinuities v on top of u)." )
            ("color-v,c", po::value<string>()->default_value( "0xff0000" ), "the color chosen for displaying the singularities v (e.g. red is 0xff0000)." )
            ("verbose,v", po::value<int>()->default_value( 0 ), "the verbose level (0: silent, 1: less silent, etc)." )
            ("step,s", po::value<double>(), "the step size." )
            ("multiresolution", po::value<string>(), "option for multiresolution." )
            ("images-size,t", po::value<string>(), "sizes to run for multiresolution." )
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
        cerr << "Usage: " << argv[0] << " -i toto.pgm -o res.pgm\n"

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
             << "Example: ./at-u2-v0-m -i ../Images/cerclesTriangle64b02.pgm -o tmp -a 0.05 -e 1 --lambda-1 0.1 --lambda-2 0.00001"
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
    double h;
    bool multires = vm.count( "multiresolution" );
    string taille = multires ? vm[ "images-size" ].as<string>() : "";


    // Copie des variables qui sont mene a bouger
    string f1_copy = f1;
    string f2_copy = f2;
    double l1_copy = l1;

    // Creating memory vectors for keep u and v for each size
    std::vector< string > file_restored;
    std::vector< string > file_contours;
    std::vector< string > filename1;
    std::vector< string > filename2;


    // Gestion des tailles pour la multiresolution
    std::stringstream iss( taille );
    std::vector<int> mySizes;
    if(multires)
    {
        // Recuperation des tailles dans un vecteur de int : mySizes
        int number;
        while ( iss >> number )
            mySizes.push_back( number );
        size_t lastindex = f1.find_last_of(".");
        for( int i = 0 ; i < mySizes.size() ; i++ )
        {
            // Creation d'un vecteur de noms pour f1
            f1 = f1_copy;
            f1.insert( lastindex , std::to_string( mySizes[i] ) );  // insertion de la taille
            filename1.push_back( f1 );
            // Creation d'un vecteur de noms pour f2
            f2 = f2_copy;
            f2.append( std::to_string( mySizes[i] ) );              // insertion de la taille
            filename2.push_back( f2 );
        }
    }

    trace.info() << endl;
    trace.beginBlock("Temps total");
    int iterator_size = 0;
    do
    {
        trace.info() << "--------------------------------------------------------------------" << endl;
        trace.info() << endl;
        trace.beginBlock("Temps pour une image");
        if(multires) trace.info() << "Taille : " << mySizes[iterator_size] << endl;

        // Initialisation
        f1 = f1_copy;
        f2 = f2_copy;
        l1 = l1_copy;

        // Gestion de la taille de fichier pour la multiresolution
        if(multires)
        {
            size_t lastindex = f1.find_last_of(".");                                // Position de l'extension
            f1.insert( lastindex , std::to_string( mySizes[ iterator_size ] ) );    // ENTRANT : Insertion de la taille avant l'extension
            f2.append( std::to_string( mySizes[ iterator_size ] ) );                // SORTANT : Insertion de la taille a la fin
        }

        trace.info()  << "Fichier entrant (f1) = " << f1 << endl;
        trace.info()  << "Fichiers sortants (f2) = " << f2 << endl << endl;

        // Determinaison du type d'image
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
        int Nx, Ny;


        //------------------------------------------------------------------------------------------------------------
        // Initialisation de g
        trace.info() << endl;
        trace.beginBlock("Initialisation de g");
        trace.info() << "Initialisation de g (addInput) : " << f1 << endl;
        if ( color_image )
        {
            if(verb > 0) trace.beginBlock("Reading PPM image");
            ColorImage image = PPMReader<ColorImage>::importPPM( f1 );
            if(verb > 0)trace.endBlock();
            if(verb > 0)trace.beginBlock("Building AT");
            domain = image.domain();
            Nx = image.domain().upperBound()[0] + 1 - image.domain().lowerBound()[0];
            Ny = image.domain().upperBound()[1] + 1 - image.domain().lowerBound()[1];
            if(verb > 0)trace.info() << "[Dimension de l'image] Nx = " << Nx << endl;
            if(verb > 0)trace.info() << "[Dimension de l'image] Ny = " << Ny << endl;

            // Definition of step
            if ( vm.count( "step" ) )
                h = vm[ "step" ].as<double>();
            else
                h = 1.0 / ( std::max(Nx, Ny) );

            K.init( domain.lowerBound(), domain.upperBound(), true );
            AT.init( K,h );
            AT.addInput( image, [] ( Color c ) -> double { return ((double) c.red())   / 255.0; } );
            AT.addInput( image, [] ( Color c ) -> double { return ((double) c.green()) / 255.0; } );
            AT.addInput( image, [] ( Color c ) -> double { return ((double) c.blue())  / 255.0; } );
            if(verb > 0)trace.endBlock();
        }
        else if ( grey_image )
        {
            if(verb > 0) trace.beginBlock("Reading PGM image");
            GreyLevelImage image = PGMReader<GreyLevelImage>::importPGM( f1 );
            if(verb > 0) trace.endBlock();
            if(verb > 0) trace.beginBlock("Building AT");
            domain = image.domain();

            Nx = image.domain().upperBound()[0] + 1 - image.domain().lowerBound()[0];
            Ny = image.domain().upperBound()[1] + 1 - image.domain().lowerBound()[1];
            if(verb > 0) trace.info() << "[Dimension de l'image] Nx = " << Nx << endl;
            if(verb > 0) trace.info() << "[Dimension de l'image] Ny = " << Ny << endl;

            // Definition of step
            if ( vm.count( "step" ) )
                h = vm[ "step" ].as<double>();
            else
                h = 1.0 / ( std::max(Nx, Ny) );


            K.init( domain.lowerBound(), domain.upperBound(), true );
            AT.init( K,h );
            AT.addInput( image, [] (unsigned char c ) { return ((double) c) / 255.0; } );
            if(verb > 0) trace.endBlock();
        }
        trace.endBlock();
        trace.info() << endl;
        // -----------------------------------------------------------------------------------------------


        //---------------------------------------------------------------------------
        if ( snr && color_image )
        {
            trace.info() << endl;
            trace.beginBlock("Reading ideal PPM image");
            ColorImage image = PPMReader<ColorImage>::importPPM( isnr );
            trace.endBlock(); trace.info() << endl;
            AT.addInput( image, [] ( Color c ) -> double { return ((double) c.red())   / 255.0; }, true );
            AT.addInput( image, [] ( Color c ) -> double { return ((double) c.green()) / 255.0; }, true );
            AT.addInput( image, [] ( Color c ) -> double { return ((double) c.blue())  / 255.0; }, true );
        }
        else if ( snr && grey_image )
        {
            trace.info() << endl;
            trace.beginBlock("Reading ideal PGM image");
            GreyLevelImage image = PGMReader<GreyLevelImage>::importPGM( isnr );
            trace.endBlock(); trace.info() << endl;
            AT.addInput( image, [] (unsigned char c ) { return ((double) c) / 255.0; }, true );
        }

        //---------------------------------------------------------------------------
        // Prepare zoomed output domain
        Domain out_domain( pix_sz * domain.lowerBound(),
                           pix_sz * domain.upperBound() + Point::diagonal( pix_sz ) );

        //---------------------------------------------------------------------------
//      // Initialisation de U
        trace.info() << endl;
        if(verb > 0) trace.beginBlock("Initialisation de u et v");
        if( multires && (iterator_size > 0) )
        {
            if(verb > 0) trace.info() << "Utilisation du resultat precedent. " << endl;
            if( grey_image )
            {
                if(verb > 0) trace.info() << "Resultat a l'etape precedente : " << endl;
                if(verb > 0) trace.info() << "Restauree : " << file_restored[iterator_size-1] << endl;
                if(verb > 0) trace.info() << "Contours : " << file_contours[iterator_size-1] << endl;
                GreyLevelImage restoredImage = PGMReader<GreyLevelImage>::importPGM(file_restored[iterator_size-1]);
                ColorImage restoredContour = PPMReader<ColorImage>::importPPM(file_contours[iterator_size-1]);
                AT.setUFromImage( restoredImage, [] (unsigned char c ) { return (double) c / 255.0; } );
                AT.setVFromImage(restoredContour ,[] ( Color c ) { return (double) c.red()   / 255.0; });
            }
            else if ( color_image )
            {
                if(verb > 0) trace.info() << "Resultat a l'etape precedente : " << endl;
                if(verb > 0) trace.info() << "Restauree : " << file_restored[iterator_size-1] << endl;
                if(verb > 0) trace.info() << "Contours : " << file_contours[iterator_size-1] << endl;
                ColorImage restoredImage = PPMReader<ColorImage>::importPPM(file_restored[iterator_size-1]);
                ColorImage restoredContour = PPMReader<ColorImage>::importPPM(file_contours[iterator_size-1]);
                AT.setUFromImage( restoredImage, [] ( Color c ) { return (double) c.red()   / 255.0; } );
                AT.setUFromImage( restoredImage, [] ( Color c ) { return (double) c.green() / 255.0; } );
                AT.setUFromImage( restoredImage, [] ( Color c ) { return (double) c.blue()  / 255.0; } );
                AT.setVFromImage( restoredContour ,[] ( Color c ) { return (double) c.red()   / 255.0; });
            }
        }
        else
        {
            if(verb > 0) trace.info() << "Utilisation de l'image donnee en entree : u = g. " << endl;
            AT.setUFromInput();
        }
        if (verb>0 )trace.endBlock();
        trace.info() << endl;




        // -------------------------------------------------------
        double g_snr = snr ? AT.computeSNR() : 0.0;

        if ( vm.count( "inpainting-mask" ) )
        {
            string fm  = vm[ "inpainting-mask" ].as<string>();
            trace.info() << endl;
            trace.beginBlock("Reading inpainting mask");
            GreyLevelImage mask = GenericReader<GreyLevelImage>::import( fm );
            trace.endBlock(); trace.info() << endl;
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

        ofstream energieFile;
        energieFile.open("./imageProcessing/Resultat/calcul_energie.txt", ofstream::out | ofstream::app);


        double n_v = 0.0;
        double eps = 0.0;
        while ( l1 >= l2 )
        {
            trace.info() << endl;
            trace.beginBlock("Minimisation : resolution alternee");   // CPU Time execution

            AT.setLambda( l1 );
            for ( eps = e1; eps >= e2; eps /= er ){
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
                trace.info() << ">> last variation = " << n_v << " " << endl;
                trace.info() << ">> number iteration = " << n << " (nbitermax="<<nbiter<<") " << endl;
                trace.info() << ">> energie =         " << AT.computeEnergy() <<  "       " << (AT.computeEnergy()*h) << endl;
            }
            trace.endBlock(); trace.info() << endl;  // CPU Time execution


            if ( grey_image )
            {
                if ( verb > 0 ) trace.beginBlock("Writing u[0] as PGM image");
                ostringstream ossU, ossV, ossW;

                ossU << boost::format("%sRestaure.pgm") % f2;
                ossV << boost::format("%sContours.ppm") % f2;
                const Calculus::PrimalForm2 u = AT.getU( 0 );
                const Calculus::PrimalForm1 v = AT.M01 * AT.getV();

                file_restored.push_back(ossU.str());
                file_contours.push_back(ossV.str());

                // Restored image
                GreyLevelImage image_u( domain );
                functions::dec::form2ToGreyLevelImage( AT.calculus, u, image_u, 0.0, 1.0, 1 );
                PGMWriter<GreyLevelImage>::exportPGM( ossU.str(), image_u );

                // Image with discontinuities (in specified color).
                ColorImage cimage( out_domain );
                functions::dec::primalForm1ToRGBColorImage( AT.calculus, v, cimage, color_v, 0.0, 1.0, pix_sz );
                PPMWriter<ColorImage, functors::Identity >::exportPPM( ossV.str(), cimage );

                if ( verb > 0 ) trace.endBlock();
            }
            else if ( color_image )
            {
                if ( verb > 0 ) trace.beginBlock("Writing u[0,1,2] as PGM image");
                ostringstream ossU, ossV;
                ossU << boost::format("%sRestaure.ppm") % f2;
                ossV << boost::format("%sContours.ppm") % f2;

                const Calculus::PrimalForm2 u0 = AT.getU( 0 );
                const Calculus::PrimalForm2 u1 = AT.getU( 1 );
                const Calculus::PrimalForm2 u2 = AT.getU( 2 );
                const Calculus::PrimalForm1 v  = AT.M01 * AT.getV();

                file_restored.push_back(ossU.str());
                file_contours.push_back(ossV.str());


                // Restored image
                ColorImage image_u( domain );
                functions::dec::threeForms2ToRGBColorImage( AT.calculus, u0, u1, u2, image_u, 0.0, 1.0, 1 );
                PPMWriter<ColorImage, functors::Identity >::exportPPM( ossU.str(), image_u );

                // Image with discontinuities (in specified color).
                ColorImage image_v( out_domain );
                functions::dec::primalForm1ToRGBColorImage( AT.calculus, v, image_v, color_v, 0.0, 1.0, pix_sz );
                PPMWriter<ColorImage, functors::Identity >::exportPPM( ossV.str(), image_v );

                if ( verb > 0 ) trace.endBlock();
            }
            
            energieFile << Nx << "\t"
                        << h << "\t"
                        << AT.computeEnergy() << "\t"
                        << AT.computePerimeter() << "\t"
                        << AT.computeVariance() << "\t"
                        << AT.computeFidelity() << "\t"
                        << AT.computeCrossTerm() << "\t"
                        << AT.computeGradV() << "\t"
                        << AT.computeConstraintV() << "\t";
            
            trace.info() << endl;
            trace.beginBlock("Calcul d'energie");
            trace.info() << "Energie calculee = " << AT.computeEnergy()  << endl;
            trace.info() << "Perimetre =        " << AT.computePerimeter() << endl;
            trace.info() << "Variance =         " << AT.computeVariance() << endl;
            trace.info() << "Fidelite =         " << AT.computeFidelity() << endl;
            trace.info() << "Cross term =       " << AT.computeCrossTerm() << endl;
            trace.info() << "Gradient de V =    " << AT.computeGradV() << endl;
            trace.info() << "Contraintes V =    " << AT.computeConstraintV()  << endl;
            trace.endBlock();
            trace.info() << endl;
            

            // Compute SNR if possible
            if ( snr )
            {
                double u_snr = AT.computeSNR();
                trace.info() << "- SNR of u = " << u_snr << "   SNR of g = " << g_snr << endl;
                energieFile << "\t" << u_snr << "\t" << g_snr;
            }
            l1 /= lr;
            energieFile << endl;
        }


        energieFile.close();

        iterator_size ++;

        trace.endBlock();
        trace.info() << endl;
    }
    while (iterator_size < mySizes.size());
    trace.endBlock();
    trace.info() << endl;

    return 0;
}
