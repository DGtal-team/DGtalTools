/*
 * Exemple de commande d'execution du code :
 * ./imageProcessing/at-u2-v0-m -i ./../Images/CarreSimple/Carre.pgm --epsilon-1 2 --epsilon-2 1 -v 0 -o ./../IPResultat/ -a 1 -l 0.1 --multiresolution true -t "8 16 32 64"
 */

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
@page DocATu2v0 imageProcessing/at-u2-v0

@brief Computes a piecewise smooth approximation of a grey-level or color image, by optimizing the Ambrosio-Tortorelli functional (with u a 2-form and v a 0-form).

@writers Marion Foare, Jacques-Olivier Lachaud

@b Usage: at-u2-v0 -i [input.pgm]

(for grey-level image restoration)

@b Usage: at-u2-v0 -i [input.ppm]

(for color image restoration)

The Ambrosio-Tortorelli functional is a classical relaxation of the
Mumford-Shah functional.

Given an input grayscale image, defined in an open bounded domain
\f$ \Omega \f$, we represent its gray levels by a function \f$ g
\in L^{\infty}(\Omega) \f$. In the Ambrosio-Tortorelli functional [1],
one wants to find a function \f$ u \in SBV(\Omega) \f$ which is a
smooth approximation of the input image \f$ g \f$.
The Ambrosio-Tortorelli functional [1] is defined by
\f[
  \displaystyle
  AT_{\varepsilon}(u,v)	= \int_\Omega \alpha |u-g|^2 + v^2 |\nabla u|^2
  + \lambda \varepsilon |\nabla v|^2 + \frac{\lambda}{4 \varepsilon} |1-v|^2 dx,
\f]
for functions \f$ u,v \in W^{1,2}(\Omega)\f$ with \f$ 0 \leq v \leq 1 \f$.


In AT functional, function \f$ v \f$ is a smooth approximation
of the set of discontinuities, and takes value close to 0 in this set,
while being close to 1 outside discontinuities. A remarkable property
of this functional is that it \f$ \Gamma \f$-converges to (a
relaxation of) MS functional as \f$ \varepsilon \f$ tends to 0 (see [1]).
The intuition is that a large \f$ \varepsilon \f$ induces a solution
with a fuzzy set of discontinuities, which is then progressively
narrowed to the crisp 1-dimensional set of discontinuites as
\f$ \varepsilon \f$ goes to 0.

We discretize AT with discrete calculus and we set \f$ u \f$ and \f$ g
\f$ to live on the faces and \f$ v \f$ to live on the vertices and
edges. Pixels are faces, so functions \f$ u \f$ and \f$ g \f$ are
2-forms since they represent the gray levels of each pixel. On the
contrary, we set \f$ v \f$ in-between cells of non null measure, so in
this case on vertices as a 0-form, and on edges by averaging with \f$
\mathbf{M} \f$. We call this formulation AT20. The DEC reformulation
is straightforward, except for the second term, where we use matrix
\f$ \mathbf{M} \f$ to transport the 0-form \f$ v \f$ onto edges :

\f[
  \displaystyle
  AT20(u,v) = \Sigma_{i=1}^n
      \alpha \langle u_i - g_i , u_i - g_i \rangle_2
    + \langle \mathbf{M} v , \bar{\mathbf{\star}} \bar{\mathbf{d_0}}
      \mathbf{\star} u_i \rangle_1 ^2 \\
    + \lambda \varepsilon \langle \mathbf{d_0} v , \mathbf{d_0} v \rangle_1
    + \frac{\lambda}{4\varepsilon} \langle 1 - v , 1 - v \rangle_0.
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
  -q [ --lambda-ratio ] arg (=1.414213) the division ratio for lambda from l1
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

@b example:

@b example:

\code
./imageProcessing/at-u2-v0 -i ../imageProcessing/Images/degrade-b04.pgm --image-snr ../imageProcessing/Images/degrade.pgm -a 0.05 --epsilon-1 4 --epsilon-2 0.25 -l 0.0075 -p 2 -c 0xff0000 -o degrade
\endcode

<center>
<table>
<tr>
<td> Input image \a g </td>
<td> Reconstructed image \a u </td>
<td> Perfect image </td>
</tr>
<tr>
<td>@image html degrade-b04.png "Input image (noise = 0.4)"</td>
<td>@image html degrade-a0.05000-l0.0075000-u2.png "AT20 alpha=0.05 lambda=0.0075 "</td>
<td>@image html degrade.png "Perfect image"</td>
</tr>
<tr>
<td> SNR of \a g = 21.9183 </td>
<td> SNR of \a u = 34.3655 </td>
<td> Perfect image </td>
</tr>
</table>
</center>

@note Other restoration examples, parameter analysis, and image
inpainting examples may be found in \ref moduleAT.

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
    double h;
    bool multires = vm.count( "multiresolution" );
    string taille = multires ? vm[ "images-size" ].as<string>() : "16 32 64 128";

    // Gestion des tailles pour la multiresolution
    std::stringstream iss( taille );
    int number;
    std::vector<int> mySizes;

    // Copie des variables qui sont mene a bouger
    string f1_copy = f1;
    string f2_copy = f2;
    double l1_copy = l1;

    // Creating memory vectors for keep u and v for each size
    std::vector< Calculus::PrimalForm2 > u_memory;
    std::vector< string > file_restored;
    std::vector< string > file_contours;
    std::vector< Calculus::PrimalForm1 > v_memory;
    std::vector< string > filename1;
    std::vector< string > filename2;


    if(multires){
        // Recuperation des tailles dans un vecteur de int : mySizes
        while ( iss >> number )
            mySizes.push_back( number );
        size_t lastindex = f1.find_last_of(".");
        for( int i = 0 ; i < mySizes.size() ; i++ ){
            // Creation d'un veceteur de noms pour f1
            f1 = f1_copy;
            f1.insert( lastindex , std::to_string( mySizes[i] ) );
            filename1.push_back( f1 );
            // Creation d'un veceteur de noms pour f2
            f2 = f2_copy;
            f2.append( std::to_string( mySizes[i] ) );
            filename2.push_back( f2 );
        }
    }

    trace.info() << endl;
    trace.beginBlock("Temps total");
    int iterator_size = 0;
    do {
        trace.info() << "-------------------------------------------------------------------------------------" << endl;
        trace.info() << endl;
        trace.beginBlock("Temps pour une image");

        // Initialisation
        f1 = f1_copy;
        f2 = f2_copy;
        l1 = l1_copy;

        // Gestion de la taille de fichier pour la multiresolution
        if(multires){
            size_t lastindex = f1.find_last_of(".");                                // Position de l'extension
            f1.insert( lastindex , std::to_string( mySizes[ iterator_size ] ) );    // ENTRANT : Insertion de la taille avant l'extension
            f2.append( std::to_string( mySizes[ iterator_size ] ) );                // SORTANT : Insertion de la taille a la fin
        }

        if (verb > 0) trace.info()  << "Iterator_size = " << iterator_size << endl
                                    << "Fichier entrant = " << f1 << endl
                                    << "Nom images resultantes = " << f2 << endl
                                    << endl;

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


        //---------------------------------------------------------------------------
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

        double one_h = 1 / h ;
        double one_h2 = 1 / (h*h) ;
        double lh = l1 * one_h ;
        double ah = a + one_h2 ;
        double e1h = e1 * h ;
        double e2h = e2 * h ;



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
        trace.beginBlock("Initialisation de u et v");
        if( multires && (iterator_size > 0) ){
            trace.info() << "Utilisation du resultat precedent. " << endl;
            if( grey_image ){
                trace.info() << "lecture de l'image resataure a l'etape precedente." << endl;
                GreyLevelImage restoredImage = PGMReader<GreyLevelImage>::importPGM(file_restored[iterator_size-1]);
                ColorImage restoredContour = PPMReader<ColorImage>::importPPM(file_contours[iterator_size-1]);
                AT.setUFromImage( restoredImage, [] (unsigned char c ) { return (double) c / 255.0; } );
                AT.setVFromImage(restoredContour ,[] ( Color c ) { return (double) c.red()   / 255.0; });
            }else if ( color_image ){
                // TODO : a faire
                trace.info() << "lecture de l'image resataure a l'etape precedente." << endl;
                ColorImage restoredImage = PPMReader<ColorImage>::importPPM(file_restored[iterator_size-1]);
                ColorImage restoredContour = PPMReader<ColorImage>::importPPM(file_contours[iterator_size-1]);
                AT.setUFromImage( restoredImage, [] ( Color c ) { return (double) c.red()   / 255.0; } );
                AT.setUFromImage( restoredImage, [] ( Color c ) { return (double) c.green() / 255.0; } );
                AT.setUFromImage( restoredImage, [] ( Color c ) { return (double) c.blue()  / 255.0; } );
                AT.setVFromImage(restoredContour ,[] ( Color c ) { return (double) c.red()   / 255.0; });
            }
        }else{
            trace.info() << "Utilisation de l'image donnee en entree : u = g. " << endl;
            AT.setUFromInput();
        }
        trace.endBlock(); trace.info() << endl;

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
            energieFile << endl;

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
                trace.info() << "[#### last variation = " << n_v << " " << endl;

            }

            trace.endBlock(); trace.info() << endl;  // CPU Time execution

            if ( grey_image )
            {
                if ( verb > 0 ) trace.beginBlock("Writing u[0] as PGM image");
                ostringstream ossU, ossV, ossW;
                /*
                          ossU << boost::format("%s-[h_%.5f]-[a_%.5f]-[l_%.7f]-u.pgm") % f2 % h % a % l1;
                          ossV << boost::format("%s-[h_%.5f]-[a_%.5f]-[l_%.7f]-u-v.pgm") % f2 % h % a % l1;
                          ossW << boost::format("%s-[h_%.5f]-[a_%.5f]-[l_%.7f]-v.ppm") % f2 % h % a % l1;
                */
                ossU << boost::format("%sRestaure.pgm") % f2;
                ossV << boost::format("%sContours.ppm") % f2;
                const Calculus::PrimalForm2 u = AT.getU( 0 );
                const Calculus::PrimalForm1 v = AT.M01 * AT.getV();

                u_memory.push_back(u);
                v_memory.push_back(v);
                file_restored.push_back(ossU.str());
                file_contours.push_back(ossV.str());

                energieFile << Nx << "\t"
                            << Ny << "\t"
                            << a << "\t"
                            << l1 << "\t"
                            << h << "\t"
                            << ah << "\t"
                            << lh << "\t"
                            << e2 << "\t"
                            << e2h << "\t"
                            << AT.computeEnergy() << "\t"
                            << AT.computeLambdaPerimeter() << "\t"
                            << AT.computePerimeter() << "\t"
                            << AT.computeFidelity() << "\t"
                            << AT.computeCrossTerm() << "\t"
                            << AT.computeGradV() << "\t"
                            << AT.computeConstraintV() << "\t"
                            << ( AT.computeEv() / lh ) << "\t"
                            << ( AT.computeCrossTerm() / lh ) << "\t"
                            << ( ( AT.computeCrossTerm() + AT.computeEv() ) / lh );

                trace.info() << endl;
                trace.beginBlock("Calcul d'energie");
                trace.info() << "Energie calculee = " << AT.computeEnergy() << endl;
                trace.info() << "Perimetre*lambda_h = " << AT.computeLambdaPerimeter() << endl;
                trace.info() << "Perimetre = " << (AT.computePerimeter()) << endl;
                trace.info() << "Fidelite = " << AT.computeFidelity() << endl;
                trace.info() << "Cross term = " << AT.computeCrossTerm() << endl;
                trace.info() << "Gradient de V = " << AT.computeGradV() << endl;
                trace.info() << "Contraintes sur V = " << AT.computeConstraintV() << endl;
                trace.endBlock(); trace.info() << endl;

                // Restored image
                GreyLevelImage image_u( domain );
                functions::dec::form2ToGreyLevelImage
                        ( AT.calculus, u, image_u, 0.0, 1.0, 1 );
                PGMWriter<GreyLevelImage>::exportPGM( ossU.str(), image_u );

                // Zoomed restored image with discontinuities (in specified color).
                ColorImage cimage( out_domain );
                functions::dec::primalForm1ToRGBColorImage
                        ( AT.calculus, v, cimage, color_v, 0.0, 1.0, pix_sz );
                PPMWriter<ColorImage, functors::Identity >::exportPPM( ossV.str(), cimage );


                if ( verb > 0 ) trace.endBlock();


/*
// NOEMIE : Tests sur les coordonnes de Khalimsky et digitales ////////////////////////////////////////////////////////
                if(multires && (iterator_size < mySizes.size()-1) ) {
                    trace.beginBlock("TEST NOEMIE : CREATION KFORME");

                    // Definition d'un nouvau KSpace
                    KSpace NewK;                               // le KSpace pour l'espace (N x N)
                    ATu2v0<KSpace> NewAT(verb);
                    Domain NewDomain;
                    GreyLevelImage NewImage = PGMReader<GreyLevelImage>::importPGM(filename1[iterator_size + 1]);
                    NewDomain = NewImage.domain();
                    NewK.init(NewDomain.lowerBound(), NewDomain.upperBound(), true);
                    double Newh = h / 2;
                    NewAT.init(NewK,Newh);
                    Domain new_out_domain( pix_sz * NewDomain.lowerBound(),
                                       pix_sz * NewDomain.upperBound() + Point::diagonal( pix_sz ) );
                    trace.info() << "---------------------------------" << endl
                                 << "filename: " << filename1[iterator_size + 1] << endl
                                 << "New KSpace " << NewK << endl
                                 << "New ATu2v0<KSpace>: " << NewAT << endl
                                 << "New Domain: " << NewDomain << endl
                                 << "New out Domain: " << new_out_domain << endl
                                 << "---------------------------------" << endl;



                    // Definition des primal forme utiles
                    //Calculus::PrimalForm2 pf_previous = u_memory[ iterator_size ];        // la primal forme de dimension (N/2 x N/2)
                    Calculus::PrimalForm2 pf_previous = u;        // la primal forme de dimension (N/2 x N/2)
                    Calculus::PrimalForm2 pf_current = Calculus::PrimalForm2(NewAT.calculus);  // la primal forme de dimension (N x N)

                    for (Calculus::Index idx_previous = 0; idx_previous < pf_previous.myContainer.rows(); idx_previous++) {
                        //for ( Calculus::Index idx_previous = 0; idx_previous < 4; idx_previous++) {
                        trace.info()    << "idx_previous: " << idx_previous << endl
                                        << "u(idx) = " << (u.myContainer(idx_previous)) << endl;


                        // TEST AVEC LES COORDONNEES KHALIMSKY //////////////////////////////////////////////

                        // Recuperation des informations de la Kform precedante
                        trace.info() << "Recuperation des informations de la Kform precedante" << endl;
                        auto spixel = pf_previous.getSCell(idx_previous);   // SCell
                        auto pk = K.sKCoords(spixel);                        // Point : coordonnees Khalimsky
                        trace.info() << "> spixel: " << spixel << endl
                                     << "> pk: " << pk << endl;

                        // Creation des 4 points incident de p
                        trace.info() << "Creation des 4 points resultants dans la nouvelle Kform" << endl;
                        Point pk1(pk[0] * 2 - 1, pk[1] * 2 - 1);
                        Point pk2(pk[0] * 2 - 1, pk[1] * 2 + 1);
                        Point pk3(pk[0] * 2 + 1, pk[1] * 2 - 1);
                        Point pk4(pk[0] * 2 + 1, pk[1] * 2 + 1);

                        trace.info() << "> pk1: " << pk1 << endl
                                     << " > pk2: " << pk2 << endl
                                     << " > pk3: " << pk3 << endl
                                     << " > pk4: " << pk4 << endl;

                        // Recuperation des Cell correspondantes
                        trace.info() << "Recuperation des SCell correspondantes" << endl;
                        Cell c1 = NewK.uCell(pk1);
                        Cell c2 = NewK.uCell(pk2);
                        Cell c3 = NewK.uCell(pk3);
                        Cell c4 = NewK.uCell(pk4);
                        trace.info() << "> c1: " << c1 << endl
                                     << " > c2: " << c2 << endl
                                     << " > c3: " << c3 << endl
                                     << " > c4: " << c4 << endl;

                        // Recuperation des index correspondants aux Cell dans la kforme
                        trace.info() << "Recuperation des index correspondants pour la nouvelle KForm" << endl;
                        auto i1 = NewAT.calculus.getCellIndex(c1);
                        auto i2 = NewAT.calculus.getCellIndex(c2);
                        auto i3 = NewAT.calculus.getCellIndex(c3);
                        auto i4 = NewAT.calculus.getCellIndex(c4);
                        trace.info() << "> i1: " << i1 << endl
                                     << " > i2: " << i2 << endl
                                     << " > i3: " << i3 << endl
                                     << " > i4: " << i4 << endl;

                        // Redefinition des index de la kforme dans l'espace de fin
                        trace.info() << "Redefinition des index de la kforme dans l'espace de fin" << endl;
                        pf_current.myContainer(i1) = pf_previous.myContainer(idx_previous);
                        pf_current.myContainer(i2) = pf_previous.myContainer(idx_previous);
                        pf_current.myContainer(i3) = pf_previous.myContainer(idx_previous);
                        pf_current.myContainer(i4) = pf_previous.myContainer(idx_previous);
                        trace.info() << "> pf_previous.myContainer( " << idx_previous << " ): "
                                     << (pf_previous.myContainer(idx_previous)) << endl
                                     << "> pf_current.myContainer( " << i1 << " ): " << (pf_current.myContainer(i1)) << endl
                                     << " > pf_current.myContainer( " << i2 << " ): " << (pf_current.myContainer(i2)) << endl
                                     << " > pf_current.myContainer( " << i3 << " ): " << (pf_current.myContainer(i3)) << endl
                                     << " > pf_current.myContainer( " << i4 << " ): " << (pf_current.myContainer(i4)) << endl;


                        trace.info() << endl
                                     << "---------------------------------" << endl
                                     << idx_previous << "   " << spixel << "   " << pk << endl
                                     << "    > " << pk1 << "  ->  " << i1 << "    ->  " << (pf_current.getSCell(i1)) << endl
                                     << "    > " << pk2 << "  ->  " << i2 << "    ->  " << (pf_current.getSCell(i2)) << endl
                                     << "    > " << pk3 << "  ->  " << i3 << "    ->  " << (pf_current.getSCell(i3)) << endl
                                     << "    > " << pk4 << "  ->  " << i4 << "    ->  " << (pf_current.getSCell(i4)) << endl
                                     << "---------------------------------" << endl;
                        trace.info() << "pf_previous.myContainer.rows() = " << pf_previous.myContainer.rows() << endl;


                    // TEST AVEC LES COORDONNEES DIGITALES ///////////////////////////////////////////////

                    // Recuperation des informations de la Kform precedante
                    trace.info() << "Recuperation des informations de la Kform precedante" << endl;
                    auto spixel = pf_previous.getSCell(idx_previous);
                    auto p = K.sCoords(spixel);
                    auto s = K.sSign(spixel);
                    trace.info() << "> spixel: " << spixel << endl
                                 << "> p: " << p << endl
                                 << "> s: " << s << endl;
                    // Creation des 4 points resultants dans la nouvelle Kform (coordonnees digitales)
                    trace.info() << "Creation des 4 points resultants dans la nouvelle Kform" << endl;
                    Point p0( p[0]*2    , p[1]*2     );
                    Point p1( p[0]*2+1  , p[1]*2     );
                    Point p2( p[0]*2    , p[1]*2+1   );
                    Point p3( p[0]*2+1  , p[1]*2+1   );
                    trace.info() << "> p0: " << p0 << endl
                                 << "> p1: " << p1 << endl
                                 << "> p2: " << p2 << endl
                                 << "> p3: " << p3 << endl;
                    // Recuperation des SCell correspondantes
                    trace.info() << "Recuperation des SCell correspondantes" << endl;
                    SCell c0 = NewK.sSpel(p0, s);
                    SCell c1 = NewK.sSpel(p1, s);
                    SCell c2 = NewK.sSpel(p2, s);
                    SCell c3 = NewK.sSpel(p3, s);
                    trace.info() << "> c0: avec sSpel : " << c0 << " avec uSpel : " << NewK.uSpel(p0) << endl
                                 << "> c1: avec sSpel : " << c1 << " avec uSpel : " << NewK.uSpel(p1) << endl
                                 << "> c2: avec sSpel : " << c2 << " avec uSpel : " << NewK.uSpel(p2) << endl
                                 << "> c3: avec sSpel : " << c3 << " avec uSpel : " << NewK.uSpel(p3) << endl;
                    // Recuperation des Cell correspondantes
                    trace.info() << "Recuperation des SCell correspondantes" << endl;
                    Cell cc0 = NewK.unsigns(c0);
                    Cell cc1 = NewK.unsigns(c1);
                    Cell cc2 = NewK.unsigns(c2);
                    Cell cc3 = NewK.unsigns(c3);
                    trace.info() << "> cc0: " << cc0  << endl
                                 << "> cc1: " << cc1  << endl
                                 << "> cc2: " << cc2  << endl
                                 << "> cc3: " << cc3  << endl;

                    // Recuperation des index correspondants pour la nouvelle KForm
                    trace.info() << "Recuperation des index correspondants pour la nouvelle KForm" << endl;
                    auto i0 = NewAT.calculus.getCellIndex(NewK.unsigns(c0));
                    auto i1 = NewAT.calculus.getCellIndex(NewK.unsigns(c1));
                    auto i2 = NewAT.calculus.getCellIndex(NewK.unsigns(c2));
                    auto i3 = NewAT.calculus.getCellIndex(NewK.unsigns(c3));
                    trace.info() << "> i0: " << i0 << endl
                                 << "> i1: " << i1 << endl
                                 << "> i2: " << i2 << endl
                                 << "> i3: " << i3 << endl;

                    trace.info()    << endl
                                    << "---------------------------------" << endl
                                    << idx_previous << "   " << spixel << "   " << p << endl
                                    << "    > " << p0 << "  ->  " << i0 << "    ->  " << (pf_current.getSCell(i0)) << "    ->  " << (NewK.sCoords(pf_current.getSCell(i0))) << endl
                                    << "    > " << p1 << "  ->  " << i1 << "    ->  " << (pf_current.getSCell(i1)) << "    ->  " << (NewK.sCoords(pf_current.getSCell(i1))) << endl
                                    << "    > " << p2 << "  ->  " << i2 << "    ->  " << (pf_current.getSCell(i2)) << "    ->  " << (NewK.sCoords(pf_current.getSCell(i2))) << endl
                                    << "    > " << p3 << "  ->  " << i3 << "    ->  " << (pf_current.getSCell(i3)) << "    ->  " << (NewK.sCoords(pf_current.getSCell(i3))) << endl
                                    << "---------------------------------" << endl << endl;


                    }
                    trace.info() << "sortie de boucle" << endl;


                    auto tabPrevious = AT.calculus.getIndexedSCells<2,PRIMAL>();
                    for(int itab = 0 ; itab < tabPrevious.size() ; itab++){
                        trace.info() << "tabPrevious("<<itab<<") = " << tabPrevious[itab] << endl;
                    }
                    trace.info() << endl;
                    auto tabCurrent = NewAT.calculus.getIndexedSCells<2,PRIMAL>();
                    for(int itab = 0 ; itab < tabCurrent.size() ; itab++){
                        trace.info() << "tabCurrent("<<itab<<") = " << tabCurrent[itab] << endl;
                    }
                    trace.info() << endl;


                    trace.info() << "NewK.lowerBound() = " << NewK.lowerBound() << endl;
                    trace.info() << NewK.lowerBound()[0] << endl;
                    trace.info() << NewK.lowerBound()[1] << endl;
                    trace.info() << "NewK.upperBound() = " << NewK.upperBound() << endl;
                    trace.info() << NewK.upperBound()[0] << endl;
                    trace.info() << NewK.upperBound()[1] << endl;

                    trace.endBlock();
                }
// NOEMIE : Tests sur les coordonnes de Khalimsky et digitales ////////////////////////////////////////////////////////
 */

// NOEMIE::TEST DE LA METHODE ADDINPUTFROMIMAGE
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * L'option "multiresolution" va prendre en compte un argument : une chaine de caractere qui donne la liste des tailles a traiter pour une image.
 * Il faut au prealable avoir une pyramide d'image avec un nom generique : <NomImage><TailleImage>.<extension>
 *
 * Exemple d'argument pour la multiresolution :
 * "8 16 32 64" avec le fichier Carre.pgm
 *      Traitera les fichiers :
 *          > Carre8.pgm    : resolution normal initialisee a 0
 *          > Carre16.pgm   : multiresolution avec initialisation selon les resultat de Carre8.pgm
 *          > Carre32.pgm   : multiresolution avec initialisation selon les resultat de Carre16.pgm
 *          > Carre64.pgm   : multiresolution avec initialisation selon les resultat de Carre32.pgm
 *
 * Exemple de commande d'execution :
 * ./imageProcessing/at-u2-v0-m -i ./../Images/CarreSimple/Carre.pgm --epsilon-1 2 --epsilon-2 1 -v 0 -o ./../IPResultat/ -a 1 -l 0.1 --multiresolution true -t "8 16 32"
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
                if( multires && (iterator_size > 0) ) {
                    GreyLevelImage restoredImage = PGMReader<GreyLevelImage>::importPGM(file_restored[iterator_size-1]);
                    //ColorImage restoredImage = PPMReader<ColorImage>::importPPM( file_contours[iterator_size-1] );
                    // Initialisation de u a partir du resultat precedent
                    AT.setUFromImage( restoredImage, [] (unsigned char c ) { return (double) c / 255.0; } );
                    //AT.setUFromImage( restoredImage, [] ( Color c ) -> double { return ((double) c.red())   / 255.0; } );
                    //AT.setUFromImage( restoredImage, [] ( Color c ) -> double { return ((double) c.green()) / 255.0; } );
                    //AT.setUFromImage( restoredImage, [] ( Color c ) -> double { return ((double) c.blue())  / 255.0; } );

                    const Calculus::PrimalForm2 u = AT.getU( 0 );
                    trace.info() << "Forme u : " << u << endl;

                    ostringstream ImU;
                    ImU << boost::format("%1%ImU.ppm") % (f2);
                    ColorImage cimageu(domain);
                    functions::dec::threeForms2ToRGBColorImage(AT.calculus, u, u, u, cimageu, 0.0, 255.0, pix_sz);
                    PPMWriter<ColorImage, functors::Identity>::exportPPM(ImU.str(), cimageu);
                }
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                /*
                          // Zoomed restored image with discontinuities (in black).
                          GreyLevelImage image_uv( out_domain );
                          functions::dec::form2ToGreyLevelImage
                            ( AT.calculus, u, image_uv, 0.0, 1.0, pix_sz );
                          functions::dec::primalForm1ToGreyLevelImage
                            ( AT.calculus, v, image_uv, 0.0, 1.0, pix_sz );
                          PGMWriter<GreyLevelImage>::exportPGM( ossV.str(), image_uv );
                */
                /*
                          functions::dec::threeForms2ToRGBColorImage
                            ( AT.calculus, u, u, u, cimage, 0.0, 1.0, pix_sz );
                          functions::dec::primalForm1ToRGBColorImage
                            ( AT.calculus, v, cimage, color_v, 0.0, 1.0, pix_sz );
                */



            }
            else if ( color_image )
            {
                if ( verb > 0 ) trace.beginBlock("Writing u[0,1,2] as PGM image");
                ostringstream ossU, ossV;
                ossU << boost::format("%sRestaure.ppm") % f2;
                ossV << boost::format("%sContours.ppm") % f2;
                /*
                          ossU << boost::format("%s-[h_%.5f]-[a_%.5f]-[l_%.7f]-u.ppm") % f2 % h % a % l1;
                          ossV << boost::format("%s-[h_%.5f]-[a_%.5f]-[l_%.7f]-u-v.ppm") % f2 % h % a % l1;
                */
                const Calculus::PrimalForm2 u0 = AT.getU( 0 );
                const Calculus::PrimalForm2 u1 = AT.getU( 1 );
                const Calculus::PrimalForm2 u2 = AT.getU( 2 );
                const Calculus::PrimalForm1 v  = AT.M01 * AT.getV();

                // TODO : entrer u et v dans le vecteur memoire

                u_memory.push_back(u0);
                u_memory.push_back(u1);
                u_memory.push_back(u2);
                v_memory.push_back(v);
                file_restored.push_back(ossU.str());
                file_contours.push_back(ossV.str());


                energieFile << Nx << "\t"
                            << Ny << "\t"
                            << a << "\t"
                            << l1 << "\t"
                            << h << "\t"
                            << ah << "\t"
                            << lh << "\t"
                            << e2 << "\t"
                            << e2h << "\t"
                            << AT.computeEnergy() << "\t"
                            << AT.computeLambdaPerimeter() << "\t"
                            << AT.computePerimeter() << "\t"
                            << AT.computeFidelity() << "\t"
                            << AT.computeCrossTerm() << "\t"
                            << AT.computeGradV() << "\t"
                            << AT.computeConstraintV() << "\t"
                            << ( AT.computeEv() / lh ) << "\t"
                            << ( AT.computeCrossTerm() / lh ) << "\t"
                            << ( ( AT.computeCrossTerm() + AT.computeEv() ) / lh );

                trace.info() << endl;
                trace.beginBlock("Calcul d'energie");
                trace.info() << "Energie calculee = " << AT.computeEnergy() << endl;
                trace.info() << "Perimetre*lambda_h = " << AT.computeLambdaPerimeter() << endl;
                trace.info() << "Perimetre = " << (AT.computePerimeter()) << endl;
                trace.info() << "Fidelite = " << AT.computeFidelity() << endl;
                trace.info() << "Cross term = " << AT.computeCrossTerm() << endl;
                trace.info() << "Gradient de V = " << AT.computeGradV() << endl;
                trace.info() << "Contraintes sur V = " << AT.computeConstraintV() << endl;
                trace.endBlock(); trace.info() << endl;


                // Restored image
                ColorImage image_u( domain );
                functions::dec::threeForms2ToRGBColorImage
                        ( AT.calculus, u0, u1, u2, image_u, 0.0, 1.0, 1 );
                PPMWriter<ColorImage, functors::Identity >::exportPPM( ossU.str(), image_u );

                ColorImage image_uv( out_domain );
                functions::dec::primalForm1ToRGBColorImage
                        ( AT.calculus, v, image_uv, color_v, 0.0, 1.0, pix_sz );
                /*          functions::dec::threeForms2ToRGBColorImage
                            ( AT.calculus, u0, u1, u2, image_uv, 0.0, 1.0, pix_sz );
                */
                PPMWriter<ColorImage, functors::Identity >::exportPPM( ossV.str(), image_uv );
                if ( verb > 0 ) trace.endBlock();
            }
            // Compute SNR if possible
            if ( snr )
            {
                double u_snr = AT.computeSNR();
                trace.info() << "- SNR of u = " << u_snr << "   SNR of g = " << g_snr << endl;
                energieFile << "\t" << u_snr << "\t" << g_snr;
            }
            l1 /= lr;
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

/*
// DEFINE METHOD ----------------------------------------------------------------------
Calculus::Index getIndexFromPoint_2PF(Calculus::PrimalForm2 pf2, Point p, KSpace K){
    Calculus::Index idx = 0;
    bool continu = true;
    while (continu && (idx < pf2.myContainer.rows()) ) {
        auto SCell = pf2.getSCell( idx );
        Point pi = K.sKCoords( SCell );
        if( pi==p )
            continu = false;
        else
            idx++;
    }
    return idx;
}

// DEFINE METHOD ----------------------------------------------------------------------
 */