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
 * @file vol2normalField.cpp
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/04/08
 *
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iterator>
#include "DGtal/base/Common.h"

#include "DGtal/topology/CanonicDigitalSurfaceEmbedder.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/DigitalSetBoundary.h"
#include "DGtal/topology/ImplicitDigitalSurface.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/ExplicitDigitalSurface.h"
#include "DGtal/topology/LightExplicitDigitalSurface.h"
#include "DGtal/graph/BreadthFirstVisitor.h"
#include "DGtal/topology/helpers/FrontierPredicate.h"
#include "DGtal/topology/helpers/BoundaryPredicate.h"
#include "DGtal/graph/CUndirectedSimpleLocalGraph.h"
#include "DGtal/graph/CUndirectedSimpleGraph.h"

#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/SimpleThresholdForegroundPredicate.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/kernel/CanonicEmbedder.h"

#include "DGtal/geometry/surfaces/estimation/DigitalSurfaceEmbedderWithNormalVectorEstimator.h"

#include "DGtal/geometry/surfaces/estimation/LocalEstimatorFromSurfelFunctorAdapter.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/ElementaryConvolutionNormalVectorEstimator.h"
#include "DGtal/geometry/volumes/distance/LpMetric.h"
#include "CLI11.hpp"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z3i;

/**
 @page vol2normalField vol2normalField
 
 @brief Generates normal vector field from a vol file using DGtal library.

 It will output the embedded vector field (Gaussian convolution on elementary normal vectors)
 an OFF file, and a TXT normal vector file (theta, phi in degree).


 @b Usage: 	vol2normalField[options] --input <volFileName> --o <outputFileName> 


 @b Allowed @b options @b are : 
 @code
  Positionals:
  1 TEXT:FILE REQUIRED                  Input vol file.
  2 TEXT REQUIRED                       Output file.

  Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         Input vol file.
  -o,--output TEXT REQUIRED             Output file.
  -l,--level UINT=0                     Iso-level for the surface construction (default 0).
  -s,--sigma FLOAT=5                    Sigma parameter of the Gaussian kernel (default 5.0).
  --exportOriginAndExtremity            exports the origin and extremity of the vector fields when exporting the vector field in TXT format (useful to be displayed in other viewer like meshViewer).
  -N,--vectorsNorm FLOAT=1              set the norm of the exported vectors in TXT format (when the extremity points are exported with --exportOriginAndExtremity). By using a negative value you will invert the direction of the vectors (default 1.0).
  -n,--neighborhood UINT=10             Size of the neighborhood for the convolution (distance on surfel graph, default 10). 
                                  (distance on surfel graph).
 @endcode

 @b Example: 

 We consider the generation of normal vector field from the Iso-level 40 and export the vectors with a norm = -3 (negative value to invert the normal direction).
 
 @code
 $ vol2normalField -i $DGtal/examples/samples/lobster.vol -o lobTreshold40 -l 40 --exportOriginAndExtremity  -N -3
 @endcode


 You can use the too meshViewer to display the resulting vector field with the Iso-level surface:
@code
$ meshViewer lobTreshold40.off -f lobTreshold40.txt  --vectorFieldIndex 2 3 4 5 6 7  -n
@endcode

 You should obtain such a result:
 @image html resVol2normalField.png "Resulting vector field visualization."
 
 @see
 @ref vol2normalField.cpp

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


int main ( int argc, char**argv )
{

    // parse command line CLI ----------------------------------------------
    CLI::App app;
    std::string filename;
    std::string outputFileName;
    unsigned int level {0};
    double sigma {5.0};
    unsigned int neighborhood {10};
    double normExport {1.0};

    app.description("Generates normal vector field from a vol file using DGtal library.\n Typical use example:\n \t vol2normalField[options] --input <volFileName> --o <outputFileName>\n");
    app.add_option("-i,--input,1",filename,"Input vol file.")->required()->check(CLI::ExistingFile);
    app.add_option("-o,--output,2",outputFileName,"Output file.")->required();
    app.add_option("--level,-l",level,"Iso-level for the surface construction (default 0).",true);
    app.add_option("--sigma,-s", sigma,"Sigma parameter of the Gaussian kernel (default 5.0).",true);
    auto expOpt = app.add_flag("--exportOriginAndExtremity", "exports the origin and extremity of the vector fields when exporting the vector field in TXT format (useful to be displayed in other viewer like meshViewer).");
    app.add_option("--vectorsNorm,-N", normExport, "set the norm of the exported vectors in TXT format (when the extremity points are exported with --exportOriginAndExtremity). By using a negative value you will invert the direction of the vectors (default 1.0).",true); 
    app.add_option("--neighborhood,-n", neighborhood,"Size of the neighborhood for the convolution (distance on surfel graph, default 10).",true);  

    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------  

    typedef ImageSelector < Z3i::Domain, unsigned char>::Type Image;
    Image image = VolReader<Image>::importVol ( filename );

    trace.info() <<image<<std::endl;

    functors::SimpleThresholdForegroundPredicate<Image> simplePredicate ( image, level );

    KSpace ks;
    bool space_ok = ks.init ( image.domain().lowerBound(),
                              image.domain().upperBound(), true );
    if ( !space_ok )
    {
        trace.error() << "Error in the Khamisky space construction."<<std::endl;
        return 2;
    }

    typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
    MySurfelAdjacency surfAdj ( true ); // interior in all directions.

    //Set up digital surface.
    typedef LightImplicitDigitalSurface<KSpace, functors::SimpleThresholdForegroundPredicate<Image>  > MyDigitalSurfaceContainer;
    typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
    SCell bel = Surfaces<KSpace>::findABel ( ks, simplePredicate );

    MyDigitalSurfaceContainer* ptrSurfContainer =
        new MyDigitalSurfaceContainer ( ks, simplePredicate, surfAdj, bel );
    MyDigitalSurface digSurf ( ptrSurfContainer ); // acquired
    MyDigitalSurface::ConstIterator it = digSurf.begin();

    // Embedder definition
    typedef CanonicDigitalSurfaceEmbedder<MyDigitalSurface> SurfaceEmbedder;
    SurfaceEmbedder surfaceEmbedder ( digSurf );

    //Convolution kernel
    typedef typename MyDigitalSurface::Surfel Surfel;
    typedef DGtal::functors::GaussianKernel GaussianFunctor;
    typedef DGtal::functors::ElementaryConvolutionNormalVectorEstimator<Surfel, CanonicSCellEmbedder<KSpace>> Functor;
    typedef LocalEstimatorFromSurfelFunctorAdapter<
        MyDigitalSurfaceContainer, 
        LpMetric<Z3i::Space>,
        Functor, GaussianFunctor> MyGaussianEstimator;

    //Estimator definition

    GaussianFunctor Gkernel(sigma);
    LpMetric<Z3i::Space> l1(1.0);
    CanonicSCellEmbedder<KSpace> embedder(digSurf.container().space());
    Functor estimator(embedder,  1.0);

    MyGaussianEstimator myNormalEstimatorG;
    myNormalEstimatorG.attach(digSurf);
    myNormalEstimatorG.setParams(l1, estimator, Gkernel, neighborhood);

    // Embedder definition
    typedef DigitalSurfaceEmbedderWithNormalVectorEstimator<SurfaceEmbedder,MyGaussianEstimator> SurfaceEmbedderWithGaussianNormal;
    SurfaceEmbedderWithGaussianNormal mySurfelEmbedderG ( surfaceEmbedder, myNormalEstimatorG );

    // Compute normal vector field and displays it.
    myNormalEstimatorG.init ( 1.0, digSurf.begin(), digSurf.end() );

    trace.info() << "Generating the NOFF surface "<< std::endl;
    ofstream out2 ( ( outputFileName + ".off" ).c_str() );
    if ( out2.good() )
      digSurf.exportAs3DNOFF ( out2 ,mySurfelEmbedderG );
    out2.close();

    trace.info() << "Generating the polar coordinates file"<< std::endl;
    ofstream out3 ( ( outputFileName + ".txt" ).c_str() );
    if ( out3.good() )
    {
        MyGaussianEstimator::Quantity res;
        for ( MyDigitalSurface::ConstIterator it =digSurf.begin(),
                itend = digSurf.end(); it != itend; ++it )
        {
            res = myNormalEstimatorG.eval ( it );
            //We output Theta - Phi
            out3<< acos ( res [2] ) *180.0/M_PI <<"  " << ( atan2 ( res [1], res [0] ) + M_PI ) *180.0/M_PI;

            if (expOpt->count()>0)
              {
                res *= normExport;
                out3 << " " << mySurfelEmbedderG(*it)[0]
                     << " " << mySurfelEmbedderG(*it)[1]
                     << " " << mySurfelEmbedderG(*it)[2] << " " 
                     <<  mySurfelEmbedderG(*it)[0]+res[0] << " "
                     << mySurfelEmbedderG(*it)[1]+res[1] <<  " "
                     << mySurfelEmbedderG(*it)[2]+res[2];
              
              }

              out3 <<std::endl;
        }
    }
    out3.close();

    return 0;
}


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

