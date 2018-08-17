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

#include "DGtal/geometry/surfaces/estimation/CNormalVectorEstimator.h"
#include "DGtal/geometry/surfaces/estimation/BasicConvolutionWeights.h"
#include "DGtal/geometry/surfaces/estimation/LocalConvolutionNormalVectorEstimator.h"
#include "DGtal/geometry/surfaces/estimation/DigitalSurfaceEmbedderWithNormalVectorEstimator.h"



#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;


/**
 @page vol2normalField vol2normalField
 
 @brief Generates normal vector field from a vol file using DGtal library.

 It will output the embedded vector field (Gaussian convolution on elementary normal vectors)
 an OFF file, and a TXT normal vector file (theta, phi in degree).


 @b Usage: 	vol2normalField[options] --input <volFileName> --o <outputFileName> 


 @b Allowed @b options @b are : 
 @code
  -h [ --help ]                   display this message.
  -i [ --input ] arg              Input vol file.
  -o [ --output ] arg             Output filename.
  -l [ --level ] arg (=0)         Iso-level for the surface construction.
  -s [ --sigma ] arg (=5)         Sigma parameter of the Gaussian kernel.
  --exportOriginAndExtremity      exports the origin and extremity of the 
                                  vector fields when exporting the vector field
                                  in TXT format (useful to be displayed in 
                                  other viewer like meshViewer).
  -N [ --vectorsNorm ] arg (=1)   set the norm of the exported vectors in TXT 
                                  format (when the extremity points are 
                                  exported with --exportOriginAndExtremity). By
                                  using a negative value you will invert the 
                                  direction of the vectors.
  -n [ --neighborhood ] arg (=10) Size of the neighborhood for the convolution 
                                  (distance on surfel graph).
 @endcode

 @b Example: 

 We consider the generation of normal vector field from the Iso-level 40 and export the vectors with a norm = -3 (negative value to invert the normal direction).
 
 @code
 $ vol2normalField -i $DGtal/examples/samples/lobster.vol -o lobTreshold40 -l 40 --exportOriginAndExtremity  -N -3
 @endcode


 You can use the too meshViewer to display the resulting vector field with the Iso-level surface:
@code
$ meshViewer -i lobTreshold40.off -f lobTreshold40.txt  --vectorFieldIndex 2 3 4 5 6 7  -n
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

    // parse command line ----------------------------------------------
    po::options_description general_opt ( "Allowed options are: " );
    general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>(), "Input vol file." )
    ( "output,o", po::value<string>(),"Output filename." )
    ( "level,l", po::value<unsigned int>()->default_value ( 0 ),"Iso-level for the surface construction." )
    ( "sigma,s", po::value<double>()->default_value ( 5.0 ),"Sigma parameter of the Gaussian kernel." )
    ("exportOriginAndExtremity", "exports the origin and extremity of the vector fields when exporting the vector field in TXT format (useful to be displayed in other viewer like meshViewer).") 
      ("vectorsNorm,N", po::value<double>()->default_value(1.0), "set the norm of the exported vectors in TXT format (when the extremity points are exported with --exportOriginAndExtremity). By using a negative value you will invert the direction of the vectors.") 
      ( "neighborhood,n", po::value<unsigned int>()->default_value ( 10 ),"Size of the neighborhood for the convolution (distance on surfel graph)." );

    bool parseOK=true;
    po::variables_map vm;
    try{
      po::store(po::parse_command_line(argc, argv, general_opt), vm);  
    }catch(const std::exception& ex){
      parseOK=false;
      trace.info()<< "Error checking program options: "<< ex.what()<< endl;
    }
    
    po::notify ( vm );
    if (!parseOK ||  vm.count ( "help" ) ||argc<=1 )
    {
        trace.info() << "Generate normal vector field from a vol file using DGtal library."<<std::endl
                     << "It will output the embedded vector field (Gaussian convolution on elementary normal vectors)"<<std::endl
                     << "an OFF file, and a TXT normal vector file (theta, phi in degree)."
                     << std::endl << "Basic usage: "<<std::endl
                     << "\tvol2normalField[options] --input <volFileName> --o <outputFileName> "<<std::endl
                     << general_opt << "\n";
        return 0;
    }

    //Parse options
    if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
    std::string filename = vm["input"].as<std::string>();
    if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
    std::string outputFileName = vm["output"].as<std::string>();

    unsigned int level = vm["level"].as<unsigned int>();
    double sigma = vm["sigma"].as<double>();
    unsigned int neighborhood = vm["neighborhood"].as<unsigned int>();
    double normExport = vm["vectorsNorm"].as<double>();
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
    deprecated::GaussianConvolutionWeights < MyDigitalSurface::Size > Gkernel ( sigma );

    //Estimator definition
    typedef deprecated::LocalConvolutionNormalVectorEstimator  < MyDigitalSurface,
                                                     deprecated::GaussianConvolutionWeights< MyDigitalSurface::Size>  > MyGaussianEstimator;
    BOOST_CONCEPT_ASSERT ( ( concepts::CNormalVectorEstimator< MyGaussianEstimator > ) );
    MyGaussianEstimator myNormalEstimatorG ( digSurf, Gkernel );

    // Embedder definition
    typedef DigitalSurfaceEmbedderWithNormalVectorEstimator<SurfaceEmbedder,MyGaussianEstimator> SurfaceEmbedderWithGaussianNormal;
    SurfaceEmbedderWithGaussianNormal mySurfelEmbedderG ( surfaceEmbedder, myNormalEstimatorG );

    // Compute normal vector field and displays it.
    myNormalEstimatorG.init ( 1.0, neighborhood );

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

            if (vm.count("exportOriginAndExtremity"))
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

