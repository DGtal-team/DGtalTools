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
 * @file volToOff.cpp
 * @ingroup Examples
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/02/06
 *
 * A simple marching cube algorithm for vol files based on digital surfaces.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/IntervalForegroundPredicate.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;



int main( int argc, char** argv )
{

   // parse command line ----------------------------------------------
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>(), "Input vol file." )
    ( "output,o", po::value<string>(),"Output OFF filename." )
    ( "minT,m", po::value<int>()->default_value(0),"minimum value of voxels to define the object (def: 0)." )
    ( "maxT,M", po::value<int>()->default_value(255),"maximum value of voxels to define the object (def: 255)." )
    ( "adjacency,a", po::value<int>()->default_value(0),"0: interior adjacency, 1: exterior adjacency.");

  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify ( vm );
  if (!parseOK || vm.count ( "help" ) ||argc<=1 )
    {
      trace.info() << "Extract dual surface of a vol file for values in ]minT,maxT]."<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\tvol2off --input <volFileName> --o <OffOutputFileName> "<<std::endl
                   << general_opt << "\n";
      return 0;
    }

  if(! vm.count("input")||! vm.count("output"))
    {
      trace.error() << " Input and output filename are required." << endl;
      return 0;
    }

  std::string inputFilename =  vm["input"].as<std::string>();
  std::string outputFilename =  vm["output"].as<std::string>();
  unsigned int minThreshold =  vm["minT"].as<int>();
  unsigned int maxThreshold =  vm["maxT"].as<int>();
  bool intAdjacency = (vm["adjacency"].as<int>() == 0);

  trace.beginBlock( "Reading vol file into an image." );
  typedef ImageSelector < Domain, int>::Type Image;
  Image image = VolReader<Image>::importVol(inputFilename);
  trace.endBlock();

  // Construct the Khalimsky space from the image domain
  KSpace K;
  bool space_ok = K.init( image.domain().lowerBound(),
                          image.domain().upperBound(), true );
  if (!space_ok)
    {
      trace.error() << "Error in the Khamisky space construction."<<std::endl;
      return 2;
    }

  typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
  MySurfelAdjacency surfAdj( intAdjacency ); // interior in all directions.

  trace.beginBlock( "Extracting boundary by scanning the space. " );
  typedef KSpace::Surfel Surfel;
  typedef KSpace::SurfelSet SurfelSet;

  //Implicit surface
  IntervalForegroundPredicate<Image> simplePredicate ( image, minThreshold, maxThreshold );
  SurfelAdjacency< KSpace::dimension > SAdj ( true );
  Surfel bel;
  try{
     bel = Surfaces<KSpace>::findABel ( K, simplePredicate, 10000 );
  }catch(const DGtal::InputException &e){
    trace.error()<< "Error findinf starting bel: "<< e.what()<< endl;
    trace.error() << "The interval you use: ]" << minThreshold << ","<< maxThreshold << "] "<< " produces probably an empty set of voxel. " << endl; 
    return 0;
  }
  LightImplicitDigitalSurface<KSpace, IntervalForegroundPredicate<Image> > LightImplDigSurf ( K, simplePredicate, SAdj, bel );
  DigitalSurface< LightImplicitDigitalSurface<KSpace, IntervalForegroundPredicate<Image> > > surf ( LightImplDigSurf );

  trace.info() << "Digital surface has " << surf.size() << " surfels."
               << std::endl;
  trace.endBlock();

  trace.beginBlock( "Extracting OFF surface" );
  std::ofstream out( outputFilename.c_str() );
  if ( out.good() )
    surf.exportSurfaceAs3DOFF( out );
  out.close();
  trace.endBlock();


  return 0;
}
