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
 * @file 3dVolMarchingCubes.cpp
 * @ingroup volumetric
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/02/29
 *
 * An example file named 3dVolMarchingCubes.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
//! [3dVolMarchingCubes-basicIncludes]
#include <iostream>
#include <queue>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <DGtal/kernel/sets/SetPredicate.h>
#include <DGtal/io/readers/VolReader.h>
#include <DGtal/images/ImageSelector.h>
#include <DGtal/images/SimpleThresholdForegroundPredicate.h>
#include <DGtal/images/ImageLinearCellEmbedder.h>
#include <DGtal/shapes/Shapes.h>
#include <DGtal/kernel/CanonicEmbedder.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/topology/DigitalSurface.h>
#include <DGtal/topology/SetOfSurfels.h>
#include <DGtal/geometry/volumes/KanungoNoise.h>
//! [3dVolMarchingCubes-basicIncludes]


#include "CLI11.hpp"

///////////////////////////////////////////////////////////////////////////////

using namespace DGtal;
using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////

/**
 @page Doc3dVolMarchingCubes 3dVolMarchingCubes
 
 @brief Outputs the isosurface of the input volume  as an OFF file
 
 @b Usage: 3dVolMarchingCubes [-i \<fileName.vol\>] [-t \<threshold\>] [-a \<adjacency\>] [-o \<output.off\>]
 
 Outputs the isosurface of value \<threshold\> of the volume
 \<fileName.vol\> as an OFF file \<output.off\>. The \<adjacency\> (0/1)
 allows to choose between interior (6,18) and exterior (18,6)
 adjacency.
 
 @b Allowed @b options @b are :
 @code

 Positionals:
   1 TEXT:FILE REQUIRED                  the volume file (.vol).
   2 TEXT=marching-cubes.off             the output OFF file that represents the geometry of the isosurface

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         the volume file (.vol).
   -t,--threshold FLOAT=1                the value that defines the isosurface in the image (an integer between 0 and 255).
   -a,--adjacency UINT=0                 0: interior adjacency, 1: exterior adjacency
   -o,--output TEXT=marching-cubes.off   the output OFF file that represents the geometry of the isosurface
   -n,--noise FLOAT                      Kanungo noise level in ]0,1[. Note that only the largest connected component is considered and that no specific embedder is used.
   
 @endcode
 
 @b Example:
 
 @code
 $ 3dVolMarchingCubes -i $DGtal/examples/samples/lobster.vol -t 30
 # we invert the default normol orientation to improve display (-n option):
 $ meshViewer -i marching-cubes.off -n
 @endcode
 
 
 You should obtain such a result:
 @image html res3dVolMarchingCubes.png "Resulting visualization."
 
 You can test on other samples like http://www.tc18.org/code_data_set/3D_greyscale/bonsai.vol.bz2
 @code
 $ 3dVolMarchingCubes -i bonsai.vol  -t 80
 @endcode
 You should obtain such a result:
 @image html res3dVolMarchingCubes2.png "Resulting visualization."
 
 
 @see
 \ref 3dVolMarchingCubes.cpp
 
 */



int main( int argc, char** argv )
{
  //! [3dVolMarchingCubes-parseCommandLine]
 
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"marching-cubes.off"};
  double noise;
  double threshold {1.0};
  unsigned int intAdjacency = 0;
  
  app.description("Outputs the isosurface of value <threshold> of the volume <fileName.vol>  as an OFF file <output.off>. The <adjacency> (0/1) allows to choose between interior (6,18) and exterior (18,6) adjacency.");
  app.add_option("-i,--input,1", inputFileName, "the volume file (.vol)." )
  ->required()
  ->check(CLI::ExistingFile);
  app.add_option("--threshold,-t",threshold, "the value that defines the isosurface in the image (an integer between 0 and 255).", true);
  app.add_option("--adjacency,-a",intAdjacency, "0: interior adjacency, 1: exterior adjacency", true);
  app.add_option("-o,--output,2", outputFileName, "the output OFF file that represents the geometry of the isosurface", true );
  auto noiseOpt = app.add_option("--noise,-n", noise, "Kanungo noise level in ]0,1[. Note that only the largest connected component is considered and that no specific embedder is used.");
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  

  bool addNoise=false;
  if (noiseOpt->count() > 0 )
  {
    addNoise=true;
  }
  //! [3dVolMarchingCubes-parseCommandLine]
  
  //! [3dVolMarchingCubes-readVol]
  trace.beginBlock( "Reading vol file into an image." );
  typedef ImageSelector < Domain, int>::Type Image;
  Image image = VolReader<Image>::importVol(inputFileName);
  
  typedef functors::SimpleThresholdForegroundPredicate<Image> ThresholdedImage;
  ThresholdedImage thresholdedImage( image, threshold );
  trace.endBlock();
  //! [3dVolMarchingCubes-readVol]
  
  
  //! [3dVolMarchingCubes-KSpace]
  trace.beginBlock( "Construct the Khalimsky space from the image domain." );
  KSpace ks;
  bool space_ok =
  ks.init( image.domain().lowerBound(), image.domain().upperBound(), true );
  if (!space_ok)
  {
    trace.error() << "Error in the Khamisky space construction."<<std::endl;
    return 2;
  }
  trace.endBlock();
  //! [3dVolMarchingCubes-KSpace]
  
  //! [3dVolMarchingCubes-SurfelAdjacency]
  typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
  MySurfelAdjacency surfAdj( intAdjacency ); // interior in all directions.
  //! [3dVolMarchingCubes-SurfelAdjacency]
  
  //! [3dVolMarchingCubes-ExtractingSurface]
  trace.beginBlock( "Extracting boundary by scanning the space. " );
  typedef KSpace::SurfelSet SurfelSet;
  typedef SetOfSurfels< KSpace, SurfelSet > MySetOfSurfels;
  typedef DigitalSurface< MySetOfSurfels > MyDigitalSurface;
  MySetOfSurfels theSetOfSurfels( ks, surfAdj );
  
  if (addNoise)
  {
    trace.info()<<"Adding some noise."<<std::endl;
    KanungoNoise<ThresholdedImage, Z3i::Domain> kanungo(thresholdedImage, image.domain(), noise);
    std::vector< std::vector<SCell > > vectConnectedSCell;
    trace.info()<<"Extracting all connected components."<<std::endl;
    Surfaces<KSpace>::extractAllConnectedSCell( vectConnectedSCell, ks, surfAdj,
                                               kanungo, false);
    if( vectConnectedSCell.size() == 0 )
    {
      trace.error()<< "No surface component exists. Please check the vol file --min and --max parameters." << std::endl;
      return 3;
    }
    
    trace.info()<<vectConnectedSCell.size()<<" components."<<std::endl;
    
    trace.info()<<"Extracting the largest one."<<std::endl;
    int cc_max_size_idx = -1;
    auto it_max         = std::max_element( vectConnectedSCell.begin(), vectConnectedSCell.end(),
                                           [] (std::vector<SCell >& v1, std::vector<SCell >& v2)
                                           { return v1.size() < v2.size(); } );
    cc_max_size_idx = std::distance( vectConnectedSCell.begin(), it_max );
    theSetOfSurfels.surfelSet().insert( vectConnectedSCell[ cc_max_size_idx ].begin(),
                                       vectConnectedSCell[ cc_max_size_idx ].end() );
  }
  else
    Surfaces<KSpace>::sMakeBoundary(
    theSetOfSurfels.surfelSet(), ks, thresholdedImage,
    image.domain().lowerBound(), image.domain().upperBound() );

  MyDigitalSurface digSurf( theSetOfSurfels );
  trace.info() << "Digital surface has " << digSurf.size() << " surfels."
  << std::endl;
  trace.endBlock();
  //! [3dVolMarchingCubes-ExtractingSurface]
  
  
  //! [3dVolMarchingCubes-makingOFF]
  trace.beginBlock( "Making OFF surface. " );
  // Describes how voxels are embedded into Euclidean space.
  typedef CanonicEmbedder< Space > MyEmbedder;
  // Describes how the centroid surface elements is placed in-between embedded voxels.
  typedef ImageLinearCellEmbedder< KSpace, Image, MyEmbedder > CellEmbedder;
  CellEmbedder cellEmbedder;
  MyEmbedder trivialEmbedder;
  
  // The +0.5 is to avoid isosurface going exactly through a voxel
  // center, especially for binary volumes.
  
  //Making sure that we probe inside the domain.
  Image largerImage( Domain( image.domain().lowerBound() - KSpace::Vector::diagonal(2),
                            image.domain().upperBound() + KSpace::Vector::diagonal(2)));
  for(auto p: image.domain())
    largerImage.setValue( p, image(p));
  
  std::ofstream out( outputFileName.c_str() );
  if (addNoise)
  {
    if ( out.good() )
      digSurf.exportSurfaceAs3DOFF( out );
    else
    {
      trace.error()<<"IO error while opening the output file"<<std::endl;
      exit(2);
    }
  }
  else
  {
    cellEmbedder.init( ks, image, trivialEmbedder,
                      ( (double) threshold ) + 0.5 );
    
    if ( out.good() )
      digSurf.exportEmbeddedSurfaceAs3DOFF( out, cellEmbedder );
    else
    {
      trace.error()<<"IO error while opening the output file"<<std::endl;
      exit(2);
    }
  }
  
  out.close();
  trace.endBlock();
  //! [3dVolMarchingCubes-makingOFF]
  return 0;
}

