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
 * @file volBoundary2obj.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * @author Bertrand Kerautret (\c bertrand.kerautret@liris.cnrs.fr )
 *
 * @date 2013/11/15
 *
 * A tool file named volBoundary2obj.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <set>
#include "CLI11.hpp"

#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/readers/GenericReader.h"
#ifdef WITH_ITK
#include "DGtal/io/readers/DicomReader.h"
#endif
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/helpers/Shortcuts.h"
#include "DGtal/helpers/ShortcutsGeometry.h"

using namespace std;
using namespace DGtal;
//using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////


/**
 @page volBoundary2obj volBoundary2obj
 @brief Export the boundary of a volume file to OBJ format. By default the resulting mesh is defined from the surfels of the surface elements, a triangulated (dual)
 
 @b Usage: converters/volBoundary2obj [OPTIONS] 1 [2]
 
 @b Allowed @b options @b are:
 
 @code


Positionals:
  1 TEXT:FILE REQUIRED                  vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
  2 TEXT                                output file (.obj or .off).

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
  -o,--output TEXT                      output file (.obj or .off).
  -m,--thresholdMin INT=128             threshold min (excluded) to define binary shape.
  -M,--thresholdMax INT=255             threshold max (included) to define binary shape.
  --rescaleInputMin INT=0               min value used to rescale the input intensity (to avoid basic cast into 8  bits image).
  --rescaleInputMax INT=255             max value used to rescale the input intensity (to avoid basic cast into 8  bits image).
  -c,--customDiffuse UINT=[230,230,230,255] x 4
                                        set the R, G, B, A components of the diffuse colors of the mesh faces.
  -t,--triangulatedSurface              save the dual triangulated surface instead instead the default digital surface.


 @endcode
 
 @b Example:
 @code
 $ volBoundary2obj -i $DGtal/examples/samples/lobster.vol -m 80 -o out.obj
 @endcode
 
 You should obtain such a visualization:
 @image html resVolBoundary2obj.png "resulting visualisation."
 
 @see
 @ref volBoundary2obj.cpp
 
 */


int main( int argc, char** argv )
{

  CLI::App app;
  
  
  // Using standard 3D digital space.
  typedef Shortcuts<Z3i::KSpace> SH3;
  typedef ShortcutsGeometry<Z3i::KSpace> SHG3;
  auto params = SH3::defaultParameters() | SHG3::defaultParameters();
  
  int thresholdMin {128};
  int thresholdMax {255};  
  string inputFilename;
  DGtal::int64_t rescaleInputMin {0};
  DGtal::int64_t rescaleInputMax {255};
  std::vector<unsigned int > vectCol=  {230, 230, 230, 255};
  bool triangulatedSurface {false};
  std::string outputFilename = "result.obj";
 
  app.description("Export the boundary of a volume file to OBJ format. By default the resulting mesh is defined from the surfels of the surface elements, a triangulated (dual)");
  
  app.add_option("-i,--input,1", inputFilename, "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("--output,-o,2",outputFilename ,"output file (.obj or .off).");
  app.add_option("--thresholdMin,-m", thresholdMin, "threshold min (excluded) to define binary shape.", true);
  app.add_option("--thresholdMax,-M", thresholdMax, "threshold max (included) to define binary shape.", true);
  app.add_option("--rescaleInputMin", rescaleInputMin, "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);
  app.add_option("--rescaleInputMax", rescaleInputMax, "max value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);
  app.add_option("--customDiffuse,-c", vectCol, "set the R, G, B, A components of the diffuse colors of the mesh faces.", true)
    ->expected(4);
  app.add_flag("--triangulatedSurface,-t", triangulatedSurface, "save the dual triangulated surface instead instead the default digital surface.");
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  


  typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
  RescalFCT f (rescaleInputMin, rescaleInputMax,0, 255);

  trace.beginBlock( "Loading file.." );
  SH3::GrayScaleImage image =
  GenericReader< SH3::GrayScaleImage >::importWithValueFunctor(inputFilename, f );
  
  auto gimage = CountedPtr<SH3::GrayScaleImage>( new SH3::GrayScaleImage( image ) );
  auto bimage = SH3::makeBinaryImage(gimage,params( "thresholdMin", thresholdMin )
                                                  ( "thresholdMax", thresholdMax ) );
  
  
  trace.info() << "Image loaded: "<<gimage<< std::endl;
  trace.endBlock();
  params( "faceSubdivision", "Centroid" )( "surfelAdjacency", 1);
  auto K         = SH3::getKSpace( bimage);

  SH3::Color cD (vectCol[0], vectCol[1], vectCol[2], vectCol[3]);
  
  
  auto surface = SH3::makeDigitalSurface( bimage, K );
  const std::string extension = outputFilename.substr( outputFilename.find_last_of(".") + 1 );  
  if (extension != "obj" && extension != "off")
  {
    trace.warning() << "File extension not recognized, saving by default in objg format"<< std::endl;
  }
  if (triangulatedSurface)
  {
    auto tr = SH3::makeTriangulatedSurface(surface);
    bool ok = true;
    if (extension!="off")
    {
      ok  = SH3::saveOBJ( tr, SH3::RealVectors(), SH3::Colors(), outputFilename, Color( 32, 32, 32 ), cD );
    }
    else 
    {
      ok  = SH3::saveOFF( tr, outputFilename, cD);
    }
    return ok ? EXIT_SUCCESS : EXIT_FAILURE ;
  }
  else
  {
    bool ok = true;
    if (extension!="off")
    {
      ok  = SH3::saveOBJ( surface, SH3::RealVectors(), SH3::Colors(), outputFilename, Color( 32, 32, 32 ), cD);
    } else
    {
      ok  = SH3::saveOFF( surface, outputFilename, cD);
    }
    return  ok ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}


