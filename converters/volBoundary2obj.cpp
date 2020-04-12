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
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

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
namespace po = boost::program_options;

/**
 @page volBoundary2obj volBoundary2obj
 @brief  Extracts digital points from 3d vol files.
 
 @b Usage: volBoundary2obj [input] [output]
 
 @b Allowed @b options @b are:
 
 @code
  -h [ --help ]                    display this message
  -i [ --input ] arg               vol file (.vol, .longvol .p3d, .pgm3d and if
                                   WITH_ITK is selected: dicom, dcm, mha, mhd).
                                   For longvol, dicom, dcm, mha or mhd formats,
                                   the input values are linearly scaled between
                                   0 and 255.
  -o [ --output ] arg              output obj file (.obj)
  -m [ --thresholdMin ] arg (=0)   threshold min (excluded) to define binary 
                                   shape
  -M [ --thresholdMax ] arg (=255) threshold max (included) to define binary 
                                   shape
  --customDiffuse arg              set the R, G, B, A components of the diffuse
                                   colors of the mesh faces.
  --rescaleInputMin arg (=0)       min value used to rescale the input 
                                   intensity (to avoid basic cast into 8  bits 
                                   image).
  --rescaleInputMax arg (=255)     max value used to rescale the input 
                                   intensity (to avoid basic cast into 8 bits 
                                   image).
  -t [ --triangulatedSurface ]     save the dual triangulated surface instead 
                                   instead the default digital surface.

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
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
  ("help,h", "display this message")
  ("input,i", po::value<std::string>(), "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
  ("output,o", po::value<std::string>(), "output obj file (.obj)" )
  ("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min (excluded) to define binary shape" )
  ("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max (included) to define binary shape" )
  ("customDiffuse",po::value<std::vector<unsigned int> >()->multitoken(), "set the R, G, B, A components of the diffuse colors of the mesh faces." )
  ("rescaleInputMin", po::value<DGtal::int64_t>()->default_value(0), "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).")
  ("rescaleInputMax", po::value<DGtal::int64_t>()->default_value(255), "max value used to rescale the input intensity (to avoid basic cast into 8 bits image).")
  ("triangulatedSurface,t","save the dual triangulated surface instead instead the default digital surface.");
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);
  if( !parseOK || vm.count("help")||argc<=1)
  {
    std::cout << "Usage: " << argv[0] << " -i [input] -o [output]\n"
    << "Export the boundary of a volume file to OBJ format. By default the resulting mesh is defined from the surfles of the surface elements, a triangulated (dual)"<< endl
    << general_opt << "\n";
    return 0;
  }
  
  if(! vm.count("input"))
  {
    trace.error() << " The file name was defined" << endl;
    return 0;
  }
  
  if(! vm.count("output"))
  {
    trace.error() << " The output filename was defined" << endl;
    return 0;
  }
  
  
  // Using standard 3D digital space.
  typedef Shortcuts<Z3i::KSpace> SH3;
  typedef ShortcutsGeometry<Z3i::KSpace> SHG3;
  auto params = SH3::defaultParameters() | SHG3::defaultParameters();
  
  
  
  string inputFilename = vm["input"].as<std::string>();
  int thresholdMin = vm["thresholdMin"].as<int>();
  int thresholdMax = vm["thresholdMax"].as<int>();
  
  DGtal::int64_t rescaleInputMin = vm["rescaleInputMin"].as<DGtal::int64_t>();
  DGtal::int64_t rescaleInputMax = vm["rescaleInputMax"].as<DGtal::int64_t>();
  
  trace.beginBlock( "Loading file.." );
  auto gimage    = SH3::makeGrayScaleImage(inputFilename);
  typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
  RescalFCT f (rescaleInputMin, rescaleInputMax,0, 255);
  for (const auto &i: gimage->domain()) gimage->setValue(i, f((*gimage)(i)));
  
  
  auto bimage = SH3::makeBinaryImage(gimage,params( "thresholdMin", rescaleInputMin ) |
                                     params( "thresholdMax", rescaleInputMax ) );
  
  
  trace.info() << "Image loaded: "<<gimage<< std::endl;
  trace.endBlock();
  params( "faceSubdivision", "Centroid" )( "surfelAdjacency", 1);
  
  auto K         = SH3::getKSpace( bimage);
  
  
  string outputFilename = vm["output"].as<std::string>();
  bool customDiffuse =  vm.count("customDiffuse");
  
  SH3::Color cD ( 30,30,30 );
  
  if(customDiffuse)
  {
    std::vector<unsigned int > vectCol = vm["customDiffuse"].as<std::vector<unsigned int> >();
    if(vectCol.size()!=4)
    {
      trace.error() << "colors specification should contain R,G,B and Alpha values"<< std::endl;
    }
    cD.setRGBi(vectCol[0], vectCol[1], vectCol[2], vectCol[3]);
  }
  
  auto surface = SH3::makeLightDigitalSurface( bimage, K, params( "thresholdMin", thresholdMin ) |
  params( "thresholdMax", thresholdMax ) );
  
  if (vm.count("triangulatedSurface"))
  {
    auto tr = SH3::makeTriangulatedSurface(surface);
    auto ok  = SH3::saveOBJ( tr, SH3::RealVectors(), SH3::Colors(), outputFilename, cD );
    return ok ? EXIT_SUCCESS : EXIT_FAILURE ;
  }
  else
  {
    auto ok  = SH3::saveOBJ( surface, SH3::RealVectors(), SH3::Colors(), outputFilename,cD);
    return  ok ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
