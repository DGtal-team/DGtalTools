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
 * @file volscope.cpp
 * @ingroup Visualisation
 * @author David Coeurjolly (\c david.coeurjolly@cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS,
 * France
 *
 * @date 2023/12/22
 *
 * Vol visualization using polyscope.
 *
 * This file is part of the DGtal library.
 */

#include <vector>
#include <array>

#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>

#include <polyscope/polyscope.h>
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"
#include "polyscope/point_cloud.h"

#include "CLI11.hpp"

// Using standard 3D digital space.
using namespace DGtal;
using namespace Z3i;
typedef Shortcuts<Z3i::KSpace> SH3;


/**
 @page volscope volscope
 
 @brief Volumetric file visualization using polyscope
 @ingroup visualizationtools
 
 @b Usage: volscope [options] --input  \<fileName\>
 
 
 @b Allowed @b options @b are :
 @code
 Usage: ./volscope [OPTIONS] 1
 
 Positionals:
 1 TEXT:FILE REQUIRED                  Input vol file.
 
 Options:
 -h,--help                             Print this help message and exit
 -i,--input TEXT:FILE REQUIRED         Input vol file.
 --volumetricMode                      Activate the volumetric mode instead of the isosurface visualization.
 --point-cloud-only                    In the volumetric mode, visualize the vol file as a point cloud instead of an hex mesh (default: false)
-m,  --min INT                             For isosurface visualization and voxel filtering, specifies the threshold min (excluded) (default: 0).
 -M, --max INT                             For isosurface visualization and voxel filtering, specifies the threshold max (included) (default: 255).
 @endcode


 @image html volscope-surface.png "Default visualization mode." width=50%

 @image html volscope-volumetric.png "Volumetric visualization using polyscope slice planes."  width=50%




**/

int main(int argc, char**argv)
{
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  app.description("Vol file vizualization using polyscope");
  
  std::string inputFileName;
  app.add_option("-i,--input,1", inputFileName, "Input vol file.")->required()->check(CLI::ExistingFile);
  
  bool volumetricMode = false;
  app.add_flag("--volumetricMode",volumetricMode, "Activate the volumetric mode instead of the isosurface visualization.");

  bool pclOnly = false;
  app.add_flag("--point-cloud-only",pclOnly, "In the volumetric mode, visualize the vol file as a point cloud instead of an hex mesh (default: false)");
  
  int thresholdMin=0;
  app.add_option("--min,--thresholdMin,-m", thresholdMin, "For isosurface visualization and voxel filtering, specifies the threshold min (excluded) (default: 0).");
  int thresholdMax=255;
  app.add_option("--max, --thresholdMax,-M", thresholdMax, "For isosurface visualization and voxel filtering, specifies the threshold max (included) (default: 255).");
  
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  
  polyscope::options::programName = "volscope " + inputFileName +  " - (DGtalTools)";  
  polyscope::init();
  
  if(!volumetricMode)
  {
    trace.beginBlock("Loading vol");
    auto params = SH3::defaultParameters();
    params("thresholdMin",thresholdMin)("thresholdMax",thresholdMax);
    auto volimage        = SH3::makeBinaryImage(inputFileName, params );
    auto K               = SH3::getKSpace( volimage );
    auto surface         = SH3::makeLightDigitalSurface( volimage, K, params );
    auto primalSurface   = SH3::makePrimalSurfaceMesh(surface);
    trace.endBlock();
    
    trace.beginBlock("Creating polyscope surface");
    std::vector<std::vector<SH3::SurfaceMesh::Vertex>> faces;
    for(auto face= 0 ; face < primalSurface->nbFaces(); ++face)
      faces.push_back(primalSurface->incidentVertices( face ));
    polyscope::registerSurfaceMesh("Vol file", primalSurface->positions(), faces);
    trace.endBlock();
  }
  else
  {
    std::vector<RealPoint> vertexPos;
    std::vector<Point> pclPos;
    std::vector<int> pclData;
    std::vector<std::array<size_t,8>> hexIndices;
    std::vector<int> hexData;
    
    trace.beginBlock("Loading vol");
    auto volimage = SH3::makeGrayScaleImage(inputFileName);
    trace.endBlock();
    trace.beginBlock("Creating polyscope HexMesh/Point Cloud");
    auto dom = volimage->domain();
    auto W = dom.upperBound() - dom.lowerBound() + Point::diagonal(1);
    auto WW = W + Point::diagonal(1); //for dual grid
    trace.info()<<W<<std::endl;
    trace.info()<<dom<<std::endl;
    
    size_t cpt=0;
    unsigned char val;
    Point p;
    std::array<size_t,8> hex;
    for(auto k=0; k <= W[2]; ++k)
      for(auto j=0; j <= W[1]; ++j)
        for(auto i=0; i <= W[0]; ++i)
        {
          p=Point(i,j,k);
          if ((i<W[0]) && (j < W[1]) &&(k<W[2]))
            val = (*volimage)(p);
          
          if (pclOnly)
          {
            if ((p < W) && (val>thresholdMin) && (val <=thresholdMax))
            {
              pclPos.push_back(p);
              pclData.push_back(val);
            }
          }
          else
          {
            vertexPos.push_back(p+dom.lowerBound() -RealPoint::diagonal(0.5) );
            hex = { cpt, cpt +1 , cpt + 1 + WW[0] , cpt +WW[0] , cpt + WW[0]*WW[1], cpt +1 + WW[0]*WW[1], cpt + 1 + WW[0]*WW[1]+WW[0] , cpt + WW[0]*WW[1]+WW[0]};
            
            if (((i+1)< WW[0]) && ((j+1)< WW[1]) && ((k+1)< WW[2])&& (val>thresholdMin) && (val <=thresholdMax))
            {
              hexData.push_back(val);
              hexIndices.push_back(hex);
            }
            ++cpt;
          }
        }
    
    if (pclOnly)
    {
      auto ps = polyscope::registerPointCloud("Vol file", pclPos);
      ps->addScalarQuantity("values", pclData);
      trace.info()<<"Nb vertices ="<<vertexPos.size()<<std::endl;
      trace.info()<<"Nb hexes ="<<hexIndices.size()<<std::endl;
    }
    else
    {
      auto ps = polyscope::registerHexMesh("Vol file", vertexPos, hexIndices);
      ps->addCellScalarQuantity("values", hexData);
      trace.info()<<"Nb points ="<<pclPos.size()<<std::endl;
    }
    trace.endBlock();
  }
  polyscope::show();
  return 0; 
}
