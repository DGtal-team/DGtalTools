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
 * @file mesh2vol.cpp
 * @ingroup converters
 *
 * @date 2018/01/11
 *
 *
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <chrono>
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/MeshVoxelizer.h"
#include "DGtal/kernel/sets/CDigitalSet.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/Display3D.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;


/**
 @page mesh2vol
 @brief Convert a mesh file into a 26-separated or 6-separated voxelization in a given resolution grid.

@b Usage: mesh2vol [input]

@b Allowed @b options @b are:

@code
ositionals:
  1 TEXT:FILE REQUIRED                  mesh file (.off).
  2 TEXT=result.vol                     filename of ouput volumetric file (vol, pgm3d, ...).

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         mesh file (.off).
  -o,--output TEXT=result.vol           filename of ouput volumetric file (vol, pgm3d, ...).
  -m,--margin UINT                      add volume margin around the mesh bounding box.
  -s,--separation UINT:{6,26}=6         voxelization 6-separated or 26-separated.
  -r,--resolution UINT=128              digitization domain size (e.g. 128). The mesh will be scaled such that its bounding box maps to [0,resolution)^3.
@endcode

@b Example:
@code
  $ mesh2vol -i ${DGtal}/examples/samples/tref.off --separation 26 --resolution 256 -o output.vol
@endcode

@see mesh2vol.cpp

*/

template< unsigned int SEP >
void voxelizeAndExport(const std::string& inputFilename,
                       const std::string& outputFilename,
                       const unsigned int resolution,
                       const unsigned int margin)
{
  using Domain   = Z3i::Domain;
  using PointR3  = Z3i::RealPoint;
  using PointZ3  = Z3i::Point;

  trace.beginBlock("Preparing the mesh");
  trace.info() << "Reading input file: " << inputFilename;
  Mesh<PointR3> inputMesh;
  MeshReader<PointR3>::importOFFFile(inputFilename.c_str(), inputMesh);
  trace.info() << " [done]" << std::endl;
  const std::pair<PointR3, PointR3> bbox = inputMesh.getBoundingBox();
  trace.info()<< "Mesh bounding box: "<<bbox.first <<" "<<bbox.second<<std::endl;

  const double smax = (bbox.second - bbox.first).max();
  const double factor = resolution / smax;
  const PointR3 translate = -bbox.first;
  trace.info() << "Scale = "<<factor<<" translate = "<<translate<<std::endl;
  for(auto it = inputMesh.vertexBegin(), itend = inputMesh.vertexEnd();
      it != itend; ++it)
  {
    //scale + translation
    *it += translate;
    *it *= factor;
  }
  trace.endBlock();
  
  trace.beginBlock("Voxelization");
  trace.info() << "Voxelization " << SEP << "-separated ; " << resolution << "^3 ";
  Domain aDomain(PointZ3().diagonal(-margin), PointZ3().diagonal(resolution+margin));
  
  //Digitization step
  Z3i::DigitalSet mySet(aDomain);
  MeshVoxelizer<Z3i::DigitalSet, SEP> aVoxelizer;
  aVoxelizer.voxelize(mySet, inputMesh, 1.0);
  trace.info() << " [done] " << std::endl;
  trace.endBlock();
  
  trace.beginBlock("Exporting");
  // Export the digital set to a vol file
  trace.info()<<aDomain<<std::endl;
  ImageContainerBySTLVector<Domain, unsigned char> image(aDomain);
  for(auto p: mySet)
    image.setValue(p, 128);
  image >> outputFilename.c_str();
  trace.endBlock();
}

int main( int argc, char** argv )
{ 

// parse command line using CLI ----------------------------------------------
   CLI::App app;
   std::string inputFileName;
   std::string outputFileName {"result.vol"};
   unsigned int margin  {0};
   unsigned int separation {6};
   unsigned int resolution {128};
 
   app.description("Convert a mesh file into a 26-separated or 6-separated volumetric voxelization in a given resolution grid. \n Example:\n mesh2vol -i ${DGtal}/examples/samples/tref.off -o output.vol --separation 26 --resolution 256 ");
   
   app.add_option("-i,--input,1", inputFileName, "mesh file (.off)." )
     ->required()
     ->check(CLI::ExistingFile);
   app.add_option("-o,--output,2", outputFileName, "filename of ouput volumetric file (vol, pgm3d, ...).",true);
   app.add_option("-m,--margin", margin, "add volume margin around the mesh bounding box.");
   app.add_option("-s,--separation", separation, "voxelization 6-separated or 26-separated.", true)
     -> check(CLI::IsMember({6, 26}));
   app.add_option("-r,--resolution", resolution,"digitization domain size (e.g. 128). The mesh will be scaled such that its bounding box maps to [0,resolution)^3.", true);
   
   
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------    
 
  if (separation==6)
    voxelizeAndExport<6>(inputFileName, outputFileName, resolution, margin);
  else
    voxelizeAndExport<26>(inputFileName, outputFileName, resolution, margin);    
 return EXIT_SUCCESS;
}

