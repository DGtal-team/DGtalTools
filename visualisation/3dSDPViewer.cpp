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
 * @file 3dSDPViewer.cpp
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2014/04/01
 *
 * An simple viewer to display 3d SDP files (sequence of discrete points).
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/readers/TableReader.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/shapes/Mesh.h"
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"

#include "CLI11.hpp"



using namespace std;
using namespace DGtal;
using namespace Z3i;

typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;

/**
 @page Doc3DSDPViewer 3DSDPViewer
 
 @brief Displays a  sequence of 3d discrete points by using QGLviewer.

 @b Usage:  3dSDPViewer [OPTIONS] 1 [f] [lineSize] [filterVectors]

 @b Allowed @b options @b are :
 
 @code

 Positionals:
   1 TEXT:FILE REQUIRED                  input file: sdp (sequence of discrete points).

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         input file: sdp (sequence of discrete points).
   --SDPindex UINT x 3                   specify the sdp index (by default 0,1,2).
   -c,--pointColor INT x 4               set the color of  points: r g b a.
   -l,--lineColor INT x 4                set the color of line: r g b a.
   -m,--addMesh TEXT                     append a mesh (off/obj) to the point set visualization.
   --customColorMesh UINT x 4            set the R, G, B, A components of the colors of the mesh faces (mesh added with option --addMesh).
   --importColors                        import point colors from the input file (R G B colors should be by default at index 3, 4, 5).
   --setColorsIndex UINT x 3 Needs: --importColors
                                         customize the index of the imported colors in the source file (used by --importColor). By default the initial indexes are 3, 4, 5.
   --importColorLabels                   import color labels from the input file (label index  should be by default at index 3).
   --setColorLabelIndex UINT=3           customize the index of the imported color labels in the source file (used by -importColorLabels).
   -f,--filter FLOAT=100                 filter input file in order to display only the [arg] percent of the input 3D points (uniformly selected).
   --noPointDisplay                      usefull for instance to only display the lines between points.
   --drawLines                           draw the line between discrete points.
   -x,--scaleX FLOAT=1                   set the scale value in the X direction
   -y,--scaleY FLOAT=1                   set the scale value in the Y direction
   -z,--scaleZ FLOAT=1                   set the scale value in the Z direction
   --sphereResolution UINT=30            defines the sphere resolution (used when the primitive is set to the sphere).
   -s,--sphereRadius FLOAT=0.2           defines the sphere radius (used when the primitive is set to the sphere).
   --sphereRadiusFromInput               takes, as sphere radius, the 4th field of the sdp input file.
   --lineSize FLOAT=0.2                  defines the line size (used when the --drawLines or --drawVectors option is selected).
   -p,--primitive TEXT:{voxel,glPoints,sphere}=voxel
                                         set the primitive to display the set of points.
   -v,--drawVectors TEXT                 SDP vector file: draw a set of vectors from the given file (each vector are determined by two consecutive point given, each point represented by its coordinates on a single line.
   -u,--unitVector FLOAT=0               specifies that the SDP vector file format (of --drawVectors option) should be interpreted as unit vectors (each vector position is be defined from the input point (with input order) with a constant norm defined by [arg]).
   --filterVectors FLOAT=100             filters vector input file in order to display only the [arg] percent of the input vectors (uniformly selected, to be used with option --drawVectors else no effect).
   --interactiveDisplayVoxCoords          by using this option the coordinates can be displayed after selection (shift+left click on voxel).
  
@endcode


 @b Basic @b example: 

 You can display a set of 3D points with sphere primitive and lines:
 @code
 $  3DSDPViewer -i $DGtal/tests/samples/sinus3D.dat -p sphere -s 0.3 --drawLines --lineSize 5 
 @endcode

 You should obtain such a result:

 @image html res3DSDPViewer.png "Resulting visualization."


@b Example @b with @b interactive @b selection :

This tool can be useful to recover coordinates from a set of voxels. To do it, you have to add the option allowing to activate the interactive selection (with --interactiveDisplayVoxCoords), for instance if you apply:
@code 
$ 3dSDPViewer -i $DGtal/tests/samples/sinus3D.dat  --interactiveDisplayVoxCoords
 
@endcode 
you should be able to select a voxel by using the SHIFT key and by clicking on a voxel:
  @image html res3DSDPViewerInteractive.png " "

 
@b Visualization @b of @b large @b  point @b set

If you need to display an important number of points, you can use the primitive @e glPoints instead @e voxel or @e sphere (-p glPoints). You will obtain such type of a visualization:

  @image html res3DSDPViewerGLPoints.png " "



 @see
 @ref 3DSDPViewer.cpp

 */



// call back function to display voxel coordinates
int
displayCoordsCallBack( void* viewer, int name, void* data )
{
  vector<Z3i::RealPoint> *vectVoxels = (vector<Z3i::RealPoint> *) data;
  std::stringstream ss;
  ss << "Selected voxel: (" << (*vectVoxels)[name][0] <<  ", ";
  ss << (*vectVoxels)[name][1] <<  ", ";
  ss << (*vectVoxels)[name][2] <<  ") ";
  ((Viewer *) viewer)->displayMessage(QString(ss.str().c_str()), 100000);
  
  return 0;
}


int main( int argc, char** argv )
{
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::vector<unsigned int > vectIndexSDP {0,1,2};
  Color lineColor(100, 100, 250);
  Color pointColor(250, 250, 250);
  std::vector<int> vectColorPt;
  std::vector<int> vectColorLine;
  std::string meshName;
  std::string typePrimitive {"voxel"};
  std::string vectorsFileName;
  std::vector<double> vectSphereRadius;
  std::vector<unsigned int> vectColMesh;
  std::vector<unsigned int> vectIndexColorImport {3,4,5};
  bool importColorLabels {false};
  bool noPointDisplay {false};
  bool drawLines {false};
  bool sphereRadiusFromInput {true};
  bool importColors {false};
  bool interactiveDisplayVoxCoords {false};
  float sx {1.0};
  float sy {1.0};
  float sz {1.0};
  unsigned int sphereResolution {30};
  double sphereRadius {0.2};
  double lineSize {0.2};
  double filterValue {100.0};
  double constantNorm {0.0};
  double percentageFilterVect {100.0};

  unsigned int colorLabelIndex = 3;
  
  app.description("Display sequence of 3d discrete points by using QGLviewer.");
  app.add_option("-i,--input,1", inputFileName, "input file: sdp (sequence of discrete points)." )
  ->required()
  ->check(CLI::ExistingFile);
  app.add_option("--SDPindex",vectIndexSDP, "specify the sdp index (by default 0,1,2).")
  ->expected(3);
  app.add_option("--pointColor,-c",vectColorPt, "set the color of  points: r g b a.")
  ->expected(4);
  app.add_option("--lineColor,-l",vectColorLine, "set the color of line: r g b a.")
  ->expected(4);
  app.add_option("--addMesh,-m",meshName,  "append a mesh (off/obj) to the point set visualization.");
  auto optMesh = app.add_option("--customColorMesh",vectColMesh, "set the R, G, B, A components of the colors of the mesh faces (mesh added with option --addMesh).")
  ->expected(4);
  auto importColOpt = app.add_flag("--importColors",importColors, "import point colors from the input file (R G B colors should be by default at index 3, 4, 5).");
  app.add_option("--setColorsIndex", vectIndexColorImport,"customize the index of the imported colors in the source file (used by --importColor). By default the initial indexes are 3, 4, 5.")
  ->expected(3)
  ->needs(importColOpt);

  app.add_flag("--importColorLabels", importColorLabels,"import color labels from the input file (label index  should be by default at index 3)." );
  app.add_option("--setColorLabelIndex", colorLabelIndex, "customize the index of the imported color labels in the source file (used by -importColorLabels).", true);
  app.add_option("--filter,-f", filterValue, "filter input file in order to display only the [arg] percent of the input 3D points (uniformly selected).", true);
  app.add_flag("--noPointDisplay",noPointDisplay,  "usefull for instance to only display the lines between points." );
  app.add_flag("--drawLines", drawLines, "draw the line between discrete points." );
  app.add_option("--scaleX,-x", sx, "set the scale value in the X direction", true );
  app.add_option("--scaleY,-y", sy, "set the scale value in the Y direction", true );
  app.add_option("--scaleZ,-z", sy, "set the scale value in the Z direction", true );
  app.add_option("--sphereResolution",sphereResolution, "defines the sphere resolution (used when the primitive is set to the sphere).", true );
  app.add_option("-s,--sphereRadius", sphereRadius, "defines the sphere radius (used when the primitive is set to the sphere).", true);
  app.add_flag("--sphereRadiusFromInput", sphereRadiusFromInput, "takes, as sphere radius, the 4th field of the sdp input file.");
  app.add_option("--lineSize", lineSize, "defines the line size (used when the --drawLines or --drawVectors option is selected).", true);
  app.add_option("--primitive,-p",typePrimitive, "set the primitive to display the set of points.", true )
     -> check(CLI::IsMember({"voxel", "glPoints", "sphere"}));
  app.add_option("--drawVectors,-v",vectorsFileName, "SDP vector file: draw a set of vectors from the given file (each vector are determined by two consecutive point given, each point represented by its coordinates on a single line.");
  
  app.add_option("--unitVector,-u",  constantNorm, "specifies that the SDP vector file format (of --drawVectors option) should be interpreted as unit vectors (each vector position is be defined from the input point (with input order) with a constant norm defined by [arg]).", true );
  
  app.add_option("--filterVectors", percentageFilterVect, "filters vector input file in order to display only the [arg] percent of the input vectors (uniformly selected, to be used with option --drawVectors else no effect). ", true );
  app.add_flag("--interactiveDisplayVoxCoords", interactiveDisplayVoxCoords, " by using this option the coordinates can be displayed after selection (shift+left click on voxel).");
  
     
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  if (vectColorLine.size() == 4){
    lineColor.setRGBi(vectColorLine[0], vectColorLine[1], vectColorLine[2], vectColorLine[3]);
  }
  if (vectColorPt.size() == 4){
    pointColor.setRGBi(vectColorPt[0], vectColorPt[1], vectColorPt[2], vectColorPt[3]);
  }

  
  QApplication application(argc,argv);
  
 
  bool useUnitVector = constantNorm != 0.0;
  
  typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;
  Z3i::KSpace K;
  Viewer viewer( K );
  viewer.setWindowTitle("3dSPD Viewer");
  viewer.show();
  viewer.setGLScale(sx, sy, sz);
  viewer.myGLLineMinWidth = lineSize;
  viewer << CustomColors3D(pointColor, pointColor);
  
  
  // Get vector of colors if imported.
  std::vector<Color> vectColors;
  if(importColors)
  {
    std::vector<unsigned int> r = TableReader<unsigned int>::getColumnElementsFromFile(inputFileName,vectIndexColorImport[0]);
    std::vector<unsigned int> g = TableReader<unsigned int>::getColumnElementsFromFile(inputFileName,vectIndexColorImport[1]);
    std::vector<unsigned int> b = TableReader<unsigned int>::getColumnElementsFromFile(inputFileName,vectIndexColorImport[2]);
    for (unsigned int i = 0; i<r.size(); i++){
      vectColors.push_back(Color(r[i], g[i], b[i]));
    }
  }

  // Get vector of colors if imported.
  std::vector< int> vectColorLabels;
  unsigned int maxLabel = 1;
  if(importColorLabels)
  {
    vectColorLabels = TableReader< int>::getColumnElementsFromFile(inputFileName, colorLabelIndex);
    maxLabel = *(std::max_element(vectColorLabels.begin(), vectColorLabels.end()));
  }
  HueShadeColorMap<unsigned int> aColorMap(0, maxLabel);
  
  
  if(sphereRadiusFromInput)
  {
    vectSphereRadius = TableReader<double>::getColumnElementsFromFile(inputFileName,3);
  }
  
  vector<Z3i::RealPoint> vectVoxels;
  vectVoxels = PointListReader<Z3i::RealPoint>::getPointsFromFile(inputFileName, vectIndexSDP);
  
  int name = 0;
  if(!noPointDisplay){
    if (typePrimitive == "glPoints")
      {
        viewer.setUseGLPointForBalls(true);
      }

    int step = max(1, (int) (100/filterValue));
    for(unsigned int i=0;i< vectVoxels.size(); i=i+step){
      if(importColors)
        {
          Color col = vectColors[i];
          viewer.setFillColor(col);
        }
      else if(importColorLabels)
        {
          unsigned int index = vectColorLabels[i];
          Color col = aColorMap(index);
          viewer.setFillColor(col);
        }
      
      if(typePrimitive=="voxel" ){
        if (interactiveDisplayVoxCoords)
        {
          viewer << SetName3D( name++ ) ;
        }
        viewer << Z3i::Point((int)vectVoxels.at(i)[0],
                             (int)vectVoxels.at(i)[1],
                             (int)vectVoxels.at(i)[2]);
      }
      else
        {
          viewer.addBall(vectVoxels.at(i), sphereRadius, sphereResolution);
        }
    }
    
    viewer << CustomColors3D(lineColor, lineColor);
    if(drawLines)
    {
      for(unsigned int i=1;i< vectVoxels.size(); i++)
      {
        viewer.addLine(vectVoxels.at(i-1), vectVoxels.at(i), lineSize);
      }
    }
    
    if(vectorsFileName != "")
    {
      std::vector<Z3i::RealPoint> vectorsPt = PointListReader<Z3i::RealPoint>::getPointsFromFile(vectorsFileName);
      if (vectorsPt.size()%2==1)
      {
        trace.info()<<"Warning the two set of points doesn't contains the same number of points, some vectors will be skipped." << std::endl;
      }

      double percentage = percentageFilterVect;
      int step = max(1, (int) (100/percentage));
      
      if(useUnitVector)
      {
        for(unsigned int i =0; i< std::min(vectVoxels.size(), vectorsPt.size()); i=i+2*step)
        {
          viewer.addLine(vectVoxels.at(i), vectVoxels.at(i)+vectorsPt.at(i)*constantNorm, lineSize);
        }
      }
      else
      {
        for(unsigned int i =0; i<vectorsPt.size()-1; i=i+2*step)
        {
          viewer.addLine(vectorsPt.at(i),vectorsPt.at(i+1), lineSize);
        }
      }
    }
    if(meshName != "")
    {
      bool customColorMesh =  vectColMesh.size() == 4;
      if(customColorMesh)
      {
        viewer.setFillColor(DGtal::Color(vectColMesh[0], vectColMesh[1], vectColMesh[2], vectColMesh[3]));
      }
      Mesh<Z3i::RealPoint> mesh(!customColorMesh);
      mesh << meshName ;
      viewer << mesh;
    }
    if (interactiveDisplayVoxCoords)
    {
      viewer << SetSelectCallback3D( displayCoordsCallBack,  &vectVoxels,  0, vectVoxels.size()-1 );
    }
    
    viewer << Viewer3D<>::updateDisplay;
    return application.exec();
  }
}

