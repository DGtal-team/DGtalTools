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


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


using namespace std;
using namespace DGtal;
using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;
typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;


/**
 @page Doc3DSDPViewer 3DSDPViewer
 
 @brief Displays a  sequence of 3d discrete points by using QGLviewer.

 @b Usage:  3DSDPViewer [options] input

 @b Allowed @b options @b are :
 
 @code
 -h [ --help ]                         display this message
  -i [ --input ] arg                    input file: sdp (sequence of discrete 
                                        points)
  --SDPindex arg                        specify the sdp index (by default 
                                        0,1,2).
  -c [ --pointColor ] arg               set the color of  points: r g b a 
  -l [ --lineColor ] arg                set the color of line: r g b a 
  -m [ --addMesh ] arg                  append a mesh (off/obj) to the point 
                                        set visualization.
  --customColorMesh arg                 set the R, G, B, A components of the 
                                        colors of the mesh faces (mesh added 
                                        with option --addMesh). 
  --importColors                        import point colors from the input file
                                        (R G B colors should be by default at 
                                        index 3, 4, 5).
  --importColorLabels                   import color labels from the input file
                                        (label index  should be by default at 
                                        index 3).
  --setColorsIndex arg                  customize the index of the imported 
                                        colors in the source file (used by 
                                        -importColor).
  --setColorLabelIndex arg (=3)         customize the index of the imported 
                                        color labels in the source file (used 
                                        by -importColorLabels).
  -f [ --filter ] arg (=100)            filter input file in order to display 
                                        only the [arg] percentage of the input 3D
                                        points (uniformly selected).
  --noPointDisplay                      usefull for instance to only display 
                                        the lines between points.
  --drawLines                           draw the line between discrete points.
  -x [ --scaleX ] arg (=1)              set the scale value in the X direction 
                                        (default 1.0)
  -y [ --scaleY ] arg (=1)              set the scale value in the Y direction 
                                        (default 1.0)
  -z [ --scaleZ ] arg (=1)              set the scale value in the Z direction 
                                        (default 1.0)
  --sphereResolution arg (=30)          defines the sphere resolution (used 
                                        when the primitive is set to the 
                                        sphere). (default resolution: 30)
  -s [ --sphereRadius ] arg (=0.20000000000000001)
                                        defines the sphere radius (used when 
                                        the primitive is set to the sphere). 
                                        (default value 0.2)
  --sphereRadiusFromInput               takes, as sphere radius, the 4th field 
                                        of the sdp input file.
  --lineSize arg (=0.20000000000000001) defines the line size (used when the 
                                        --drawLines or --drawVectors option is 
                                        selected). (default value 0.2))
  -p [ --primitive ] arg (=voxel)       set the primitive to display the set of
                                        points (can be sphere, voxel (default),
                                        or glPoints (opengl points).
  -v [ --drawVectors ] arg              SDP vector file: draw a set of vectors 
                                        from the given file (each vector are 
                                        determined by two consecutive point 
                                        given, each point represented by its 
                                        coordinates on a single line.
  -u [ --unitVector ] arg (=1)          specifies that the SDP vector file 
                                        format (of --drawVectors option) should
                                        be interpreted as unit vectors (each 
                                        vector position is be defined from the 
                                        input point (with input order) with a 
                                        constant norm defined by [arg]).
  --filterVectors arg (=100)            filters vector input file in order to 
                                        display only the [arg] percentage of the 
                                        input vectors (uniformly selected, to 
                                        be used with option --drawVectors otherwise 
                                        no effect). 
  --interactiveDisplayVoxCoords         by using this option the pixel 
                                        coordinates can be displayed after 
                                        selection (shift+left click on voxel). 
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
$ 3DSDPViewer -i $DGtal/tests/samples/Al100.sdp  --interactiveDisplayVoxCoords 
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
  
  
  // parse command line ----------------------------------------------
  po::options_description general_opt(" Allowed options are");
  general_opt.add_options()
  ("help,h", "display this message")
  ("input,i", po::value<std::string>(), "input file: sdp (sequence of discrete points)" )
  ("SDPindex", po::value<std::vector <unsigned int> >()->multitoken(), "specify the sdp index (by default 0,1,2).")
  ("pointColor,c", po::value<std::vector <int> >()->multitoken(), "set the color of  points: r g b a " )
  ("lineColor,l",po::value<std::vector <int> >()->multitoken(), "set the color of line: r g b a " )
  ("addMesh,m", po::value<std::string>(), "append a mesh (off/obj) to the point set visualization.")
  ("customColorMesh",po::value<std::vector<unsigned int> >()->multitoken(), "set the R, G, B, A components of the colors of the mesh faces (mesh added with option --addMesh). " )   
  ("importColors", "import point colors from the input file (R G B colors should be by default at index 3, 4, 5).")
  ("importColorLabels", "import color labels from the input file (label index  should be by default at index 3).")
  ("setColorsIndex", po::value<std::vector<unsigned int> >()->multitoken(), "customize the index of the imported colors in the source file (used by -importColor).")
  ("setColorLabelIndex", po::value<unsigned int >()->default_value(3), "customize the index of the imported color labels in the source file (used by -importColorLabels).")
  ("filter,f",po::value<double>()->default_value(100.0), "filter input file in order to display only the [arg] percent of the input 3D points (uniformly selected)." )
  ("noPointDisplay", "usefull for instance to only display the lines between points.")
  ("drawLines", "draw the line between discrete points." )
  ("scaleX,x",  po::value<float>()->default_value(1.0), "set the scale value in the X direction (default 1.0)" )
  ("scaleY,y",  po::value<float>()->default_value(1.0), "set the scale value in the Y direction (default 1.0)" )
  ("scaleZ,z",  po::value<float>()->default_value(1.0), "set the scale value in the Z direction (default 1.0)")
  ("sphereResolution",  po::value<unsigned int>()->default_value(30), "defines the sphere resolution (used when the primitive is set to the sphere). (default resolution: 30)")
  ("sphereRadius,s",  po::value<double>()->default_value(0.2), "defines the sphere radius (used when the primitive is set to the sphere). (default value 0.2)")
  ("sphereRadiusFromInput", "takes, as sphere radius, the 4th field of the sdp input file.")
  ("lineSize",  po::value<double>()->default_value(0.2), "defines the line size (used when the --drawLines or --drawVectors option is selected). (default value 0.2))")
  ("primitive,p", po::value<std::string>()->default_value("voxel"), "set the primitive to display the set of points (can be sphere, voxel (default), or glPoints (opengl points).")
  ("drawVectors,v", po::value<std::string>(), "SDP vector file: draw a set of vectors from the given file (each vector are determined by two consecutive point given, each point represented by its coordinates on a single line.")
  ("unitVector,u", po::value<double>()->default_value(1.0), "specifies that the SDP vector file format (of --drawVectors option) should be interpreted as unit vectors (each vector position is be defined from the input point (with input order) with a constant norm defined by [arg]).")

  ("filterVectors",po::value<double>()->default_value(100.0), "filters vector input file in order to display only the [arg] percent of the input vectors (uniformly selected, to be used with option --drawVectors else no effect). " )
 
    ("interactiveDisplayVoxCoords", "by using this option the pixel coordinates can be displayed after selection (shift+left click on voxel)." );
  
  
  bool parseOK=true;
  bool cannotStart= false;
  
  
  
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.error()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);
  std::string typePrimitive;
  double sphereRadius = 0.2;
  std::vector<double> vectSphereRadius;
  unsigned int sphereResolution = vm["sphereResolution"].as<unsigned int>();
  double lineSize =0.2;
  bool useMultiRad = vm.count("sphereRadiusFromInput");
  
  Color lineColor(100, 100, 250);
  Color pointColor(250, 250, 250);
  if(parseOK)
  {
    typePrimitive = vm["primitive"].as<std::string>();
    sphereRadius = vm["sphereRadius"].as<double>();
    lineSize = vm["lineSize"].as<double>();
  }
  
  if (parseOK && typePrimitive !="voxel"
      && typePrimitive !="glPoints" 
      && typePrimitive  != "sphere" )
  {
    trace.error() << " The primitive should be sphere or voxel (primitive: "
    << typePrimitive << " not implemented)" << std::endl;
    cannotStart = true;
  }
  
  if(parseOK && vm.count("lineColor"))
  {
    std::vector<int> vcol= vm["lineColor"].as<std::vector<int > >();
    if(vcol.size()<4)
    {
      trace.error() << " Not enough parameter: color specification should contains four elements: red, green, blue and alpha values "
      << "(Option --lineColor ignored). "  << std::endl;
    }
    lineColor.setRGBi(vcol[0], vcol[1], vcol[2], vcol[3]);
  }
  if(parseOK && vm.count("pointColor"))
  {
    std::vector<int> vcol= vm["pointColor"].as<std::vector<int > >();
    if(vcol.size()<4)
    {
      trace.error() << " Not enough parameter: color specification should contains four elements: red, green, blue and alpha values "
      << "(Option --pointColor ignored)."  << std::endl;
    }
    pointColor.setRGBi(vcol[0], vcol[1], vcol[2], vcol[3]);
  }
  
  if( !parseOK || cannotStart ||  vm.count("help")||argc<=1)
  {
    trace.info() << "Usage: " << argv[0] << " [input]\n"
    << "Display sequence of 3d discrete points by using QGLviewer."
    << general_opt << "\n";
    return 0;
  }
  
  string inputFilename = vm["input"].as<std::string>();
  
  
  QApplication application(argc,argv);
  
  float sx = vm["scaleX"].as<float>();
  float sy = vm["scaleY"].as<float>();
  float sz = vm["scaleZ"].as<float>();

  bool importColorLabels = vm.count("importColorLabels");
  bool importColors = vm.count("importColors");
  bool interactiveDisplayVoxCoords = vm.count("interactiveDisplayVoxCoords");
  bool useUnitVector = vm.count("unitVector");
  double constantNorm = vm["unitVector"].as<double>();
  
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
  if(vm.count("importColors"))
  {
    std::vector<unsigned int > vectIndex;
    if(vm.count("setColorsIndex"))
    {
      vectIndex = vm["setColorsIndex"].as<std::vector<unsigned int > >();
      if(vectIndex.size()!=3)
      {
        trace.error() << "you need to specify the three indexes of color." << std::endl;
        return 0;
      }
    }
    else
    {
      vectIndex.push_back(3);
      vectIndex.push_back(4);
      vectIndex.push_back(5);
    }
    
    
    std::vector<unsigned int> r = TableReader<unsigned int>::getColumnElementsFromFile(inputFilename,vectIndex[0]);
    std::vector<unsigned int> g = TableReader<unsigned int>::getColumnElementsFromFile(inputFilename,vectIndex[1]);
    std::vector<unsigned int> b = TableReader<unsigned int>::getColumnElementsFromFile(inputFilename,vectIndex[2]);
    for (unsigned int i = 0; i<r.size(); i++){
      vectColors.push_back(Color(r[i], g[i], b[i]));
    }
  }

  // Get vector of colors if imported.
  std::vector< int> vectColorLabels;
  unsigned int maxLabel = 1;
  if(vm.count("importColorLabels"))
  {
    unsigned int index = vm["setColorLabelIndex"].as<unsigned int >();
    vectColorLabels = TableReader< int>::getColumnElementsFromFile(inputFilename,index);
    maxLabel = *(std::max_element(vectColorLabels.begin(), vectColorLabels.end()));
  }
  HueShadeColorMap<unsigned int> aColorMap(0, maxLabel);
  
  
  if(useMultiRad)
  {
    vectSphereRadius = TableReader<double>::getColumnElementsFromFile(inputFilename,3);
  }
  
  vector<Z3i::RealPoint> vectVoxels;
  if(vm.count("SDPindex"))
  {
    std::vector<unsigned int > vectIndex = vm["SDPindex"].as<std::vector<unsigned int > >();
    if(vectIndex.size()!=3)
    {
      trace.error() << "you need to specify the three indexes of vertex." << std::endl;
      return 0;
    }
    vectVoxels = PointListReader<Z3i::RealPoint>::getPointsFromFile(inputFilename, vectIndex);
    
  }else{
    vectVoxels = PointListReader<Z3i::RealPoint>::getPointsFromFile(inputFilename);
  }
  int name = 0;
  if(!vm.count("noPointDisplay")){
    if (typePrimitive == "glPoints")
      {
        viewer.setUseGLPointForBalls(true);
      }

    double percent = vm["filter"].as<double>();
    int step = max(1, (int) (100/percent));
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
    if(vm.count("drawLines"))
    {
      for(unsigned int i=1;i< vectVoxels.size(); i++)
      {
        viewer.addLine(vectVoxels.at(i-1), vectVoxels.at(i), lineSize);
      }
    }
    
    
    if(vm.count("drawVectors"))
    {
      std::string vectorsFileName = vm["drawVectors"].as<std::string>();
      std::vector<Z3i::RealPoint> vectorsPt = PointListReader<Z3i::RealPoint>::getPointsFromFile(vectorsFileName);
      if (vectorsPt.size()%2==1)
      {
        trace.info()<<"Warning the two set of points doesn't contains the same number of points, some vectors will be skipped." << std::endl;
      }
      int step=1;
      if(vm.count("filterVectors"))
        {
          double percentage = vm["filterVectors"].as<double>();
          step = max(1, (int) (100/percentage));
        }
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
    if(vm.count("addMesh"))
    {
      bool customColorMesh =  vm.count("customColorMesh");
      if(customColorMesh)
      {
        std::vector<unsigned int > vectCol = vm["customColorMesh"].as<std::vector<unsigned int> >();
        if(vectCol.size()!=4)
        {
          trace.error() << "colors specification should contain R,G,B and Alpha values"<< std::endl;
        }
        viewer.setFillColor(DGtal::Color(vectCol[0], vectCol[1], vectCol[2], vectCol[3]));
      }
      std::string meshName = vm["addMesh"].as<std::string>();
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

