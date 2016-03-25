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
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/shapes/Mesh.h"
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
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
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "input file: sdp (sequence of discrete points)" )
    ("SDPindex", po::value<std::vector <unsigned int> >()->multitoken(), "specify the sdp index (by default 0,1,2).")
    ("pointColor,c", po::value<std::vector <int> >()->multitoken(), "set the color of  points: r g b a " )
    ("addMesh,m", po::value<std::string>(), "append a mesh (off/obj) to the point set visualization.")
    ("lineColor,l",po::value<std::vector <int> >()->multitoken(), "set the color of line: r g b a " )
    ("colorFromLabels", "use the color indexed from labels in the file.")
    ("labelsIndex", po::value<unsigned int>(), "define the index of the label in the source file (used by --LabelsIndex) ")
    ("filter,f",po::value<double>()->default_value(100.0), "filter input file in order to display only the [arg] pourcent of the input 3D points (uniformly selected)." )
    ("noPointDisplay", "usefull for instance to only display the lines between points.")
    ("drawLines", "draw the line between discrete points." )
    ("scaleX,x",  po::value<float>()->default_value(1.0), "set the scale value in the X direction (default 1.0)" )
    ("scaleY,y",  po::value<float>()->default_value(1.0), "set the scale value in the Y direction (default 1.0)" )
    ("scaleZ,z",  po::value<float>()->default_value(1.0), "set the scale value in the Z direction (default 1.0)")
    ("sphereRadius,s",  po::value<double>()->default_value(0.2), "defines the sphere radius (used when the primitive is set to the sphere). (default value 0.2)")
    ("sphereResolution",  po::value<unsigned int>()->default_value(30), "defines the sphere resolution (used when the primitive is set to the sphere). (default resolution: 30)")
    ("lineSize",  po::value<double>()->default_value(0.2), "defines the line size (used when the --drawLines or --drawVectors option is selected). (default value 0.2))")
    ("primitive,p", po::value<std::string>()->default_value("voxel"), "set the primitive to display the set of points (can be sphere or voxel (default)")
    ("drawVectors,v", po::value<std::string>(), "SDP vector file: draw a set of vectors from the given file (each vector are determined by two consecutive point given, each point represented by its coordinates on a single line.")
    ("interactiveDisplayVoxCoords", "by using this option the pixel coordinates can be displayed after selection (shift+left click on voxel)." );


  bool parseOK=true;
  bool cannotStart= false;

  typedef PointVector<1, int> Point1D;


  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.error()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);
  if(parseOK && ! vm.count("input"))
    {
      trace.error() << " The input file name was not defined" << endl;
      cannotStart = true;
    }
  std::string typePrimitive;
  double sphereRadius = 0.2;
  unsigned int sphereResolution = vm["sphereResolution"].as<unsigned int>();
  double lineSize =0.2;
  
  Color lineColor(100, 100, 250);
  Color pointColor(250, 250, 250);
  if(parseOK){
    typePrimitive = vm["primitive"].as<std::string>();
    sphereRadius = vm["sphereRadius"].as<double>();
    lineSize = vm["lineSize"].as<double>();
  }

  if (parseOK && typePrimitive !="voxel" && typePrimitive  != "sphere" ){
    trace.error() << " The primitive should be sphere or voxel (primitive: "
                  << typePrimitive << " not implemented)" << std::endl;
    cannotStart = true;
  }

  if(parseOK && vm.count("lineColor")){
    std::vector<int> vcol= vm["lineColor"].as<std::vector<int > >();
    if(vcol.size()<4){
      trace.error() << " Not enough parameter: color specification should contains four elements: red, green, blue and alpha values "
                    << "(Option --lineColor ignored). "  << std::endl;
    }
    lineColor.setRGBi(vcol[0], vcol[1], vcol[2], vcol[3]);
  }
  if(parseOK && vm.count("pointColor")){
    std::vector<int> vcol= vm["pointColor"].as<std::vector<int > >();
    if(vcol.size()<4){
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

  bool colorFromLabels = vm.count("colorFromLabels");
  bool interactiveDisplayVoxCoords = vm.count("interactiveDisplayVoxCoords");

  typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;
  Z3i::KSpace K;
  Viewer viewer( K );
  viewer.setWindowTitle("3dSPD Viewer");
  viewer.show();
  viewer.setGLScale(sx, sy, sz);
  viewer.myGLLineMinWidth = lineSize;
  viewer << CustomColors3D(pointColor, pointColor);


  // Get vector of labels if exists.
  std::vector<unsigned int> vectLabels;
  if(colorFromLabels){
    std::vector<Point1D> vectVal;
    std::vector<unsigned int> vectIndex;
    if(vm.count("labelsIndex")){
      vectIndex.push_back(vm["labelsIndex"].as<unsigned int>());
    }else{
      vectIndex.push_back(3);
    }
    vectVal = PointListReader<Point1D>::getPointsFromFile(inputFilename, vectIndex);
    for(std::vector<Point1D>::iterator it= vectVal.begin(); it!=vectVal.end(); it++){
      vectLabels.push_back((unsigned int)(*it)[0]);
    }
  }


  GradientColorMap< int > gradientColorMap( 1, (!colorFromLabels)? 1:  * std::max_element(vectLabels.begin(), vectLabels.end()));
  if(colorFromLabels){
    gradientColorMap.addColor( Color(255,100,100 ) );
    gradientColorMap.addColor( Color(1,100, 255 ) );
  }
  vector<Z3i::RealPoint> vectVoxels;
  if(vm.count("SDPindex")) {
    std::vector<unsigned int > vectIndex = vm["SDPindex"].as<std::vector<unsigned int > >();
    if(vectIndex.size()!=3){
      trace.error() << "you need to specify the three indexes of vertex." << std::endl;
      return 0;
    }
    vectVoxels = PointListReader<Z3i::RealPoint>::getPointsFromFile(inputFilename, vectIndex);
  }else{
    vectVoxels = PointListReader<Z3i::RealPoint>::getPointsFromFile(inputFilename);
  }
  int name = 0;
  if(!vm.count("noPointDisplay")){
    double percent = vm["filter"].as<double>();
    int step = max(1, (int) (100/percent));
    for(unsigned int i=0;i< vectVoxels.size(); i=i+step){
      if(colorFromLabels){
        Color col = gradientColorMap((int) vectLabels.at(i));
        viewer.setFillColor(col);
      }

      if(typePrimitive=="voxel"){
        if (interactiveDisplayVoxCoords)
          {
            viewer << SetName3D( name++ ) ;
          }
        viewer << Z3i::Point((int)vectVoxels.at(i)[0],
                             (int)vectVoxels.at(i)[1],
                             (int)vectVoxels.at(i)[2]);
      }else{
        viewer.addBall(vectVoxels.at(i), sphereRadius, sphereResolution);
      }
    }
  }

  viewer << CustomColors3D(lineColor, lineColor);
  if(vm.count("drawLines")){
    for(unsigned int i=1;i< vectVoxels.size(); i++){
      viewer.addLine(vectVoxels.at(i-1), vectVoxels.at(i), lineSize);
    }
  }


  if(vm.count("drawVectors")){
    std::string vectorsFileName = vm["drawVectors"].as<std::string>();
    std::vector<Z3i::RealPoint> vectorsPt = PointListReader<Z3i::RealPoint>::getPointsFromFile(vectorsFileName);
    if (vectorsPt.size()%2==1){
      trace.info()<<"Warning the two set of points doesn't contains the same number of points, some vectors will be skipped." << std::endl;
    }
    for(unsigned int i =0; i<vectorsPt.size()-1; i=i+2){
      viewer.addLine(vectorsPt.at(i),vectorsPt.at(i+1), lineSize);
    }

  }
  if(vm.count("addMesh")){
    viewer.setFillColor(DGtal::Color::White);
    std::string meshName = vm["addMesh"].as<std::string>();
    Mesh<Z3i::RealPoint> mesh;
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

