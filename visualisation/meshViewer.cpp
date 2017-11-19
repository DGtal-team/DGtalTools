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
 * @file meshViewer.cpp
 * @ingroup visualisation
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2012/07/08
 *
 *
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <sstream>
#include "DGtal/base/Common.h"

#include "DGtal/io/Display3D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/PointListReader.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


using namespace std;
using namespace DGtal;


/**
 @page meshViewer meshViewer
 
 @brief Displays OFF mesh file by using QGLviewer.

 @b Usage:   meshViewer [input]

 @b Allowed @b options @b are :
 
 @code
  -h [ --help ]                    display this message
  -i [ --input ] arg               off files (.off), or OFS file (.ofs) 
  -x [ --scaleX ] arg (=1)         set the scale value in the X direction 
                                   (default 1.0)
  -y [ --scaleY ] arg (=1)         set the scale value in the Y direction 
                                   (default 1.0)
  -z [ --scaleZ ] arg (=1)         set the scale value in the Z direction 
                                   (default 1.0)
  -w [ --minLineWidth ] arg (=1.5) set the min line width of the mesh faces 
                                   (default 1.5)
  --customColorMesh arg            set the R, G, B, A components of the colors 
                                   of the mesh faces and eventually the color 
                                   R, G, B, A of the mesh edge lines (set by 
                                   default to black). 
  --customColorSDP arg             set the R, G, B, A components of the colors 
                                   of  the sdp view
  -f [ --displayVectorField ] arg  display a vector field from a simple sdp 
                                   file (two points per line)
  --vectorFieldIndex arg           specify special indices for the two point 
                                   coordinates (instead usinf the default 
                                   indices: 0 1, 2, 3, 4, 5)
  --customLineColor arg            set the R, G, B components of the colors of 
                                   the lines displayed from the 
                                   --displayVectorField option (red by 
                                   default). 
  -s [ --displaySDP ] arg          Add the display of a set of discrete points 
                                   as ball of radius 0.5.
  --SDPradius arg (=0.5)           change the ball radius to display a set of 
                                   discrete points (used with displaySDP 
                                   option)
  -n [ --invertNormal ]            threshold min to define binary shape
  -A [ --addAmbientLight ] arg(=0) add an ambient light for better display 
                                   (between 0 and 1).

  -v [ --drawVertex ]              draw the vertex of the mesh
  -d [ --doSnapShotAndExit]  arg,  save display snapshot into file. Notes that the camera setting is set by default according the last saved configuration (use SHIFT+Key_M to save current camera setting in the Viewer3D). If the camera setting was not saved it will use the default camera setting.

 @endcode


 @b Example: 


 @code
 $ meshViewer -i bunny.off
    
 @endcode

 You should obtain such a result:

 @image html resMeshViewer.png "Resulting visualization."
 

 @see
 @ref meshViewer.cpp

 */


/**
 * Custom Viewer3D to override KeyPressEvent method and handle new key display.
 * It also desactivate the double Rendering mode for more efficiency.
 **/
class CustomViewer3D: public Viewer3D<>
{
protected:

  virtual void init()
  {
    Viewer3D<>::init();
    Viewer3D<>::setKeyDescription ( Qt::Key_I, "Display mesh informations about #faces, #vertices" );
    Viewer3D<>::setGLDoubleRenderingMode(false);
    if(mySaveSnap){
      QObject::connect(this, SIGNAL(drawFinished(bool)), this, SLOT(saveSnapshot(bool)));
    }
  }
  virtual void keyPressEvent(QKeyEvent *e){
    bool handled = false;
    if( e->key() == Qt::Key_I)
    {
      handled=true;
      myIsDisplayingInfoMode = !myIsDisplayingInfoMode;
      stringstream ss;
      qglviewer::Vec camPos = camera()->position();
      DGtal::Z3i::RealPoint c (camPos[0], camPos[1], camPos[2]);
      ss << myInfoDisplay << " distance to camera: " << (c-centerMesh).norm();
      Viewer3D<>::displayMessage(QString(myIsDisplayingInfoMode ?
                                         ss.str().c_str() : " "), 1000000);
      
      Viewer3D<>::update();
    }
    
    if(!handled)
      {
        Viewer3D<>::keyPressEvent(e);        
      }
  };

public: 
  std::string myInfoDisplay = "No information loaded...";
  bool myIsDisplayingInfoMode = false;
  bool mySaveSnap = false;
  DGtal::Z3i::RealPoint centerMesh;
};


  

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
  ("help,h", "display this message")
  ("input,i", po::value<std::vector<string> >()->multitoken(), "off files (.off), or OFS file (.ofs) " )
  ("scaleX,x",  po::value<float>()->default_value(1.0), "set the scale value in the X direction (default 1.0)" )
  ("scaleY,y",  po::value<float>()->default_value(1.0), "set the scale value in the Y direction (default 1.0)" )
  ("scaleZ,z",  po:: value<float>()->default_value(1.0), "set the scale value in the Z direction (default 1.0)")
  ("minLineWidth,w",  po:: value<float>()->default_value(1.5), "set the min line width of the mesh faces (default 1.5)")
  ("customColorMesh",po::value<std::vector<unsigned int> >()->multitoken(), "set the R, G, B, A components of the colors of the mesh faces and eventually the color R, G, B, A of the mesh edge lines (set by default to black). " )
  ("customColorSDP",po::value<std::vector<unsigned int> >()->multitoken(), "set the R, G, B, A components of the colors of  the sdp view" )
  ("displayVectorField,f",po::value<std::string>(),  "display a vector field from a simple sdp file (two points per line)" )
  ("vectorFieldIndex",po::value<std::vector<unsigned int> >()->multitoken(), "specify special indices for the two point coordinates (instead usinf the default indices: 0 1, 2, 3, 4, 5)" )
  ("customLineColor",po::value<std::vector<unsigned int> >()->multitoken(), "set the R, G, B components of the colors of the lines displayed from the --displayVectorField option (red by default). " )
  ("displaySDP,s", po::value<std::string>(), "add the display of a set of discrete points as ball of radius 0.5.")
  ("SDPradius", po::value<double>()->default_value(0.5), "change the ball radius to display a set of discrete points (used with displaySDP option)")
  ("invertNormal,n", "invert face normal vectors." )
  ("addAmbientLight,A",po:: value<float>()->default_value(0.0), "add an ambient light for better display (between 0 and 1)." )
  ("drawVertex,v", "draw the vertex of the mesh" )
  ("doSnapShotAndExit,d", po::value<std::string>(), "save display snapshot into file. Notes that the camera setting is set by default according the last saved configuration (use SHIFT+Key_M to save current camera setting in the Viewer3D). If the camera setting was not saved it will use the default camera setting." );
  
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
    std::cout << "Usage: " << argv[0] << " [input]\n"
    << "Display OFF mesh file by using QGLviewer"
    << general_opt << "\n";
    return 0;
  }
  
  if(! vm.count("input"))
  {
    trace.error() << " The file name was defined" << endl;
    return 0;
  }
  
  
  
  std::vector<std::string> inputFilenameVect = vm["input"].as<std::vector<std::string > >();
  float sx = vm["scaleX"].as<float>();
  float sy = vm["scaleY"].as<float>();
  float sz = vm["scaleZ"].as<float>();
  
  unsigned int  meshColorR = 240;
  unsigned int  meshColorG = 240;
  unsigned int  meshColorB = 240;
  unsigned int  meshColorA = 255;
  
  unsigned int  meshColorRLine = 0;
  unsigned int  meshColorGLine = 0;
  unsigned int  meshColorBLine = 0;
  unsigned int  meshColorALine = 255;
  
  
  unsigned int  sdpColorR = 240;
  unsigned int  sdpColorG = 240;
  unsigned int  sdpColorB = 240;
  unsigned int  sdpColorA = 255;
  
  
  bool displayVectorField = vm.count("displayVectorField");
  std::vector<unsigned int> vectFieldIndices = {0,1,2,3,4,5};
  
  if (displayVectorField) {
    if(vm.count("vectorFieldIndex")){
      vectFieldIndices = vm["vectorFieldIndex"].as<std::vector<unsigned int> >();
      if (vectFieldIndices.size() != 6) {
        trace.warning() << "you should specify indices for each of the 6 fields of the two coordinates." << std::endl;
        vectFieldIndices = {0,1,2,3,4,5};
      }
    }
  }
  
  float lineWidth = vm["minLineWidth"].as<float>();
  
  DGtal::Color vFieldLineColor = DGtal::Color::Red;
  if(vm.count("customLineColor")){
    std::vector<unsigned int > vectCol = vm["customLineColor"].as<std::vector<unsigned int> >();
    if(vectCol.size()!=3 ){
      trace.error() << "colors specification should contain R,G,B values (using default red)."<< std::endl;
    }
    vFieldLineColor.setRGBi(vectCol[0], vectCol[1], vectCol[2], 255);
  }
  
  if(vm.count("customColorMesh")){
    std::vector<unsigned int > vectCol = vm["customColorMesh"].as<std::vector<unsigned int> >();
    if(vectCol.size()!=4 && vectCol.size()!=8 ){
      trace.error() << "colors specification should contain R,G,B and Alpha values"<< std::endl;
    }
    meshColorR = vectCol[0];
    meshColorG = vectCol[1];
    meshColorB = vectCol[2];
    meshColorA = vectCol[3];
    if(vectCol.size() == 8){
      meshColorRLine = vectCol[4];
      meshColorGLine = vectCol[5];
      meshColorBLine = vectCol[6];
      meshColorALine = vectCol[7];
      
    }
    
  }
  if(vm.count("customColorSDP")){
    std::vector<unsigned int > vectCol = vm["customColorSDP"].as<std::vector<unsigned int> >();
    if(vectCol.size()!=4){
      trace.error() << "colors specification should contain R,G,B and Alpha values"<< std::endl;
    }
    sdpColorR = vectCol[0];
    sdpColorG = vectCol[1];
    sdpColorB = vectCol[2];
    sdpColorA = vectCol[3];
  }
  
  
  
  QApplication application(argc,argv);
  CustomViewer3D viewer;
  viewer.mySaveSnap = vm.count("doSnapShotAndExit");
  if(vm.count("doSnapShotAndExit")){
    viewer.setSnapshotFileName(QString(vm["doSnapShotAndExit"].as<std::string>().c_str()));
  }

  std::stringstream title;
  title  << "Simple Mesh Viewer: " << inputFilenameVect[0];
  viewer.setWindowTitle(title.str().c_str());
  viewer.show();
  viewer.myGLLineMinWidth = lineWidth;
  viewer.setGLScale(sx, sy, sz);
  bool invertNormal= vm.count("invertNormal");  
  double ballRadius = vm["SDPradius"].as<double>();
  if(vm.count("addAmbientLight"))
    {
      float val = vm["addAmbientLight"].as<float>();
      GLfloat lightAmbientCoeffs [4] = {val,val, val, 1.0f};
      viewer.setGLLightAmbientCoefficients(lightAmbientCoeffs);
    }
  
  trace.info() << "Importing mesh... ";
  
  std::vector<Mesh<DGtal::Z3i::RealPoint> >  vectMesh;
  for(unsigned int i = 0; i< inputFilenameVect.size(); i++){
    Mesh<DGtal::Z3i::RealPoint> aMesh(!vm.count("customColorMesh"));
    aMesh << inputFilenameVect[i];
    vectMesh.push_back(aMesh);
  }
  DGtal::Z3i::RealPoint centerMeshes;
  unsigned int tot=0;
  for(const auto & m: vectMesh)
    {
      for( auto p = m.vertexBegin(); p!=m.vertexEnd(); ++p)
        centerMeshes += *p;
      tot+=m.nbVertex();
    }
  centerMeshes /= tot;
  viewer.centerMesh = centerMeshes;
  bool import = vectMesh.size()==inputFilenameVect.size();
  if(!import){
    trace.info() << "File import failed. " << std::endl;
    return 0;
  }
  
  trace.info() << "[done]. "<< std::endl;
  if(vm.count("displaySDP")){
    std::string filenameSDP = vm["displaySDP"].as<std::string>();
    vector<Z3i::RealPoint> vectPoints;
    vectPoints = PointListReader<Z3i::RealPoint>::getPointsFromFile(filenameSDP);
    viewer << CustomColors3D(Color(sdpColorR, sdpColorG, sdpColorB, sdpColorA),
                             Color(sdpColorR, sdpColorG, sdpColorB, sdpColorA));
    for(unsigned int i=0;i< vectPoints.size(); i++){
      viewer.addBall(vectPoints.at(i), ballRadius);
    }
  }
  if(invertNormal){
    for(unsigned int i=0; i<vectMesh.size(); i++){
      vectMesh[i].invertVertexFaceOrder();
    }
  }
  
  viewer << CustomColors3D(Color(meshColorRLine, meshColorGLine, meshColorBLine, meshColorALine),
                           Color(meshColorR, meshColorG, meshColorB, meshColorA));
  for(unsigned int i=0; i<vectMesh.size(); i++){
    viewer << vectMesh[i];
  }
  
  if(vm.count("drawVertex")){
    for(unsigned int i=0; i<vectMesh.size(); i++){
      for( Mesh<DGtal::Z3i::RealPoint>::VertexStorage::const_iterator it = vectMesh[i].vertexBegin();
          it!=vectMesh[i].vertexEnd(); ++it){
        DGtal::Z3i::Point pt;
        pt[0]=(*it)[0]; pt[1]=(*it)[1]; pt[2]=(*it)[2];
        viewer << pt;
      }
    }
  }

  
  if (displayVectorField) {
    std::vector<unsigned int > vectFieldIndices1 = {vectFieldIndices[0],vectFieldIndices[1], vectFieldIndices[2]};
    std::vector<unsigned int > vectFieldIndices2 = {vectFieldIndices[3],vectFieldIndices[4], vectFieldIndices[5]};
    
    std::vector<DGtal::Z3i::RealPoint> vectPt1 = PointListReader<DGtal::Z3i::RealPoint>::getPointsFromFile(vm["displayVectorField"].as<std::string>(), vectFieldIndices1);
    std::vector<DGtal::Z3i::RealPoint> vectPt2 = PointListReader<DGtal::Z3i::RealPoint>::getPointsFromFile(vm["displayVectorField"].as<std::string>(), vectFieldIndices2);
    viewer.createNewLineList();
    for (unsigned int i = 0; i < vectPt1.size(); i++) {
    
      viewer.setLineColor(vFieldLineColor);
      viewer.addLine(vectPt1[i], vectPt2[i]);
    }
  }
  unsigned int nbVertex = 0;
  unsigned int nbFaces = 0;
  for(auto const &m:  vectMesh)
    {
      nbVertex += m.nbVertex();
      nbFaces +=m.nbFaces();
    }
  stringstream ss;
  ss << "# faces: " << std::fixed << nbFaces << "    #vertex: " <<  nbVertex ;
  viewer.myInfoDisplay = ss.str();
  viewer  << CustomViewer3D::updateDisplay;
  if(vm.count("doSnapShotAndExit")){
    // Appy cleaning just save the last snap
    if(!viewer.restoreStateFromFile())
      {
        viewer.update();
      }    
    std::string name = vm["doSnapShotAndExit"].as<std::string>();
    std::string extension = name.substr(name.find_last_of(".") + 1);
    std::string basename = name.substr(0, name.find_last_of("."));
    for(int i=0; i< viewer.snapshotCounter()-1; i++){
      std::stringstream s;
      s << basename << "-"<< setfill('0') << setw(4)<<  i << "." << extension;
      trace.info() << "erase temp file: " << s.str() << std::endl;
      remove(s.str().c_str());
    }

    std::stringstream s;
    s << basename << "-"<< setfill('0') << setw(4)<<  viewer.snapshotCounter()-1 << "." << extension;
    rename(s.str().c_str(), name.c_str());
    return 0;
  }

  return application.exec();
}
