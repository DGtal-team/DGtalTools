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

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;

/**
 @page meshViewer meshViewer

 @brief Displays OFF mesh file by using QGLviewer.

 @b Usage:   meshViewer [input]

 @b Allowed @b options @b are :

 @code

 Positionals:
 1 TEXT:FILE ... REQUIRED              inputFileNames.off files (.off), or OFS file (.ofs)

 Options:
 -h,--help                             Print this help message and exit
 -i,--input TEXT:FILE ... REQUIRED     inputFileNames.off files (.off), or OFS file (.ofs)
 -x,--scaleX FLOAT                     set the scale value in the X direction (default 1.0)
 -y,--scaleY FLOAT                     set the scale value in the y direction (default 1.0)
 -z,--scaleZ FLOAT                     set the scale value in the z direction (default 1.0)
 --minLineWidth FLOAT=1.5              set the min line width of the mesh faces (default 1.5)
 --customColorMesh UINT ...            set the R, G, B, A components of the colors of the mesh faces and eventually the color R, G, B, A of the mesh edge lines (set by default to black).
 --customAlphaMesh UINT ...            set the alpha components of the colors of the mesh faces (can be applied for each mesh).
 --customColorSDP UINT x 4             set the R, G, B, A components of the colors of  the sdp view
 -f,--displayVectorField TEXT          display a vector field from a simple sdp file (two points per line)
 --vectorFieldIndex UINT x 6           specify special indices for the two point coordinates (instead usinf the default indices: 0 1, 2, 3, 4, 5)
 --customLineColor UINT x 4            set the R, G, B components of the colors of the lines displayed from the --displayVectorField option (red by default).
 --SDPradius FLOAT=0.5                 change the ball radius to display a set of discrete points (used with displaySDP option)
 -s,--displaySDP TEXT                  add the display of a set of discrete points as ball of radius 0.5.
 -A,--addAmbientLight FLOAT            add an ambient light for better display (between 0 and 1).
 -b,--customBGColor UINT x 3           set the R, G, B components of the colors of the background color.
 -d,--doSnapShotAndExit TEXT           save display snapshot into file. Notes that the camera                                             setting is set by default according the last saved                                                 configuration (use SHIFT+Key_M to save current camera                                              setting in the Viewer3D). If the camera setting was not                                            saved it will use the default camera setting.
 -c,--useLastCameraSetting             use the last camera setting of the user (i.e if a                                                  .qglviewer.xml file is present in the current directory)
 -l,--fixLightToScene                  Fix light source to scence instead to camera
 -n,--invertNormal                     invert face normal vectors.
 -v,--drawVertex                       draw the vertex of the mesh

 @endcode


 @b Example:

 @code
 $ meshViewer bunny.off

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
class CustomViewer3D : public Viewer3D<>
{
protected:
  virtual void init()
  {
    Viewer3D<>::init();
    Viewer3D<>::setKeyDescription(Qt::Key_I, "Display mesh informations about #faces, #vertices");
    Viewer3D<>::setGLDoubleRenderingMode(false);
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    bool handled = false;
    if (e->key() == Qt::Key_I)
    {
      handled = true;
      myIsDisplayingInfoMode = !myIsDisplayingInfoMode;
      stringstream ss;
      qglviewer::Vec camPos = camera()->position();
      DGtal::Z3i::RealPoint c(camPos[0], camPos[1], camPos[2]);
      ss << myInfoDisplay << " distance to camera: " << (c - centerMesh).norm();
      Viewer3D<>::displayMessage(QString(myIsDisplayingInfoMode ? ss.str().c_str() : " "), 1000000);

      Viewer3D<>::update();
    }

    if (!handled)
    {
      Viewer3D<>::keyPressEvent(e);
    }
  };

public:
  void changeDefaultBGColor(const DGtal::Color &col)
  {
    myDefaultBackgroundColor = col;
    Viewer3D<>::update();
    Viewer3D<>::draw();
  }
  std::string myInfoDisplay = "No information loaded...";
  bool myIsDisplayingInfoMode = false;
  DGtal::Z3i::RealPoint centerMesh;
};

int main(int argc, char **argv)
{
  float sx{1.0};
  float sy{1.0};
  float sz{1.0};

  unsigned int meshColorR{240};
  unsigned int meshColorG{240};
  unsigned int meshColorB{240};
  unsigned int meshColorA{255};

  unsigned int meshColorRLine{0};
  unsigned int meshColorGLine{0};
  unsigned int meshColorBLine{0};
  unsigned int meshColorALine{255};

  unsigned int sdpColorR{240};
  unsigned int sdpColorG{240};
  unsigned int sdpColorB{240};
  unsigned int sdpColorA{255};

  float lineWidth{1.5};

  std::vector<unsigned int> customColorMesh;
  std::vector<unsigned int> customColorSDP;
  std::vector<unsigned int> customLineColor;
  std::vector<unsigned int> customBGColor;
  std::vector<unsigned int> customAlphaMesh;
  std::vector<unsigned int> vectFieldIndices = {0, 1, 2, 3, 4, 5};
  std::string displayVectorField;

  std::string snapshotFile;
  std::string filenameSDP;
  double ballRadius{0.5};
  bool invertNormal{false};
  bool drawVertex{false};
  bool useLastCamSet{false};
  bool fixLightToScene{false};
  float ambiantLight{0.0};

  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::vector<std::string> inputFileNames;
  std::string outputFileName{"result.raw"};
  app.description("Display OFF mesh file by using QGLviewer");
  app.add_option("-i,--input,1", inputFileNames, "inputFileNames.off files (.off), or OFS file (.ofs)")
      ->check(CLI::ExistingFile)
      ->required();
  app.add_option("-x,--scaleX", sx, "set the scale value in the X direction (default 1.0)");
  app.add_option("-y,--scaleY", sy, "set the scale value in the y direction (default 1.0)");
  app.add_option("-z,--scaleZ", sz, "set the scale value in the z direction (default 1.0)");
  app.add_option("--minLineWidth", lineWidth, "set the min line width of the mesh faces (default 1.5)", true);
  app.add_option("--customColorMesh", customColorMesh, "set the R, G, B, A components of the colors of the mesh faces and eventually the color R, G, B, A of the mesh edge lines (set by default to black).");
  app.add_option("--customAlphaMesh", customAlphaMesh, "set the alpha(A) components of the colors of the mesh faces (can be applied for each mesh).");
  app.add_option("--customColorSDP", customColorSDP, "set the R, G, B, A components of the colors of the sdp view")
      ->expected(4);
  app.add_option("--displayVectorField,-f", displayVectorField, "display a vector field from a simple sdp file (two points per line)");
  app.add_option("--vectorFieldIndex", vectFieldIndices, "specify special indices for the two point coordinates (instead usinf the default indices: 0 1, 2, 3, 4, 5)")
      ->expected(6);
  app.add_option("--customLineColor", customLineColor, "set the R, G, B components of the colors of the lines displayed from the --displayVectorField option (red by default).")
      ->expected(4);
  app.add_option("--SDPradius", ballRadius, "change the ball radius to display a set of discrete points (used with displaySDP option)", true);
  app.add_option("--displaySDP,-s", filenameSDP, "add the display of a set of discrete points as ball of radius 0.5.");
  app.add_option("--addAmbientLight,-A", ambiantLight, "add an ambient light for better display (between 0 and 1).");
  app.add_option("--customBGColor,-b", customBGColor, "set the R, G, B components of the colors of the background color.")
      ->expected(3);
  app.add_option("--doSnapShotAndExit,-d", snapshotFile, "save display snapshot into file. Notes that the camera setting is set by default according the last saved configuration (use SHIFT+Key_M to save current camera setting in the Viewer3D). If the camera setting was not saved it will use the default camera setting.");
  app.add_flag("--useLastCameraSetting,-c", useLastCamSet, "use the last camera setting of the user (i.e if a .qglviewer.xml file is present in the current directory)");
  app.add_flag("--fixLightToScene,-l", fixLightToScene, "Fix light source to scence instead to camera");
  app.add_flag("--invertNormal,-n", invertNormal, "invert face normal vectors.");
  app.add_flag("--drawVertex,-v", drawVertex, "draw the vertex of the mesh");

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  DGtal::Color vFieldLineColor = DGtal::Color::Red;
  if (customLineColor.size() == 4)
  {
    vFieldLineColor.setRGBi(customLineColor[0], customLineColor[1], customLineColor[2], 255);
  }

  if (customColorMesh.size() != 0)
  {
    if (customColorMesh.size() < 4)
    {
      trace.error() << "colors specification should contain R,G,B and Alpha values" << std::endl;
    }
  }

  if (customColorSDP.size() == 4)
  {
    sdpColorR = customColorSDP[0];
    sdpColorG = customColorSDP[1];
    sdpColorB = customColorSDP[2];
    sdpColorA = customColorSDP[3];
  }

  QApplication application(argc, argv);
  CustomViewer3D viewer;
  if (snapshotFile != "")
  {
    viewer.setSnapshotFileName(QString(snapshotFile.c_str()));
  }

  std::stringstream title;
  title << "Simple Mesh Viewer: " << inputFileNames[0];
  viewer.setWindowTitle(title.str().c_str());
  viewer.show();
  viewer.myGLLineMinWidth = lineWidth;
  viewer.setGLScale(sx, sy, sz);
  if (customBGColor.size() == 3)
  {
    viewer.changeDefaultBGColor(DGtal::Color(customBGColor[0],
                                             customBGColor[1],
                                             customBGColor[2], 255));
  }
  if (ambiantLight != 0.0)
  {
    GLfloat lightAmbientCoeffs[4] = {ambiantLight, ambiantLight, ambiantLight, 1.0f};
    viewer.setGLLightAmbientCoefficients(lightAmbientCoeffs);
  }

  trace.info() << "Importing mesh... ";

  std::vector<Mesh<DGtal::Z3i::RealPoint>> vectMesh;
  for (unsigned int i = 0; i < inputFileNames.size(); i++)
  {
    Mesh<DGtal::Z3i::RealPoint> aMesh(customColorMesh.size() != 4 && customColorMesh.size() != 8);
    aMesh << inputFileNames[i];
    // for obj mesh by default the mesh color face are not necessary uniform.
    if (aMesh.isStoringFaceColors() && customColorMesh.size() >= 4)
    {
      if (i * 8 < customColorMesh.size())
      {
        meshColorR = customColorMesh[i * 8];
      }
      if (i * 8 + 1 < customColorMesh.size())
      {
        meshColorG = customColorMesh[i * 8 + 1];
      }
      if (i * 8 + 2 < customColorMesh.size())
      {
        meshColorB = customColorMesh[i * 8 + 2];
      }
      if (i * 8 + 3 < customColorMesh.size())
      {
        meshColorA = customColorMesh[i * 8 + 3];
      }
      for (unsigned int j = 0; j < aMesh.nbFaces(); j++)
      {
        aMesh.setFaceColor(j, Color(meshColorR, meshColorG, meshColorB, meshColorA));
      }
    }
    else if (customAlphaMesh.size() > 0)
    {
      for (unsigned int j = 0; j < aMesh.nbFaces(); j++)
      {
        auto c = aMesh.getFaceColor(j);
        aMesh.setFaceColor(j, Color(c.red(), c.green(), c.blue(), customAlphaMesh.at(i < customAlphaMesh.size() ? i : 0)));
      }
    }

    vectMesh.push_back(aMesh);
  }
  DGtal::Z3i::RealPoint centerMeshes;
  unsigned int tot = 0;
  for (const auto &m : vectMesh)
  {
    for (auto p = m.vertexBegin(); p != m.vertexEnd(); ++p)
      centerMeshes += *p;
    tot += m.nbVertex();
  }
  centerMeshes /= tot;
  viewer.centerMesh = centerMeshes;
  bool import = vectMesh.size() == inputFileNames.size();
  if (!import)
  {
    trace.info() << "File import failed. " << std::endl;
    return 0;
  }

  trace.info() << "[done]. " << std::endl;
  if (filenameSDP != "")
  {
    if (customAlphaMesh.size() > 0)
    {
      trace.info() << "New meshViewer" << std::endl;
      auto vOrigins = PointListReader<PointVector<1, DGtal::int32_t>>::getPolygonsFromFile(filenameSDP);
      viewer << CustomColors3D(Color(sdpColorR, sdpColorG, sdpColorB, sdpColorA),
                               Color(sdpColorR, sdpColorG, sdpColorB, sdpColorA));
      for (auto l : vOrigins)
      {
        DGtal::Z3i::Point pt(l[0][0], l[1][0], l[2][0]);
        DGtal::Color cl(l[3][0], l[4][0], l[5][0], sdpColorA);
        viewer.setFillColor(cl);
        viewer.addBall(pt, ballRadius);
      }
    }
    else
    {
      vector<Z3i::RealPoint> vectPoints;
      vectPoints = PointListReader<Z3i::RealPoint>::getPointsFromFile(filenameSDP);
      viewer << CustomColors3D(Color(sdpColorR, sdpColorG, sdpColorB, sdpColorA),
                               Color(sdpColorR, sdpColorG, sdpColorB, sdpColorA));
      for (unsigned int i = 0; i < vectPoints.size(); i++)
      {
        viewer.addBall(vectPoints.at(i), ballRadius);
      }
    }
  }
  if (invertNormal)
  {
    for (unsigned int i = 0; i < vectMesh.size(); i++)
    {
      vectMesh[i].invertVertexFaceOrder();
    }
  }

  for (unsigned int i = 0; i < vectMesh.size(); i++)
  {
    if (i * 8 < customColorMesh.size())
    {
      meshColorR = customColorMesh[i * 8];
    }
    if (i * 8 + 1 < customColorMesh.size())
    {
      meshColorG = customColorMesh[i * 8 + 1];
    }
    if (i * 8 + 2 < customColorMesh.size())
    {
      meshColorB = customColorMesh[i * 8 + 2];
    }
    if (i * 8 + 3 < customColorMesh.size())
    {
      meshColorA = customColorMesh[i * 8 + 3];
    }
    if (i * 8 + 4 < customColorMesh.size())
    {
      meshColorALine = customColorMesh[i * 8 + 4];
    }
    if (i * 8 + 5 < customColorMesh.size())
    {
      meshColorBLine = customColorMesh[i * 8 + 5];
    }
    if (i * 8 + 6 < customColorMesh.size())
    {
      meshColorRLine = customColorMesh[i * 8 + 6];
    }
    if (i * 8 + 7 < customColorMesh.size())
    {
      meshColorALine = customColorMesh[i * 8 + 7];
    }

    viewer << CustomColors3D(Color(meshColorRLine, meshColorGLine, meshColorBLine,
                                   meshColorALine),
                             Color(meshColorR, meshColorG, meshColorB, meshColorA));

    viewer << vectMesh[i];
  }

  if (drawVertex)
  {
    for (unsigned int i = 0; i < vectMesh.size(); i++)
    {
      for (Mesh<DGtal::Z3i::RealPoint>::VertexStorage::const_iterator it = vectMesh[i].vertexBegin();
           it != vectMesh[i].vertexEnd(); ++it)
      {
        DGtal::Z3i::Point pt;
        pt[0] = (*it)[0];
        pt[1] = (*it)[1];
        pt[2] = (*it)[2];
        viewer << pt;
      }
    }
  }

  if (displayVectorField != "")
  {
    std::vector<unsigned int> vectFieldIndices1 = {vectFieldIndices[0], vectFieldIndices[1], vectFieldIndices[2]};
    std::vector<unsigned int> vectFieldIndices2 = {vectFieldIndices[3], vectFieldIndices[4], vectFieldIndices[5]};
    std::vector<DGtal::Z3i::RealPoint> vectPt1 = PointListReader<DGtal::Z3i::RealPoint>::getPointsFromFile(displayVectorField, vectFieldIndices1);
    std::vector<DGtal::Z3i::RealPoint> vectPt2 = PointListReader<DGtal::Z3i::RealPoint>::getPointsFromFile(displayVectorField, vectFieldIndices2);
    viewer.createNewLineList();
    for (unsigned int i = 0; i < vectPt1.size(); i++)
    {
      viewer.setLineColor(vFieldLineColor);
      viewer.addLine(vectPt1[i], vectPt2[i]);
    }
  }
  unsigned int nbVertex = 0;
  unsigned int nbFaces = 0;
  for (auto const &m : vectMesh)
  {
    nbVertex += m.nbVertex();
    nbFaces += m.nbFaces();
  }
  stringstream ss;
  ss << "# faces: " << std::fixed << nbFaces << "    #vertex: " << nbVertex;
  viewer.myInfoDisplay = ss.str();
  viewer << CustomViewer3D::updateDisplay;
  if (fixLightToScene)
  {
    viewer.setLightModeFixToCamera(false, false);
  }

  if (useLastCamSet)
  {
    viewer.restoreStateFromFile();
  }
  else
  {
    // useful in non interactive case in order to retain the default camera settings (that are not saved in case of process kill).
    viewer.saveStateToFile();
  }
  // First display transparency improve
  viewer.sortTriangleFromCamera();
  viewer.sortQuadFromCamera();
  viewer.sortSurfelFromCamera();
  viewer.sortPolygonFromCamera();
  viewer << CustomViewer3D::updateDisplay;

  if (snapshotFile != "")
  {
    // Recover mesh position
    viewer.restoreStateFromFile();
    viewer.saveSnapshot(QString(snapshotFile.c_str()), true);
    return 0;
  }
  trace.info() << "[display ready]" << std::endl;

  return application.exec();
}
