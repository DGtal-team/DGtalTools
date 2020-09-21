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
 * @file 3dImageViewer.cpp
 * @ingroup visualisation
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2011/01/04
 *
 * An example file named qglViewer.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/SurfelAdjacency.h"

#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#ifdef WITH_ITK
#include "DGtal/io/readers/DicomReader.h"
#endif

#include "DGtal/images/ImageSelector.h"



#include "specificClasses/Viewer3DImage.cpp"

#include "CLI11.hpp"



using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
   @page Doc3dImageViewer 3dImageViewer
 
   @brief Displays volume file as a voxel set by using QGLviewer.

   @b Usage:  3dImageViewer [OPTIONS] 1 [s]

   @b Allowed @b options @b are :
 
   @code

   Positionals:
     1 TEXT:FILE REQUIRED                  vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
     s TEXT                                display a set of discrete points (.sdp)

   Options:
     -h,--help                             Print this help message and exit
     -i,--input TEXT:FILE REQUIRED         vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
     --grid                                draw slice images using grid mode.
     --intergrid                           draw slice images using inter grid mode.
     --emptyMode                           remove the default boundingbox display.
     --thresholdImage                      threshold the image to define binary shape
     --thresholdMin INT=0                  threshold min to define binary shape
     --thresholdMax INT=255                threshold maw to define binary shape
     --displaySDP TEXT                     display a set of discrete points (.sdp)
     --SDPindex UINT x 3                   specify the sdp index.
     --SDPball FLOAT=0                     use balls to display a set of discrete points (if not set to 0 and used with displaySDP option).
     --displayMesh TEXT                    display a Mesh given in OFF or OFS format.
     --displayDigitalSurface               display the digital surface instead of display all the set of voxels (used with thresholdImage or displaySDP options)
     --colorizeCC                          colorize each Connected Components of the surface displayed by displayDigitalSurface option.
     -c,--colorSDP UINT x 4                set the color  discrete points: r g b a
     --colorMesh UINT x 4                  set the color of Mesh (given from displayMesh option) : r g b a
     -x,--scaleX FLOAT=1                   set the scale value in the X direction
     -y,--scaleY FLOAT=1                   set the scale value in the Y direction
     -z,--scaleZ FLOAT=1                   set the scale value in the Z direction
     --rescaleInputMin INT=0               min value used to rescale the input intensity (to avoid basic cast into 8  bits image).
     --rescaleInputMax INT=255             max value used to rescale the input intensity (to avoid basic cast into 8  bits image).
     -t,--transparency UINT=?              change the default transparency
  @endcode


   @b Example: 
   With the image display you can also threshold the image and display a set of voxel:  
   @code
   3dImageViewer -i $DGtal/examples/samples/lobster.vol --thresholdImage -m 180
   @endcode

   You should obtain such a result:

   @image html res3dImageViewer.png "resulting visualisation of 3d image with thresholded set of voxels."
 
   @see
   @ref 3dImageViewer.cpp

*/


int main( int argc, char** argv )
{

  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain,  unsigned char > Image3D;
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain,  unsigned char > Image2D;

  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string inputFileNameSDP;
  std::string inputFileNameMesh;

  DGtal::int64_t rescaleInputMin {0};
  DGtal::int64_t rescaleInputMax {255};
  bool grid {false};
  bool intergrid {false};
  bool emptyMode {false};
  bool displayDigitalSurface {false};
  bool thresholdImage {false};
  bool colorizeCC {false};
  int thresholdMin {0};
  int thresholdMax {255};
  std::vector<unsigned int> vectSDPIndex {0,1,2};
  std::vector<unsigned int> colorSDP;
  std::vector<unsigned int> colorMesh;

  float sx {1.0};
  float sy {1.0};
  float sz {1.0};
  double ballRadius = {0.0};
  unsigned char transp {255};

  
  app.description("Displays volume file as a voxel set by using QGLviewer\n 3dImageViewer -i $DGtal/examples/samples/lobster.vol --thresholdImage -m 180");
  
  app.add_option("-i,--input,1", inputFileName, "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
  ->required()
  ->check(CLI::ExistingFile);
  
  app.add_flag("--grid", grid , "draw slice images using grid mode.");
  app.add_flag("--intergrid", grid , "draw slice images using inter grid mode.");
  app.add_flag("--emptyMode", emptyMode,"remove the default boundingbox display.");
  app.add_flag("--thresholdImage", thresholdImage,"threshold the image to define binary shape");
  app.add_option("--thresholdMin", thresholdMin, "threshold min to define binary shape", true);
  app.add_option("--thresholdMax", thresholdMax, "threshold maw to define binary shape", true);
  app.add_option("--displaySDP,s",inputFileNameSDP, "display a set of discrete points (.sdp)" );
  app.add_option("--SDPindex", vectSDPIndex, "specify the sdp index.")
   ->expected(3);
  app.add_option("--SDPball",ballRadius, "use balls to display a set of discrete points (if not set to 0 and used with displaySDP option).", true);
  
  app.add_option("--displayMesh", inputFileNameMesh, "display a Mesh given in OFF or OFS format.");
  app.add_flag("--displayDigitalSurface",displayDigitalSurface, "display the digital surface instead of display all the set of voxels (used with thresholdImage or displaySDP options)" );
  
  app.add_flag("--colorizeCC", colorizeCC, "colorize each Connected Components of the surface displayed by displayDigitalSurface option.");
  app.add_option("--colorSDP,-c", colorSDP, "set the color  discrete points: r g b a ")
   ->expected(4);
  app.add_option("--colorMesh", colorMesh, "set the color of Mesh (given from displayMesh option) : r g b a ")
   ->expected(4);

  app.add_option("--scaleX,-x", sx, "set the scale value in the X direction", true );
  app.add_option("--scaleY,-y", sy, "set the scale value in the Y direction", true );
  app.add_option("--scaleZ,-z", sy, "set the scale value in the Z direction", true );
  app.add_option("--rescaleInputMin",rescaleInputMin, "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true );
  app.add_option("--rescaleInputMax",rescaleInputMax, "max value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true );
  app.add_option("--transparency,-t",transp, "change the default transparency", true );

  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------


  QApplication application(argc,argv);
  
  
  string extension = inputFileName.substr(inputFileName.find_last_of(".") + 1);
 
  Viewer3DImage<>::ModeVisu mode;
  if(emptyMode)
    mode=Viewer3DImage<>::Empty;
  else if(grid)
    mode=Viewer3DImage<>::Grid;
  else if(intergrid)
    mode=Viewer3DImage<>::InterGrid;
  else
    mode=Viewer3DImage<>::BoundingBox;

  Viewer3DImage<> viewer(mode);
  viewer.setWindowTitle("simple Volume Viewer");
  viewer.show();
  viewer.setGLScale(sx, sy, sz);

  typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
  Image3D image =  GenericReader< Image3D >::importWithValueFunctor( inputFileName,RescalFCT(rescaleInputMin,
                                                                                             rescaleInputMax,
                                                                                             0, 255) );
  Domain domain = image.domain();

  trace.info() << "Image loaded: "<<image<< std::endl;
  viewer.setVolImage(&image);
  
  // Used to display 3D surface
  Z3i::DigitalSet set3d(domain);

  viewer << Viewer3D<>::updateDisplay;
  if(thresholdImage){
    GradientColorMap<long> gradient( thresholdMin, thresholdMax);
    gradient.addColor(Color::Blue);
    gradient.addColor(Color::Green);
    gradient.addColor(Color::Yellow);
    gradient.addColor(Color::Red);
    for(Domain::ConstIterator it = domain.begin(), itend=domain.end(); it!=itend; ++it){
      unsigned char  val= image( (*it) );
      Color c= gradient(val);
      if(val<=thresholdMax && val >=thresholdMin)
      {
        if(!displayDigitalSurface)
        {
          viewer <<  CustomColors3D(Color((float)(c.red()), (float)(c.green()),(float)(c.blue()), transp),
                                    Color((float)(c.red()), (float)(c.green()),(float)(c.blue()), transp));
          viewer << *it;
        }
      }else
      {
        set3d.insert(*it);
      }
    }
  }

  if(inputFileNameSDP != "" ){
    if(colorSDP.size()==4){
      Color c(colorSDP[0], colorSDP[1], colorSDP[2], colorSDP[3]);
      viewer << CustomColors3D(c, c);
    }

    vector<Z3i::Point> vectVoxels;
    if(vectSDPIndex.size()==3)
    {
      vectVoxels = PointListReader<Z3i::Point>::getPointsFromFile(inputFileNameSDP, vectSDPIndex);
    }else
    {
      vectVoxels = PointListReader<Z3i::Point>::getPointsFromFile(inputFileNameSDP);
    }
    for(unsigned int i=0;i< vectVoxels.size(); i++)
    {
      if(!displayDigitalSurface)
      {
        if(ballRadius != 0.0)
        {
          viewer.addBall (vectVoxels.at(i), ballRadius);
        }
        else
        {
          viewer << vectVoxels.at(i);
        }
      }
      else
      {
        set3d.insert(vectVoxels.at(i));
      }
    }
  }

  if(inputFileNameMesh != "")
  {
    if(colorMesh.size() != 0)
    {
      Color c(colorMesh[0], colorMesh[1], colorMesh[2], colorMesh[3]);
      viewer.setFillColor(c);
    }
    DGtal::Mesh<Z3i::RealPoint> aMesh(colorMesh.size() == 0);
    MeshReader<Z3i::RealPoint>::importOFFFile(inputFileNameMesh, aMesh);
    viewer << aMesh;
  }

  if(displayDigitalSurface)
  {
    KSpace K;
    Point low = domain.lowerBound(); low[0]=low[0]-1; low[1]=low[1]-1; low[2]=low[2]-1;
    Point upp = domain.upperBound(); upp[0]=upp[0]+1; upp[1]=upp[1]+1; upp[2]=upp[2]+1;
    K.init(low, upp , true);
    SurfelAdjacency<3> SAdj( true );
    vector<vector<SCell> > vectConnectedSCell;
    trace.info() << "Extracting surface  set ... " ;
    Surfaces<KSpace>::extractAllConnectedSCell(vectConnectedSCell,K, SAdj, set3d, true);
    trace.info()<< " [done] " <<std::endl;
    GradientColorMap<long> gradient( 0, vectConnectedSCell.size());
    gradient.addColor(DGtal::Color::Red);
    gradient.addColor(DGtal::Color::Yellow);
    gradient.addColor(DGtal::Color::Green);
    gradient.addColor(DGtal::Color::Cyan);
    gradient.addColor(DGtal::Color::Blue);
    gradient.addColor(DGtal::Color::Magenta);
    gradient.addColor(DGtal::Color::Red);

    viewer << DGtal::SetMode3D(vectConnectedSCell.at(0).at(0).className(), "Basic");
    for(unsigned int i= 0; i <vectConnectedSCell.size(); i++)
    {
      for(unsigned int j= 0; j <vectConnectedSCell.at(i).size(); j++)
      {
        if(colorizeCC)
        {
          DGtal::Color c= gradient(i);
          viewer << CustomColors3D(Color(250, 0,0, transp), Color(c.red(),
                                                                  c.green(),
                                                                  c.blue(), transp));
        }else  if(colorSDP.size() != 0)
        {
          Color c(colorSDP[0], colorSDP[1], colorSDP[2], colorSDP[3]);
          viewer << CustomColors3D(c, c);
        }
        viewer << vectConnectedSCell.at(i).at(j);
      }
    }
  }
  
  viewer << Viewer3D<>::updateDisplay;
  DGtal::Z3i::Point size = image.domain().upperBound() - image.domain().lowerBound();
  DGtal::Z3i::Point center = image.domain().lowerBound()+size/2;
  unsigned int maxDist = std::max(std::max(size[2], size[1]), size[0]);
  viewer.camera()->setPosition(qglviewer::Vec(center[0],center[1], 
                                              center[2] + 2.0*maxDist));
  viewer.camera()->setSceneCenter(qglviewer::Vec(center[0],center[1],center[2]));
  return application.exec();
}
