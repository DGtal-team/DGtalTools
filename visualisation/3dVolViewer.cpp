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
 * @file 3dVolViewer.cpp
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Lorraine, France
 *
 * @date 2011/01/04
 *
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/readers/MeshReader.h"

#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/images/ImageSelector.h"

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
 @page Doc3dVolViewer 3dVolViewer
 
 @brief Displays volume file as a voxel set by using QGLviewer.
 
 The mode  specifies if you wish to see surface elements (BDRY), the inner
 voxels (INNER) or the outer voxels (OUTER) that touch the boundary.
 
 @b Usage:   3dVolViewer [input]
 
 @b Allowed @b options @b are :
 
 @code
 

 Positionals:
   1 TEXT:FILE REQUIRED                  vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
   -m,--thresholdMin INT=0               threshold min (excluded) to define binary shape.
   -M,--thresholdMax INT=255             threshold max (included) to define binary shape.
   --rescaleInputMin INT=0               min value used to rescale the input intensity (to avoid basic cast into 8  bits image).
   --rescaleInputMax INT=255             max value used to rescale the input intensity (to avoid basic cast into 8  bits image).
   -n,--numMaxVoxel UINT                 set the maximal voxel number to be displayed.
   --displayMesh TEXT                    display a Mesh given in OFF or OFS format.
   --colorMesh UINT ...                  set the color of Mesh (given from displayMesh option) : r g b a
   -d,--doSnapShotAndExit                save display snapshot into file. Notes that the camera setting is set by default according the last saved configuration (use SHIFT+Key_M to save current camera setting in the Viewer3D). If the camera setting was not saved it will use the default camera setting.
   -t,--transparency UINT=255            change the defaukt transparency
   
 
 @endcode
 
 
 @b Example:
 
 
 @code
 $ 3dVolViewer -i $DGtal/examples/samples/lobster.vol -m 60 -t 10
 @endcode
 
 You should obtain such a result:
 
 @image html res3dVolViewer.png "Resulting visualization."
 
 
 @see
 @ref 3dVolViewer.cpp
 
 */




template < typename Space = DGtal::Z3i::Space, typename KSpace = DGtal::Z3i::KSpace>
struct ViewerSnap: DGtal::Viewer3D <Space, KSpace>
{
  
  ViewerSnap(bool saveSnap): Viewer3D<Space, KSpace>(), mySaveSnap(saveSnap){
  };
  
  virtual  void
  init(){
    DGtal::Viewer3D<>::init();
    if(mySaveSnap){
      QObject::connect(this, SIGNAL(drawFinished(bool)), this, SLOT(saveSnapshot(bool)));
    }
  };
  bool mySaveSnap;
};





int main( int argc, char** argv )
{
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  app.description("Display volume file as a voxel set by using QGLviewer. \n Example: \n \t 3dVolViewer -i $DGtal/examples/samples/lobster.vol -m 60 -t 10");
  std::string inputFileName;
  DGtal::int64_t rescaleInputMin {0};
  DGtal::int64_t rescaleInputMax {255};
  int thresholdMin {0};
  int thresholdMax {255};
  unsigned int transparency {255};
  unsigned int numDisplayedMax {500000};
  std::string displayMesh;
  std::string snapShotFile;
  std::vector<unsigned int> colorMesh;
  app.add_option("-i,--input,1", inputFileName, "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
  ->required()
  ->check(CLI::ExistingFile);
  
  app.add_option("--thresholdMin,-m", thresholdMin, "threshold min (excluded) to define binary shape.", true);
  app.add_option("--thresholdMax,-M", thresholdMax, "threshold max (included) to define binary shape.", true);
  app.add_option("--rescaleInputMin", rescaleInputMin, "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);
  app.add_option("--rescaleInputMax", rescaleInputMax, "max value used to rescale the input intensity (to avoid basic cast into 8  bits image).", true);
  app.add_option("--numMaxVoxel,-n", numDisplayedMax, "set the maximal voxel number to be displayed.");
  app.add_option("--displayMesh", displayMesh, "display a Mesh given in OFF or OFS format.");
  app.add_option("--colorMesh", colorMesh, "set the color of Mesh (given from displayMesh option) : r g b a ")
   ->expected(4);
  app.add_flag("--doSnapShotAndExit,-d",snapShotFile, "save display snapshot into file. Notes that the camera setting is set by default according the last saved configuration (use SHIFT+Key_M to save current camera setting in the Viewer3D). If the camera setting was not saved it will use the default camera setting." );
  
  app.add_option("--transparency,-t", transparency, "change the defaukt transparency", true);
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------
  
  
  QApplication application(argc,argv);
  typedef ViewerSnap<> Viewer;
  
  Viewer viewer(snapShotFile != "");
  if(snapShotFile != ""){
    viewer.setSnapshotFileName(QString(snapShotFile.c_str()));
  }
  
  viewer.setWindowTitle("simple Volume Viewer");
  viewer.show();
  
  typedef ImageSelector<Domain, unsigned char>::Type Image;
  string extension = inputFileName.substr(inputFileName.find_last_of(".") + 1);
  if(extension != "sdp")
  {
    unsigned int numDisplayed=0;
    
    typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
    Image image =  GenericReader< Image >::importWithValueFunctor( inputFileName,RescalFCT(rescaleInputMin,
                                                                                           rescaleInputMax,
                                                                                           0, 255) );
    
    trace.info() << "Image loaded: "<<image<< std::endl;
    Domain domain = image.domain();
    GradientColorMap<long> gradient( thresholdMin, thresholdMax);
    gradient.addColor(Color::Blue);
    gradient.addColor(Color::Green);
    gradient.addColor(Color::Yellow);
    gradient.addColor(Color::Red);
    for(Domain::ConstIterator it = domain.begin(), itend=domain.end(); it!=itend; ++it){
      unsigned char  val= image( (*it) );
      if(numDisplayed > numDisplayedMax)
        break;
      Color c= gradient(val);
      if(val<=thresholdMax && val >=thresholdMin)
      {
        viewer <<  CustomColors3D(Color((float)(c.red()), (float)(c.green()),(float)(c.blue()), transparency),
                                  Color((float)(c.red()), (float)(c.green()),(float)(c.blue()), transparency));
        viewer << *it;
        numDisplayed++;
      }
    }
  }else if(extension=="sdp")
  {
    vector<Z3i::RealPoint> vectVoxels = PointListReader<Z3i::RealPoint>::getPointsFromFile(inputFileName);
    for(unsigned int i=0;i< vectVoxels.size(); i++){
      viewer << Z3i::Point(vectVoxels.at(i), functors::Round<>());
    }
  }
  if(displayMesh != "")
  {
    if(colorMesh.size() != 0)
    {
      Color c(colorMesh[0], colorMesh[1], colorMesh[2], colorMesh[3]);
      viewer.setFillColor(c);
    }
    
    DGtal::Mesh<Z3i::RealPoint> aMesh(colorMesh.size() == 0);
    MeshReader<Z3i::RealPoint>::importOFFFile(displayMesh, aMesh);
    viewer << aMesh;
  }
  
  viewer << Viewer3D<>::updateDisplay;
  if(snapShotFile != "")
  {
    // Appy cleaning just save the last snap
    if(!viewer.restoreStateFromFile())
    {
      viewer.update();
    }
    std::string extension = snapShotFile.substr(snapShotFile.find_last_of(".") + 1);
    std::string basename = snapShotFile.substr(0, snapShotFile.find_last_of("."));
    for(int i=0; i< viewer.snapshotCounter()-1; i++){
      std::stringstream s;
      s << basename << "-"<< setfill('0') << setw(4)<<  i << "." << extension;
      trace.info() << "erase temp file: " << s.str() << std::endl;
      remove(s.str().c_str());
    }
    std::stringstream s;
    s << basename << "-"<< setfill('0') << setw(4)<<  viewer.snapshotCounter()-1 << "." << extension;
    rename(s.str().c_str(), snapShotFile.c_str());
    return 0;
  }
  
  return application.exec();
}
