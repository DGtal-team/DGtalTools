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
   --useLastCameraSetting                use the last camera setting of the user (i.e if a    .qglviewer.xml file is present in the current directory)

 
 @endcode
 
 
 @b Example:
 
 
 @code
 $ 3dVolViewer $DGtal/examples/samples/lobster.vol -m 60 -t 10
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

typedef ViewerSnap<> Viewer;



// call back function to display voxel coordinates
template<typename TImage>
int
displayCoordsCallBack( void* viewer, int name, void* data )
{
  TImage *image = (TImage *) data;
  std::stringstream ss;
  Z3i::Point p  =  DGtal::Linearizer<typename TImage::Domain>::getPoint(name, image->domain());
  // Check needed since simetimes the point appears outside (only in non debug mode). 
  if (image->domain().isInside(p)){
    ss << "Selected intensity:  " << (*image)(p) << "p " << p[0] << " "<< p[1] << " " << p[2] ;
    ((Viewer *) viewer)->displayMessage(QString(ss.str().c_str()), 100000);
  }
  return 0;
}


template <typename TImage>
void
processDisplay(ViewerSnap<> &viewer,  TImage &image,
               const typename TImage::Value &thresholdMin,
               const typename TImage::Value &thresholdMax,
               unsigned int numDisplayedMax,
               unsigned int transparency,
               bool useTransIntensity=false,
               bool interDisplay=false,
               bool squaredTransp=false)
{
  Domain domain = image.domain();
  GradientColorMap<typename TImage::Value> gradient( thresholdMin, thresholdMax);
  unsigned int numDisplayed = 0;
  gradient.addColor(Color::Blue);
  gradient.addColor(Color::Green);
  gradient.addColor(Color::Yellow);
  gradient.addColor(Color::Red);
  typename TImage::Value varInt = 0;
  if (useTransIntensity){
    varInt = thresholdMax-thresholdMin;
  }
  for(Domain::ConstIterator it = domain.begin(), itend=domain.end(); it!=itend; ++it){
    double trInt = 1.0;
    if (useTransIntensity){
      trInt = ((float)(image( (*it) )-thresholdMin)/(float)varInt);
    }
    typename TImage::Value val= image( (*it) );
    if(numDisplayed > numDisplayedMax)
      break;
    
    Color c= gradient(val);
    if(val<=thresholdMax && val >=thresholdMin)
    {
      
      viewer <<  CustomColors3D(Color(c.red(), c.green(),c.blue(),
                                      static_cast<unsigned int>(transparency*trInt*(squaredTransp? trInt : 1.0))),
                                Color(c.red(), c.green(),c.blue(),
                                      static_cast<unsigned int>(transparency*trInt*(squaredTransp? trInt: 1.0))));
      if (interDisplay)
      {
        auto p = *it;
        auto index = DGtal::Linearizer<typename TImage::Domain>::getIndex( p, domain);
        viewer << SetName3D( index ) ;
      }
      viewer << *it;
      numDisplayed++;
    }
  }
  if (interDisplay){
    viewer << SetSelectCallback3D(displayCoordsCallBack<TImage>,
                                  &image );
  }
}




int main( int argc, char** argv )
{
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  app.description("Display volume file as a voxel set by using QGLviewer. \n Example: \n \t 3dVolViewer  $DGtal/examples/samples/lobster.vol -m 60 -t 10");
  std::string inputFileName;
  DGtal::int64_t rescaleInputMin {0};
  DGtal::int64_t rescaleInputMax {255};
  double thresholdMin {0};
  double thresholdMax {255};
  unsigned int transparency {255};
  unsigned int numDisplayedMax {500000};
  std::string displayMesh;
  std::string snapShotFile;
  std::vector<unsigned int> colorMesh;
  string inputType {""};
  bool interactiveDisplayVoxCoords {false};
  bool transIntensity {false};
  bool transIntensitySq {false};
  bool useLastCamSet {false};

  app.add_option("-i,--input,1", inputFileName, "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
  ->required()
  ->check(CLI::ExistingFile);
#ifdef WITH_ITK
  app.add_option("--inputType", inputType, "to specify the input image type (int or double).")
    -> check(CLI::IsMember({"int", "double"}));
#endif 
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
  app.add_flag("--transIntensity",transIntensity , "Used vocel intensity to define transparency valeue");
  app.add_flag("--transIntensitySq",transIntensitySq , "Used squared vocel intensity to define transparency valeue");
  app.add_flag("--interactiveDisplayVoxCoords,-c", interactiveDisplayVoxCoords, " by using this option the coordinates can be displayed after selection (shift+left click on voxel).");
 app.add_flag("--useLastCameraSetting", useLastCamSet, "use the last camera setting of the user (i.e if a .qglviewer.xml file is present in the current directory)");
    app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------
  
  QApplication application(argc,argv);
  
  Viewer viewer(snapShotFile != "");
  if(snapShotFile != ""){
    viewer.setSnapshotFileName(QString(snapShotFile.c_str()));
  }
  viewer.setWindowTitle("simple Volume Viewer");
  viewer.show();
  typedef ImageContainerBySTLVector < Z3i::Domain,  double > Image3D_D;
  typedef ImageContainerBySTLVector < Z3i::Domain,  int > Image3D_I;
  typedef ImageSelector<Domain, unsigned char>::Type Image;

  string extension = inputFileName.substr(inputFileName.find_last_of(".") + 1);
  
  std::vector<double> vectValD;
  std::vector<int> vectValI;
  std::vector<unsigned char> vectValUC;
  // Image of different types are pre constructed here else it will be deleted after the type selection
  // (and it is used in display callback)
  Z3i::Domain d;
  Image3D_D  imageD = Image3D_D(d);
  Image3D_I  imageI = Image3D_I(d);
  Image  image = Image(d);

  if(extension != "sdp")
  {
    unsigned int numDisplayed=0;

#ifdef WITH_ITK
  if (inputType=="double")
  {
    imageD = DGtal::GenericReader<Image3D_D>::import(inputFileName);
    trace.info() << "[done]"<< std::endl;
    trace.info() << "Image loaded:  D "<<imageD<< std::endl;
    processDisplay(viewer, imageD, thresholdMin, thresholdMax, numDisplayedMax, transparency,
                   transIntensity||transIntensitySq, interactiveDisplayVoxCoords, transIntensitySq);
  }
  else if (inputType=="int")
  {
    imageI= DGtal::GenericReader<Image3D_I>::import(inputFileName);
    trace.info() << "Image loaded: "<<image<< std::endl;
    processDisplay(viewer, imageI, (int)thresholdMin, (int)thresholdMax, numDisplayedMax, transparency,
                   transIntensity||transIntensitySq, interactiveDisplayVoxCoords, transIntensitySq);
  } else {
    typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
    image =  GenericReader< Image >::importWithValueFunctor( inputFileName,RescalFCT(rescaleInputMin,
                                                                                     rescaleInputMax,
                                                                                     0, 255) );
    trace.info() << "Image loaded: "<<image<< std::endl;
    processDisplay(viewer, image, thresholdMin, thresholdMax, numDisplayedMax, transparency,
                   transIntensity||transIntensitySq, interactiveDisplayVoxCoords, transIntensitySq);
  }
#else
    typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
    image =  GenericReader< Image >::importWithValueFunctor( inputFileName,RescalFCT(rescaleInputMin,
                                                                                     rescaleInputMax,
                                                                                     0, 255) );
    trace.info() << "Image loaded: "<<image<< std::endl;
    processDisplay(viewer, image, thresholdMin, thresholdMax, numDisplayedMax, transparency,
                   transIntensity||transIntensitySq, interactiveDisplayVoxCoords, transIntensitySq);
#endif
  }
  else if(extension=="sdp")
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
    
  if (useLastCamSet)
  {
        viewer.restoreStateFromFile();
  }
  // First display transparency improve
  viewer.sortTriangleFromCamera();
  viewer.sortQuadFromCamera();
  viewer.sortSurfelFromCamera();
  viewer.sortPolygonFromCamera();
  viewer  << Viewer::updateDisplay;
  trace.info() << "[display ready]"<< std::endl;
  return application.exec();
}
