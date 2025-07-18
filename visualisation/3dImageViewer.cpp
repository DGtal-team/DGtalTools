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
 * An example file named 3dImageViewer.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/viewers/PolyscopeViewer.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/SurfelAdjacency.h"

#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#ifdef DGTAL_WITH_ITK
#include "DGtal/io/readers/DicomReader.h"
#endif

#include "DGtal/images/ImageSelector.h"



// #include "specificClasses/Viewer3DImage.cpp"

#include "CLI11.hpp"



using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
   @page Doc3dImageViewer 3dImageViewer
 
   @brief Displays volume file as a voxel set by using PolyscopeViewer.

   @b Usage:  3dImageViewer [OPTIONS] 1 [s]

   @b Allowed @b options @b are :
 
   @code

   Positionals:
     1 TEXT:FILE REQUIRED                  vol file (.vol, .longvol .p3d, .pgm3d and if DGTAL_WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
     s TEXT                                display a set of discrete points (.sdp)

   Options:
     -h,--help                             Print this help message and exit
     -i,--input TEXT:FILE REQUIRED         vol file (.vol, .longvol .p3d, .pgm3d and if DGTAL_WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255.
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
   3dImageViewer $DGtal/examples/samples/lobster.vol --thresholdImage -m 180
   @endcode

   You should obtain such a result:

   @image html res3dImageViewer.png "resulting visualisation of 3d image with thresholded set of voxels."
 
   @see
   @ref 3dImageViewer.cpp

*/




enum TypeSlice {SliceX, SliceY, SliceZ};


typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain,  unsigned char > Image3D;
typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain,  unsigned char > Image2D;
typedef DGtal::ConstImageAdapter<Image3D, Image2D::Domain, DGtal::functors::Projector< DGtal::Z3i::Space>,
                                Image3D::Value,  DGtal::functors::Identity >  SliceImageAdapter;
typedef DGtal::ConstImageAdapter<Image3D, DGtal::Z2i::Domain, DGtal::functors::SliceRotator2D< DGtal::Z3i::Domain >,
Image3D::Value,  DGtal::functors::Identity >  MyRotatorSliceImageAdapter;

static polyscope::SurfaceTextureScalarQuantity* texSliceX = nullptr;
static polyscope::SurfaceTextureScalarQuantity* texSliceY = nullptr;
static polyscope::SurfaceTextureScalarQuantity* texSliceZ = nullptr;

static int sliceXNum = 100;
static int sliceYNum = 100;
static int sliceZNum = 25;

static PolyscopeViewer<> *refViewer;
static polyscope::SurfaceTextureScalarQuantity* g_qScalarSliceX = nullptr;
static polyscope::SurfaceTextureScalarQuantity* g_qScalarSliceY = nullptr;
static polyscope::SurfaceTextureScalarQuantity* g_qScalarSliceZ = nullptr;

static polyscope::SurfaceTextureScalarQuantity* g_qScalar[] = {g_qScalarSliceX, g_qScalarSliceY, g_qScalarSliceZ};
Image3D image  = Image3D(DGtal::Z3i::Domain());
DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctorX;
DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctorY;
DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctorZ;


polyscope::SurfaceMesh *slicePlaneX;
polyscope::SurfaceMesh *slicePlaneY;
polyscope::SurfaceMesh *slicePlaneZ;


std::vector<glm::vec3>
getVertices(TypeSlice aTypeSlice, const Z3i::Point &ptInf, const Z3i::Point &ptSup){
    switch (aTypeSlice) {
        case TypeSlice::SliceX:
            return  {
                {ptInf[0],    ptInf[1],   ptInf[2]   }, // bas-gauche
                {ptInf[0],    ptSup[1],   ptInf[2]   }, // bas-droit
                {ptInf[0],    ptSup[1],   ptSup[2]   }, // haut-droit
                {ptInf[0],    ptInf[1],   ptSup[2]   }, // haut-gauche
            };
        case TypeSlice::SliceY:
            return  {
                {ptInf[0],    ptInf[1],   ptInf[2]   }, // bas-gauche
                {ptSup[0],    ptInf[1],   ptInf[2]   }, // bas-droit
                {ptSup[0],    ptInf[1],   ptSup[2]   }, // haut-droit
                {ptInf[0],    ptInf[1],   ptSup[2]   }, // haut-gauche
            };
        case TypeSlice::SliceZ:
            return  {
                {ptInf[0],    ptInf[1],   ptInf[2]   }, // bas-gauche
                {ptSup[0],    ptInf[1],   ptInf[2]   }, // bas-droit
                {ptSup[0],    ptSup[1],   ptInf[2]   }, // haut-droit
                {ptInf[0],    ptSup[1],   ptInf[2]   }, // haut-gauche
            };
    }
    return  {
        {ptInf[0],    ptInf[1],   ptInf[2]   }, // bas-gauche
        {ptInf[0],    ptInf[1],   ptSup[2]   }, // bas-droit
        {ptSup[0],    ptInf[1],   ptSup[2]   }, // haut-droit
        {ptSup[0],    ptInf[1],   ptInf[2]   }, // haut-gauche
    };
            
}


polyscope::SurfaceMesh *
initSlices(const MyRotatorSliceImageAdapter &sliceIm, string name, TypeSlice aTypeSlice) {
    polyscope::SurfaceMesh * res;
    auto ptInf = sliceIm.sourceDomainPoint(sliceIm.domain().lowerBound());
    auto ptSup = sliceIm.sourceDomainPoint(sliceIm.domain().upperBound());
    auto dim = sliceIm.domain().upperBound() - sliceIm.domain().lowerBound();
    auto dimX = dim[0];
    auto dimY = dim[1];
    auto maxDim = max(dimX, dimY);
    std::vector<glm::vec3> vertices = getVertices(aTypeSlice, ptInf, ptSup);
    std::vector<std::vector<size_t>> faces = {{0, 1, 2, 3}};
    res = polyscope::registerSurfaceMesh(name, vertices, faces);
    std::vector<glm::vec2> param = {
        {0.0f, 0.0f},
        {1.0f, 0.0f},
        {1.0f, 1.0f},
        {0.0f, 1.0f}
    };

    auto qParam = res->addParameterizationQuantity("param", param);
    std::vector<float> valuesTex;
    for(unsigned int y =  0; y< dimY; y++){
        for(unsigned int x =  0; x< dimX; x++){
            valuesTex.push_back(((float)sliceIm(Z2i::Point(x,y))));
        }
    }
     
    float minV = *std::min_element(valuesTex.begin(), valuesTex.end());
    float maxV = *std::max_element(valuesTex.begin(), valuesTex.end());
    for (float& v : valuesTex) {
        v = (v - minV) / (maxV - minV );
    }
    g_qScalar[aTypeSlice] = res->addTextureScalarQuantity("tScalar", *qParam, dimX, dimY,
                                              valuesTex, polyscope::ImageOrigin::LowerLeft);
    g_qScalar[aTypeSlice]->setFilterMode(polyscope::FilterMode::Nearest); // change filter for sampling
    g_qScalar[aTypeSlice]->setEnabled(true);
    return res;
}



void
updateSlices(const MyRotatorSliceImageAdapter &sliceIm, string name, TypeSlice aTypeSlice) {
    polyscope::SurfaceMesh * sm = polyscope::getSurfaceMesh(name);
    auto ptInf = sliceIm.sourceDomainPoint(sliceIm.domain().lowerBound());
    auto ptSup = sliceIm.sourceDomainPoint(sliceIm.domain().upperBound());
    auto dim = sliceIm.domain().upperBound() - sliceIm.domain().lowerBound();
    auto dimX = dim[0];
    auto dimY = dim[1];
    std::vector<float> valuesTex;
    for(unsigned int y =  0; y< dimY; y++){
        for(unsigned int x =  0; x< dimX; x++){
            valuesTex.push_back(((float)sliceIm(Z2i::Point(x,y))));
        }
    }
    
    float minV = *std::min_element(valuesTex.begin(), valuesTex.end());
    float maxV = *std::max_element(valuesTex.begin(), valuesTex.end());
    for (float& v : valuesTex) {
        v = (v - minV) / (maxV - minV );
    }
    g_qScalar[aTypeSlice] -> updateData(valuesTex);
    sm->updateVertexPositions(getVertices(aTypeSlice, ptInf, ptSup));

}

void callbackFaceID() {
    srand((unsigned) time(NULL));
    ImGui::Begin("Editing tools");
    ImGui::Text("Setting selection size:");
    if (ImGui::Button("+X"))
        {
         
        }
    if (ImGui::SliderInt("slice x", &sliceXNum, image.domain().lowerBound()[0], image.domain().upperBound()[0], "size = %i")){
        std::cout << "val "<<sliceXNum << std::endl;
       
        DGtal::functors::SliceRotator2D<DGtal::Z3i::Domain> aSliceFunctorX(0, image.domain(), sliceXNum, 2, 0 );
                const DGtal::functors::Identity identityFunctor{};
        DGtal::Z2i::Domain domain2DX(invFunctorX(image.domain().lowerBound()),
                                     invFunctorX(image.domain().upperBound()));

        MyRotatorSliceImageAdapter sliceImageX( image, domain2DX, aSliceFunctorX, identityFunctor );
        updateSlices(sliceImageX, "slicex", TypeSlice::SliceX);
   

    }
    if (ImGui::SliderInt("slice y", &sliceYNum, image.domain().lowerBound()[1], image.domain().upperBound()[1], "size = %i")){
        std::cout << "val "<<sliceYNum << std::endl;
        DGtal::functors::SliceRotator2D<DGtal::Z3i::Domain> aSliceFunctorY(1, image.domain(), sliceYNum, 2, 0 );
        DGtal::Z2i::Domain domain2DY(invFunctorY(image.domain().lowerBound()),
                                     invFunctorY(image.domain().upperBound()));
        const DGtal::functors::Identity identityFunctor{};
        MyRotatorSliceImageAdapter sliceImageY( image, domain2DY, aSliceFunctorY, identityFunctor );
        updateSlices(sliceImageY, "slicey", TypeSlice::SliceY);
    }
    if (ImGui::SliderInt("slice z", &sliceZNum, image.domain().lowerBound()[2], image.domain().upperBound()[2], "size = %i")){
        std::cout << "val "<<sliceZNum << std::endl;
        DGtal::functors::SliceRotator2D<DGtal::Z3i::Domain> aSliceFunctorZ(2, image.domain(), sliceZNum, 0, 0 );
        DGtal::Z2i::Domain domain2DZ(invFunctorZ(image.domain().lowerBound()),
                                     invFunctorZ(image.domain().upperBound()));
        const DGtal::functors::Identity identityFunctor{};
        MyRotatorSliceImageAdapter sliceImageZ( image, domain2DZ, aSliceFunctorZ, identityFunctor );
        updateSlices(sliceImageZ, "slicez", TypeSlice::SliceZ);

        
   

    }

    ImGui::End();
}



int main( int argc, char** argv )
{

 
    
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

  
  app.description("Displays volume file as a voxel set by using PolyscopeViewer\n 3dImageViewer  $DGtal/examples/samples/lobster.vol --thresholdImage -m 180");
  
  app.add_option("-i,--input,1", inputFileName, "vol file (.vol, .longvol .p3d, .pgm3d and if DGTAL_WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
  ->required()
  ->check(CLI::ExistingFile);
  
  app.add_flag("--grid", grid , "draw slice images using grid mode.");
  app.add_flag("--intergrid", grid , "draw slice images using inter grid mode.");
  app.add_flag("--emptyMode", emptyMode,"remove the default boundingbox display.");
  app.add_flag("--thresholdImage", thresholdImage,"threshold the image to define binary shape");
  app.add_option("--thresholdMin", thresholdMin, "threshold min to define binary shape");
  app.add_option("--thresholdMax", thresholdMax, "threshold maw to define binary shape");
  app.add_option("--displaySDP,s",inputFileNameSDP, "display a set of discrete points (.sdp)" );
  app.add_option("--SDPindex", vectSDPIndex, "specify the sdp index.")
   ->expected(3);
  app.add_option("--SDPball",ballRadius, "use balls to display a set of discrete points (if not set to 0 and used with displaySDP option).");
  
  app.add_option("--displayMesh", inputFileNameMesh, "display a Mesh given in OFF or OFS format.");
  app.add_flag("--displayDigitalSurface",displayDigitalSurface, "display the digital surface instead of display all the set of voxels (used with thresholdImage or displaySDP options)" );
  
  app.add_flag("--colorizeCC", colorizeCC, "colorize each Connected Components of the surface displayed by displayDigitalSurface option.");
  app.add_option("--colorSDP,-c", colorSDP, "set the color  discrete points: r g b a ")
   ->expected(4);
  app.add_option("--colorMesh", colorMesh, "set the color of Mesh (given from displayMesh option) : r g b a ")
   ->expected(4);

  app.add_option("--scaleX,-x", sx, "set the scale value in the X direction" );
  app.add_option("--scaleY,-y", sy, "set the scale value in the Y direction" );
  app.add_option("--scaleZ,-z", sy, "set the scale value in the Z direction" );
  app.add_option("--rescaleInputMin",rescaleInputMin, "min value used to rescale the input intensity (to avoid basic cast into 8  bits image)." );
  app.add_option("--rescaleInputMax",rescaleInputMax, "max value used to rescale the input intensity (to avoid basic cast into 8  bits image)." );
  app.add_option("--transparency,-t",transp, "change the default transparency" );

  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

    polyscope::options::programName = "3dImageViewer - (DGtalTools)";
    polyscope::options::buildGui=false;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    
    PolyscopeViewer<> viewer;
    refViewer = &viewer;
    
    string extension = inputFileName.substr(inputFileName.find_last_of(".") + 1);
   

    invFunctorX.initRemoveOneDim(0);
    invFunctorY.initRemoveOneDim(1);
    invFunctorZ.initRemoveOneDim(2);

  image =  GenericReader< Image3D >::importWithValueFunctor( inputFileName,RescalFCT(rescaleInputMin,
                                                                                             rescaleInputMax,
                                                                                             0,
                                                                                             255));
  Domain domain = image.domain();

  trace.info() << "Image loaded: "<<image<< std::endl;
  //viewer << image;

    // Used to display 3D surface
  Z3i::DigitalSet set3d(domain);

  if(thresholdImage){
    viewer.newCubeList("Threshold image");
    viewer.allowReuseList = true;
    for(Domain::ConstIterator it = domain.begin(), itend=domain.end(); it!=itend; ++it){
      unsigned char  val= image( (*it) );
      if(val<=thresholdMax && val >=thresholdMin)
      {
        if(!displayDigitalSurface)
        {
          viewer << WithQuantity(*it, "value", val);
        }
      }else
      {
        set3d.insert(*it);
      }
    }
    viewer.endCurrentGroup();
  }
    
    
    
  auto  myImageOrigin = image.domain().lowerBound();
    
 auto   mySliceXPos=myImageOrigin[0];
 auto   mySliceYPos=myImageOrigin[1];
 auto   mySliceZPos=myImageOrigin[2];

   // Adding X slice in the viewer.
    DGtal::Z2i::Domain domain2DX(invFunctorX(image.domain().lowerBound()),
                     invFunctorX(image.domain().upperBound()));
    DGtal::functors::SliceRotator2D<DGtal::Z3i::Domain> aSliceFunctorX(0, image.domain(), sliceXNum, 2, 0 );
    const DGtal::functors::Identity identityFunctor{};
    MyRotatorSliceImageAdapter sliceImageX( image, domain2DX, aSliceFunctorX, identityFunctor );
    slicePlaneX = initSlices(sliceImageX, "slicex", TypeSlice::SliceX);

    // Adding Y slice in the viewer.
     DGtal::Z2i::Domain domain2DY(invFunctorY(image.domain().lowerBound()),
                      invFunctorY(image.domain().upperBound()));
     DGtal::functors::SliceRotator2D<DGtal::Z3i::Domain> aSliceFunctorY(1, image.domain(), sliceYNum, 0, 0 );
     MyRotatorSliceImageAdapter sliceImageY( image, domain2DY, aSliceFunctorY, identityFunctor );
     slicePlaneY = initSlices(sliceImageY, "slicey", TypeSlice::SliceY);
   
    // Adding Z slice in the viewer.
     DGtal::Z2i::Domain domain2DZ(invFunctorZ(image.domain().lowerBound()),
                      invFunctorZ(image.domain().upperBound()));
     DGtal::functors::SliceRotator2D<DGtal::Z3i::Domain> aSliceFunctorZ(2, image.domain(), sliceZNum, 1, 0 );
     MyRotatorSliceImageAdapter sliceImageZ( image, domain2DZ, aSliceFunctorZ, identityFunctor );
     slicePlaneZ = initSlices(sliceImageZ, "slicez", TypeSlice::SliceZ);
    
  if(inputFileNameSDP != "" ){
    if(colorSDP.size()==4){
      Color c(colorSDP[0], colorSDP[1], colorSDP[2], colorSDP[3]);
      viewer << c;
    }

    vector<Z3i::Point> vectVoxels;
    if(vectSDPIndex.size()==3)
    {
      vectVoxels = PointListReader<Z3i::Point>::getPointsFromFile(inputFileNameSDP, vectSDPIndex);
    }else
    {
      vectVoxels = PointListReader<Z3i::Point>::getPointsFromFile(inputFileNameSDP);
    }

    if (ballRadius != 0.0) viewer.drawAsBalls();
    for(unsigned int i=0;i< vectVoxels.size(); i++)
    {
      if(!displayDigitalSurface)
      {
          viewer << vectVoxels.at(i);
      }
      else
      {
        set3d.insert(vectVoxels.at(i));
      }
    }
    viewer.drawAsPaving();
  }

  if(inputFileNameMesh != "")
  {
    if(colorMesh.size() != 0)
    {
      Color c(colorMesh[0], colorMesh[1], colorMesh[2], colorMesh[3]);
      viewer.drawColor(c);
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

    viewer.drawAsSimplified();
    for(unsigned int i= 0; i <vectConnectedSCell.size(); i++)
    {
      for(unsigned int j= 0; j <vectConnectedSCell.at(i).size(); j++)
      {
        const auto& toDraw = vectConnectedSCell.at(i).at(j);
        if(colorizeCC) 
        {
          viewer << WithQuantity(toDraw, "index", i);
        } 
        else if(colorSDP.size() != 0)
        {
          Color c(colorSDP[0], colorSDP[1], colorSDP[2], colorSDP[3]);
          viewer << WithQuantity(toDraw, "color", c);
        } 
        else 
        {
          viewer << toDraw; 
        }
      }
    }
  }
  
  DGtal::Z3i::Point size = image.domain().upperBound() - image.domain().lowerBound();
  DGtal::Z3i::Point center = image.domain().lowerBound()+size/2;
  unsigned int maxDist = std::max(std::max(size[2], size[1]), size[0]);
  polyscope::state::userCallback = callbackFaceID;

  viewer.show();
  return 0;
}
