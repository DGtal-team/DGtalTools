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
 * @file mesh2heightfield.cpp
 *  @ingroup Converters
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2015/04/08
 *
 *
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include <DGtal/base/BasicFunctors.h>
#include "DGtal/math/linalg/SimpleMatrix.h"

#include <DGtal/math/linalg/SimpleMatrix.h>
#include <DGtal/math/linalg/EigenDecomposition.h>



#include "CLI11.hpp"



using namespace std;
using namespace DGtal;



/**
 @page mesh2heightfield mesh2heightfield
 @brief  Converts a mesh file into a projected 2D image given from a normal direction N and from a starting point P.
 @ingroup convertertools
 
 The 3D mesh is discretized and scanned in the normal direction N, starting from P with a step 1.


@b Usage: mesh2heightfield [input] [output]

@b Allowed @b options @b are:

@code

 Positionals:
   1 TEXT:FILE                           input heightfield file (2D image).
   2 TEXT                                output image.

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE                  input heightfield file (2D image).
   -d,--domain UINT x 2                  specify the domain (required when normal are imported and if --inout is not given).
   -o,--output TEXT                      output image.
   --importNormal TEXT                   import normals from file.
   --orderedNormalsImport                Use ordered normals.
   --lightDir,--lDir,--ld FLOAT x 3      light source direction: lx ly lz.
   --lightPos,--lPos,--lp FLOAT x 3      light source position: px py pz.
   -s,--specularModel FLOAT x 3          use specular Nayar model with 3 param Kdiff, Kspec, sigma.
   -r,--reflectanceMap TEXT:FILE         specify a image as reflectance map.
   --hsvShading                          use shading with HSV shading (given from the normal vector)
   --normalMap                           generates normal map.
   -v,--invertNormals                    invert normal orientations.


@endcode

@b Example:
@code
  $ mesh2heightfield -i ${DGtal}/examples/samples/tref.off  heighfield.pgm

@endcode
You will obtain such image:
@image html  resMesh2heightfield.png "Resulting heightfield."
@see mesh2heightfield.cpp

 */
typedef PointVector<9, double> Matrix3x3Point;
typedef SimpleMatrix<double, 3, 3 > CoVarianceMat;

/**
 * Order natural (book reading from upper left to lower right)
 **/
static CoVarianceMat
getCoVarianceMatFrom(const Matrix3x3Point &aMatrixPt){
  CoVarianceMat res;
  for(unsigned int i = 0; i<3; i++){
    for(unsigned int j = 0; j<3; j++){
      res.setComponent(i, j, aMatrixPt[j+i*3] );
    }
  }
  return res;
}


std::pair<DGtal::Z3i::RealPoint, DGtal::Z3i::RealPoint>
getMainDirsCoVar(Mesh<Z3i::RealPoint>::Iterator begin, Mesh<Z3i::RealPoint>::Iterator end )
{
  std::pair<DGtal::Z3i::RealPoint, DGtal::Z3i::RealPoint> res;
  Matrix3x3Point c ;
  unsigned int nb = 0;
  //compute centroïd of point cloud
  Z3i::RealPoint centroid;
  int cp=0;
  for(auto it=begin; it != end; it++){
      Z3i::RealPoint pt=*it;
      centroid[0] += pt[0];
      centroid[1] += pt[1];
      centroid[2] += pt[2];
       ++cp;
  }
  centroid /=cp;
  //compute matrix like PCL : https://pointclouds.org/documentation/moment__of__inertia__estimation_8hpp_source.html
  for(auto it=begin; it != end; it++){

      Z3i::RealPoint pt=*it;
      pt[0] -= centroid[0];
      pt[1] -= centroid[1];
      pt[2] -= centroid[2];

      trace.info()<<pt<< std::endl;
      c[4] += pt[1] * pt[1];
      c[5] += pt[1] * pt[2];
      c[8] += pt[2] * pt[2];
      pt *= pt[0];
      c[0] += pt[0];
      c[1] += pt[1];
      c[2] += pt[2];

      nb++;

  }
  c[3] = c[1];
  c[6] = c[2];
  c[7] = c[5];
  //normalize matrix
  c = c/nb;
  CoVarianceMat covar = getCoVarianceMatFrom(c);
  SimpleMatrix<double, 3, 3 > eVects;
  PointVector<3, double> eVals;
  DGtal::EigenDecomposition<3, double, CoVarianceMat>::getEigenDecomposition (covar, eVects, eVals);
  unsigned int temp = 0;
  unsigned int major_index = 0;
  unsigned int middle_index = 1;
  unsigned int minor_index = 2;
  //reorder first and second axis according to their eigen value
  if (eVals(major_index) < eVals(middle_index))
  {
     temp = major_index;
     major_index = middle_index;
     middle_index = temp;
  }
  if (eVals(major_index) < eVals(minor_index))
  {
    temp = major_index;
    major_index = minor_index;
    minor_index = temp;
  }
  if (eVals(middle_index) < eVals(minor_index))
  {
    temp = minor_index;
    minor_index = middle_index;
    middle_index = temp;
  }

  res.first = eVects.column(major_index);
  res.second =  eVects.column(middle_index);

  return res;
}

template<typename TPoint, typename TEmbeder>
TPoint
getNormal(const TPoint &vect, const TEmbeder &emb ){
  Z3i::RealPoint p0 = emb(Z2i::RealPoint(0.0,0.0), false);
  Z3i::RealPoint px = emb(Z2i::RealPoint(10.0,0.0), false);
  Z3i::RealPoint py = emb(Z2i::RealPoint(0.0,10.0), false);
  Z3i::RealPoint u = px-p0;  u /= u.norm();
  Z3i::RealPoint v = py-p0; v /= v.norm();
  Z3i::RealPoint w = u.crossProduct(v); w /= w.norm();
  SimpleMatrix<double, 3, 3> t;
  t.setComponent(0,0, u[0]);  t.setComponent(0,1, v[0]); t.setComponent(0, 2, w[0]);
  t.setComponent(1,0, u[1]);  t.setComponent(1,1, v[1]); t.setComponent(1, 2, w[1]);
  t.setComponent(2,0, u[2]);  t.setComponent(2,1, v[2]); t.setComponent(2, 2, w[2]);
  return ((t.inverse())*vect);
}




template<typename TImage, typename TVectorField, typename TPoint>
void
fillPointArea(TImage &anImage, TVectorField &imageVectorField,
              const TPoint aPoint, const unsigned int size, const unsigned int value, const Z3i::RealPoint &n){
  typename TImage::Domain aDom(aPoint-TPoint::diagonal(size), aPoint+TPoint::diagonal(size));
  for(typename TImage::Domain::ConstIterator it= aDom.begin(); it != aDom.end(); it++){
    if (anImage.domain().isInside(*it)){
      anImage.setValue(*it, value);
      imageVectorField.setValue(*it, n);
    }
  }
}


int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, Z3i::RealPoint > VectorFieldImage3D;
  typedef ImageContainerBySTLVector < Z2i::Domain, Z3i::RealPoint > VectorFieldImage2D;
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char> Image2D;
  typedef DGtal::ConstImageAdapter<Image3D, Z2i::Domain, DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>,
                                   Image3D::Value,  DGtal::functors::Identity >  ImageAdapterExtractor;


// parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.pgm"};
  double meshScale = 1.0;
  double triangleAreaUnit = 0.01;
  double widthImageScan = {500};
  double heightImageScan = {500};
  bool orientAutoFrontX = false;
  bool orientAutoFrontY = false;
  bool orientAutoFrontZ = false;
  bool orientAuto = true;
  bool invertNormal = false;
  bool orientBack = false;
  bool exportNormals = false;
  bool setBackgroundLastDepth = false;
  bool backgroundNormalBack = false;
  int maxScan {255};
  int centerX {0};
  int centerY {0};
  int centerZ {200};


  double nx{0};
  double ny{0};
  double nz{1};
  DGtal::Z3i::RealPoint secDir; //orientation to be used when flag orientAuto is used
  unsigned int minDiagVolSize {400}; //min size used to automatically resize mesh if included in unit box.
  app.description("Convert a mesh file into a projected 2D image given from a normal direction N and from a starting point P. The 3D mesh discretized and scanned in the normal direction N, starting from P with a step 1.\n   Example:\n mesh2heightfield -i ${DGtal}/examples/samples/tref.off  heighfield.pgm  \n");
  app.add_option("-i,--input,1", inputFileName, "mesh file (.off)" )
    ->required()
    ->check(CLI::ExistingFile);
  app.add_option("-o,--output,2", outputFileName, "sequence of discrete point file (.sdp) " );
  app.add_option("--meshScale,-s", meshScale, "change the default mesh scale (each vertex multiplied by the scale) ");
  app.add_option("--remeshMinArea,-a", triangleAreaUnit, "ajust the remeshing min triangle are used to avoid empty areas");

  app.add_option("--heightFieldMaxScan", maxScan, "set the maximal scan deep." );
  app.add_option("-x,--centerX",
                              centerX, "choose x center of the projected image.");
  app.add_option("-y,--centerY",
                              centerY, "choose y center of the projected image.");
  app.add_option("-z,--centerZ",
                              centerZ, "choose z center of the projected image.");
  auto optNx = app.add_option("--nx",
                              nx, "set the x component of the projection direction.");
  auto optNy = app.add_option("--ny",
                              ny, "set the y component of the projection direction.");
  auto optNz = app.add_option("--nz",
                              nz, "set the z component of the projection direction.");
  app.add_flag("--invertNormals,-v", invertNormal, "invert normal vector of the mesh");
  app.add_option("--width", widthImageScan, "set the width of the area to be extracted as an height field image. (note that the resulting image width also depends of the scale parameter (option --meshScale))" );
  app.add_option("--height", heightImageScan, "set the height of the area to extracted  as an height field image. (note that the resulting image height also depends of the scale parameter (option --meshScale))" );
  app.add_flag("--orientAutoFrontX", orientAutoFrontX,"automatically orients the camera in front according the x axis." );
  app.add_flag("--orientAutoFrontY", orientAutoFrontY,"automatically orients the camera in front according the y axis." );
  app.add_flag("--orientAutoFrontZ", orientAutoFrontZ,"automatically orients the camera in front according the z axis." );
  app.add_flag("--orientBack", orientBack, "change the camera direction to back instead front as given in  orientAutoFront{X,Y,Z} options.");
  app.add_flag("--exportNormals", exportNormals, "export mesh normal vectors (given in the image height field basis).");
  app.add_flag("--backgroundNormalBack", backgroundNormalBack, "set the normals of background in camera opposite direction (to obtain a black background in rendering).");
  app.add_flag("--setBackgroundLastDepth", setBackgroundLastDepth, "change the default background (black with the last filled intensity).");

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  orientAuto = optNx->count() == 0 && optNy->count() == 0 && optNy->count() == 0;

  trace.info() << "Reading input file " << inputFileName ;
  Mesh<Z3i::RealPoint> inputMesh(true);
  inputMesh << inputFileName;
  std::pair<Z3i::RealPoint, Z3i::RealPoint> b = inputMesh.getBoundingBox();
  double diagDist = (b.first-b.second).norm();
  if(diagDist<minDiagVolSize){
    meshScale = minDiagVolSize/diagDist;
    inputMesh.rescale(meshScale);
  }

  triangleAreaUnit *= meshScale;
  // get vertex
  inputMesh.vertexBegin();

  inputMesh.quadToTriangularFaces();
  trace.info() << " [done] " << std::endl ;
  int maxArea = triangleAreaUnit+1.0 ;
  while (maxArea> triangleAreaUnit)
    {
      trace.info()<< "Iterating mesh subdivision ... "<< maxArea;
      maxArea = inputMesh.subDivideTriangularFaces(triangleAreaUnit);
      trace.info() << " [done]"<< std::endl;
    }

  std::pair<Z3i::RealPoint, Z3i::RealPoint> bb = inputMesh.getBoundingBox();
  Image3D::Domain meshDomain( bb.first-Z3i::Point::diagonal(1),
                              bb.second+Z3i::Point::diagonal(1));

  //  vol image filled from the mesh vertex.
  Image3D meshVolImage(meshDomain);
  VectorFieldImage3D meshNormalImage(meshDomain);
  Z3i::RealPoint z(0.0,0.0,0.0);
  for(Image3D::Domain::ConstIterator it = meshVolImage.domain().begin();
      it != meshVolImage.domain().end(); it++)
    {
      meshVolImage.setValue(*it, 0);
      meshNormalImage.setValue(*it, z);
    }

  // Filling mesh faces in volume.
  for(unsigned int i =0; i< inputMesh.nbFaces(); i++)
    {
      trace.progressBar(i, inputMesh.nbFaces());
      typename Mesh<Z3i::RealPoint>::MeshFace aFace = inputMesh.getFace(i);
      if(aFace.size()==3)
        {
          Z3i::RealPoint p1 = inputMesh.getVertex(aFace[0]);
          Z3i::RealPoint p2 = inputMesh.getVertex(aFace[1]);
          Z3i::RealPoint p3 = inputMesh.getVertex(aFace[2]);
          Z3i::RealPoint n = (p2-p1).crossProduct(p3-p1);
          n /= invertNormal? -n.norm():  n.norm();
          Z3i::RealPoint c = (p1+p2+p3)/3.0;
          fillPointArea(meshVolImage, meshNormalImage, p1, 1, 1, n);
          fillPointArea(meshVolImage, meshNormalImage, p2, 1, 1, n);
          fillPointArea(meshVolImage, meshNormalImage, p3, 1, 1, n);
          fillPointArea(meshVolImage, meshNormalImage, c, 1, 1, n);
        }
    }


  if(orientAutoFrontX || orientAutoFrontY || orientAutoFrontZ)
    {
      Z3i::Point ptL = meshVolImage.domain().lowerBound();
      Z3i::Point ptU = meshVolImage.domain().upperBound();
      Z3i::Point ptC = (ptL+ptU)/2;
      centerX=ptC[0]; centerY=ptC[1]; centerZ=ptC[2];
      nx=0; ny=0; nz=0;
    }

  if(orientAutoFrontX)
    {
      nx=(orientBack?-1.0:1.0);
      maxScan = meshVolImage.domain().upperBound()[0]- meshVolImage.domain().lowerBound()[0];
      centerX = centerX + (orientBack? maxScan/2: -maxScan/2) ;
    }
  if(orientAutoFrontY)
    {
      ny=(orientBack?-1.0:1.0);
      maxScan = meshVolImage.domain().upperBound()[1]- meshVolImage.domain().lowerBound()[1];
      centerY = centerY + (orientBack? maxScan/2: -maxScan/2);
    }
  if(orientAutoFrontZ)
    {
      nz=(orientBack?-1.0:1.0);
      maxScan = meshVolImage.domain().upperBound()[2]-meshVolImage.domain().lowerBound()[2];
      centerZ = centerZ + (orientBack? maxScan/2: -maxScan/2);
    }
  if (orientAuto)
  {
    Z3i::Point ptL = meshVolImage.domain().lowerBound();
    Z3i::Point ptU = meshVolImage.domain().upperBound();
    Z3i::Point ptC = (ptL+ptU)/2;
    int maxScanZ = meshVolImage.domain().upperBound()[2]-meshVolImage.domain().lowerBound()[2];
    int maxScanY = meshVolImage.domain().upperBound()[1]-meshVolImage.domain().lowerBound()[1];
    int maxScanX = meshVolImage.domain().upperBound()[0]-meshVolImage.domain().lowerBound()[0];
    maxScan = std::max(std::max(maxScanX, maxScanY), maxScanZ);
    trace.info() << "Automatic orientation on volume of bounds: " << ptL << " " << ptU << std::endl;
    trace.info() << "Computing main dir...";
    auto dir = getMainDirsCoVar(inputMesh.vertexBegin(), inputMesh.vertexEnd());
    trace.info() << "[done]"<<std::endl;
    nx = dir.first[0];
    ny = dir.first[1];
    nz = dir.first[2];
    secDir = dir.second;
    centerX=ptC[0]-(double)maxScan*nx/2;
    centerY=ptC[1]-(double)maxScan*ny/2;
    centerZ=ptC[2]-(double)maxScan*nz/2;


  }
  functors::Rescaling<unsigned int, Image2D::Value> scaleFctDepth(0, maxScan, 0, 255);
  if(maxScan > std::numeric_limits<Image2D::Value>::max())
    {
      trace.info()<< "Max depth value outside image intensity range: " << maxScan
                  << " (use a rescaling functor which implies a loss of precision)"  << std::endl;
    }
  trace.info() << "Processing image to output file " << outputFileName;

  Image2D::Domain aDomain2D(DGtal::Z2i::Point(0,0),
                            DGtal::Z2i::Point(widthImageScan, heightImageScan));
  Z3i::Point ptCenter (centerX, centerY, centerZ);

  Z3i::RealPoint normalDir (nx, ny, nz);
  Image2D resultingImage(aDomain2D);
  VectorFieldImage2D resultingVectorField(aDomain2D);


  for(Image2D::Domain::ConstIterator it = resultingImage.domain().begin();
      it != resultingImage.domain().end(); it++){
    resultingImage.setValue(*it, 0);
    resultingVectorField.setValue(*it, z);
  }
  DGtal::functors::Identity idV;

  unsigned int maxDepthFound = 0;
  Z3i::Point c (ptCenter-normalDir, DGtal::functors::Round<>());
  DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedder(meshVolImage.domain(),
                                                                      c,
                                                                      normalDir,
                                                                      widthImageScan);
  if (orientAuto){
   embedder =  DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >(meshVolImage.domain(),
                                                                        c,
                                                                        normalDir,
                                                                        secDir,
                                                                        widthImageScan);
  }
  int k = 0;
  bool firstFound = false;
  while (k < maxScan || !firstFound){
      embedder.shiftOriginPoint(normalDir);
      trace.progressBar(k, maxScan);
      ImageAdapterExtractor extractedImage(meshVolImage, aDomain2D, embedder, idV);
      for(Image2D::Domain::ConstIterator it = extractedImage.domain().begin();
          it != extractedImage.domain().end(); it++)
        {
          if(resultingImage(*it)== 0 &&  extractedImage(*it)!=0)
            {

              if (!firstFound) {
                firstFound = true;
              }
              maxDepthFound = k;
              resultingImage.setValue(*it, scaleFctDepth(maxScan-k));
              resultingVectorField.setValue(*it, getNormal(meshNormalImage(embedder(*it)), embedder));
            }
        }
    if (firstFound){
      k++;
    }
  }
  if (setBackgroundLastDepth)
    {
      for(Image2D::Domain::ConstIterator it = resultingImage.domain().begin();
          it != resultingImage.domain().end(); it++){
        if( resultingImage(*it)== 0 )
          {
            resultingImage.setValue(*it, scaleFctDepth(maxScan-maxDepthFound));
          }
      }
    }
  bool inverBgNormal = backgroundNormalBack;
  for(Image2D::Domain::ConstIterator it = resultingImage.domain().begin();
      it != resultingImage.domain().end(); it++){
    if(resultingVectorField(*it)==Z3i::RealPoint(0.0, 0.0, 0.0))
      {
        resultingVectorField.setValue(*it,Z3i::RealPoint(0, 0, inverBgNormal?-1: 1));
      }
  }
  resultingImage >> outputFileName;
  if(exportNormals){
    std::stringstream ss;
    ss << outputFileName << ".normals";
    std::ofstream outN;
    outN.open(ss.str().c_str(), std::ofstream::out);
    for(Image2D::Domain::ConstIterator it = resultingImage.domain().begin();
        it != resultingImage.domain().end(); it++){
      outN << (*it)[0] << " " << (*it)[1] << " " <<  0 << std::endl;
      outN <<(*it)[0]+ resultingVectorField(*it)[0] << " "
           <<(*it)[1]+resultingVectorField(*it)[1] << " "
           << resultingVectorField(*it)[2];
      outN << std::endl;
    }
  }
  return EXIT_SUCCESS;;
}
