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
 * @ingroup converters
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
#include "DGtal/math/linalg/SimpleMatrix.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;



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
  
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "mesh file (.off) " )
    ("meshScale,s", po::value<double>()->default_value(1.0),
     "change the default mesh scale (each vertex multiplied by the scale) " )
    ("output,o", po::value<std::string>(), "sequence of discrete point file (.sdp) ") 
    ("orientAutoFrontX",  "automatically orients the camera in front according the x axis." )
    ("orientAutoFrontY",  "automatically orients the camera in front according the y axis." )
    ("orientAutoFrontZ",  "automatically orients the camera in front according the z axis." )
    ("orientBack",  "change the camera direction to back instead front as given in  orientAutoFront{X,Y,Z} options." )
    ("nx", po::value<double>()->default_value(0), "set the x component of the projection direction." )
    ("ny", po::value<double>()->default_value(0), "set the y component of the projection direction." )
    ("nz", po::value<double>()->default_value(1), "set the z component of the projection direction." )
    ("centerX,x", po::value<unsigned int>()->default_value(0), "choose x center of the projected image." )
    ("centerY,y", po::value<unsigned int>()->default_value(0), "choose y center of the projected image." )
    ("centerZ,z", po::value<unsigned int>()->default_value(1), "choose z center of the projected image." )
    ("width", po::value<unsigned int>()->default_value(100), "set the width of the area to be extracted as an height field image." )
    ("height", po::value<unsigned int>()->default_value(100), "set the height of the area to extracted  as an height field image." )
    ("heightFieldMaxScan", po::value<unsigned int>()->default_value(255), "set the maximal scan deep." )
    ("exportNormals",  "export mesh normal vectors (given in the image height field basis)." )
    ("backgroundNormalBack",  "set the normals of background in camera opposite direction (to obtain a black background in rendering). " )
  
    ("setBackgroundLastDepth", "change the default background (black with the last filled intensity).");
  
  
  double triangleAreaUnit = 0.5;
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
      std::cout << "Usage: " << argv[0] << " [input] [output]\n"
		<< "Convert a mesh file into a projected 2D image given from a normal direction N and from a starting point P. The 3D mesh discretized and scanned in the normal direction N, starting from P with a step 1."
		<< general_opt << "\n";
      std::cout << "Example:\n"
		<< "mesh2heightfield -i ${DGtal}/examples/samples/tref.off --orientAutoFrontZ --width 25 --height 25 -o heighMap.pgm -s 10  \n";
      return 0;
    }
  
  if(! vm.count("input") ||! vm.count("output"))
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;      
      return 0;
    }
  
  string inputFilename = vm["input"].as<std::string>();
  string outputFilename = vm["output"].as<std::string>();
  double meshScale = vm["meshScale"].as<double>();
  trace.info() << "Reading input file " << inputFilename ; 
  Mesh<Z3i::RealPoint> inputMesh(true);
      
  DGtal::MeshReader<Z3i::RealPoint>::importOFFFile(inputFilename, inputMesh);
  std::pair<Z3i::RealPoint, Z3i::RealPoint> b = inputMesh.getBoundingBox();
  double diagDist = (b.first-b.second).norm();
  if(diagDist<2.0*sqrt(2.0)){
    inputMesh.changeScale(2.0*sqrt(2.0)/diagDist);  
  }
 
  inputMesh.quadToTriangularFaces();
  inputMesh.changeScale(meshScale);  
  trace.info() << " [done] " << std::endl ; 
  double maxArea = triangleAreaUnit+1.0 ;
  while(maxArea> triangleAreaUnit)
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
          n /= n.norm();
          Z3i::RealPoint c = (p1+p2+p3)/3.0;
          fillPointArea(meshVolImage, meshNormalImage, p1, 1, 1, n);          
          fillPointArea(meshVolImage, meshNormalImage, p2, 1, 1, n);          
          fillPointArea(meshVolImage, meshNormalImage, p3, 1, 1, n);          
          fillPointArea(meshVolImage, meshNormalImage, c, 1, 1, n);          
        }
    }
    
  unsigned int widthImageScan = vm["height"].as<unsigned int>()*meshScale;
  unsigned int heightImageScan = vm["width"].as<unsigned int>()*meshScale;

  int maxScan = vm["heightFieldMaxScan"].as<unsigned int>()*meshScale;  
  int centerX = vm["centerX"].as<unsigned int>()*meshScale;
  int centerY = vm["centerY"].as<unsigned int>()*meshScale;
  int centerZ = vm["centerZ"].as<unsigned int>()*meshScale;
  
  double nx = vm["nx"].as<double>();
  double ny = vm["ny"].as<double>();
  double nz = vm["nz"].as<double>();

  if(vm.count("orientAutoFrontX")|| vm.count("orientAutoFrontY") ||
     vm.count("orientAutoFrontZ"))
    {
      Z3i::Point ptL = meshVolImage.domain().lowerBound();
      Z3i::Point ptU = meshVolImage.domain().upperBound();     
      Z3i::Point ptC = (ptL+ptU)/2;
      centerX=ptC[0]; centerY=ptC[1]; centerZ=ptC[2];
      nx=0; ny=0; nz=0;
    }  
  bool opp = vm.count("orientBack");
  if(vm.count("orientAutoFrontX"))
    {
      nx=(opp?-1.0:1.0); 
      maxScan = meshVolImage.domain().upperBound()[0]- meshVolImage.domain().lowerBound()[0];
      centerX = centerX + (opp? maxScan/2: -maxScan/2) ;
    }
  if(vm.count("orientAutoFrontY"))
    {
      ny=(opp?-1.0:1.0); 
      maxScan = meshVolImage.domain().upperBound()[1]- meshVolImage.domain().lowerBound()[1];
      centerY = centerY + (opp? maxScan/2: -maxScan/2);
    }
  if(vm.count("orientAutoFrontZ"))
    {
      nz=(opp?-1.0:1.0); 
      maxScan = meshVolImage.domain().upperBound()[2]-meshVolImage.domain().lowerBound()[2];
      centerZ = centerZ + (opp? maxScan/2: -maxScan/2);
    }
  
  functors::Rescaling<unsigned int, Image2D::Value> scaleFctDepth(0, maxScan, 0, 255);  
  if(maxScan > std::numeric_limits<Image2D::Value>::max())
    {
      trace.info()<< "Max depth value outside image intensity range: " << maxScan 
                  << " (use a rescaling functor which implies a loss of precision)"  << std::endl; 
    }  
  trace.info() << "Processing image to output file " << outputFilename; 
    
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
  for(unsigned int k=0; k < maxScan; k++)
    {
      trace.progressBar(k, maxScan);
      DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedder(meshVolImage.domain(), 
                                                                          ptCenter+normalDir*k,
                                                                          normalDir,
                                                                          widthImageScan);
      ImageAdapterExtractor extractedImage(meshVolImage, aDomain2D, embedder, idV);
      for(Image2D::Domain::ConstIterator it = extractedImage.domain().begin(); 
          it != extractedImage.domain().end(); it++)
        {
          if(resultingImage(*it)== 0 &&  extractedImage(*it)!=0)
            {
              maxDepthFound = k;
              resultingImage.setValue(*it, scaleFctDepth(maxScan-k));
              resultingVectorField.setValue(*it, getNormal(meshNormalImage(embedder(*it)), embedder));
            }
        }    
    }
  if (vm.count("setBackgroundLastDepth"))
    {
      for(Image2D::Domain::ConstIterator it = resultingImage.domain().begin(); 
          it != resultingImage.domain().end(); it++){
        if( resultingImage(*it)== 0 )
          {
            resultingImage.setValue(*it, scaleFctDepth(maxScan-maxDepthFound));
          }
      }
    } 
  bool inverBgNormal = vm.count("backgroundNormalBack");
  for(Image2D::Domain::ConstIterator it = resultingImage.domain().begin(); 
      it != resultingImage.domain().end(); it++){
    if(resultingVectorField(*it)==Z3i::RealPoint(0.0, 0.0, 0.0))      
      {
        resultingVectorField.setValue(*it,Z3i::RealPoint(0, 0, inverBgNormal?-1: 1));         
      }
  }
  resultingImage >> outputFilename;
  if(vm.count("exportNormals")){ 
    std::stringstream ss;
    ss << outputFilename << ".normals";
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

  return 0;  
}

