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
 * @file heightfield2shading.cpp
 * @ingroup converters
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2015/04/11
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
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/writers/PPMWriter.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"

#include "DGtal/io/colormaps/GrayscaleColorMap.h"

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;

/**
 @page heightfield2shading heightfield2shading
 @brief Renders a 2D heightfield image into a shading image. 
 
 You can choose between lambertian model (diffuse reflectance) and
 specular model (Nayar reflectance model). You can also choose between
 a single directional light source (using -l{x,y,z} options) or use
 light source which emits in all direction (by specifying the light
 source position with -p{x,y,z} option). Another rendering mode is
 given from a bitmap reflectance map which represents the rendering for
 a normal vector value (mapped according the x/y coordinates).

@b Usage: heightfield2shading [input] [output]

@b Allowed @b options @b are:

@code

Positionals:
  1 TEXT:FILE REQUIRED                  input heightfield file (2D image).
  2 TEXT                                output image.

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         input heightfield file (2D image).
  -o,--output TEXT                      output image.
  --importNormal TEXT                   import normals from file.
  --lightDir,--lDir,--ld FLOAT x 3      light source direction: lx ly lz.
  --lightPos,--lPos,--lp FLOAT x 3      light source position: px py pz.
  -s,--specularModel FLOAT x 3          use specular Nayar model with 3 param Kdiff, Kspec, sigma.
  -r,--reflectanceMap TEXT:FILE         specify a image as reflectance map.

@endcode

@b Example:
@code
  $ heightfield2shading   ${DGtal}/examples/samples/bunnyHeightField.pgm  shading.pgm --lDir 0.0  1.0 1.0 --importNormal ${DGtal}/examples/samples/bunnyHeightField_normals.sdp -s 1.0 0.2 0.8
@endcode
You will obtain such image:
@image html  resHeightfield2shading.png "Resulting image with a 90° ccw rotation"
@see heightfield2shading.cpp

*/



template<typename TImage, typename TImageVector>
void 
computerBasicNormalsFromHeightField(const TImage &anHeightMap, TImageVector &vectorField)
{
  for(typename TImage::Domain::ConstIterator it = anHeightMap.domain().begin(); 
      it != anHeightMap.domain().end(); it++){
    if(anHeightMap.domain().isInside(*it+Z2i::Point::diagonal(1))&&
       anHeightMap.domain().isInside(*it-Z2i::Point::diagonal(1))){
      double dx = (anHeightMap(*it-Z2i::Point(1,0))-anHeightMap(*it+Z2i::Point(1,0)))/2.0;
      double dy = (anHeightMap(*it-Z2i::Point(0,1))-anHeightMap(*it+Z2i::Point(0,1)))/2.0;
      Z3i::RealPoint n (dx, dy, 1);
      n /= n.norm();
      vectorField.setValue(*it,n);
    }
  }
}



template<typename TImageVector>
void 
importNormals(std::string file, TImageVector &vectorField)
{
  std::vector<Z3i::RealPoint> vp = PointListReader<Z3i::RealPoint>::getPointsFromFile(file);
   for(unsigned int i = 0; i< vp.size(); i=i+2){
        Z3i::RealPoint p = vp.at(i);
        Z3i::RealPoint q = vp.at(i+1);
        Z3i::RealPoint n = (q-p)/(p-q).norm();
        vectorField.setValue(Z2i::Point(p[0], p[1]),n);
   }
}

// Basic Lambertian reflectance model only based on light source direction. 
template<typename TImage2D, typename TPoint3D >
struct LambertianShadindFunctor{
  LambertianShadindFunctor(const TPoint3D &aLightSourceDir ):
    myLightSourceDirection(aLightSourceDir/aLightSourceDir.norm()){}
  
  inline
  unsigned int operator()(const TPoint3D &aNormal)  const
  {
    int intensity = aNormal.dot(myLightSourceDirection)*std::numeric_limits<typename TImage2D::Value>::max();
    return intensity>0? intensity:0;
  } 
  TPoint3D  myLightSourceDirection;
};


// Basic Lambertian reflectance model from one light source positiion emeting in all direction.
template<typename TImage2D, typename TPoint3D >
struct LambertianShadindFunctorAllDirections{
  LambertianShadindFunctorAllDirections(const TPoint3D &aLightSourcePosition ):
    myLightSourcePosition(aLightSourcePosition){}
  
  inline
  unsigned int operator()(const TPoint3D &aNormal, const Z2i::Point &aPoint, const double h)  const
  {
    TPoint3D l;
    
    Z3i::RealPoint posL (aPoint[0], aPoint[1], h);
    l = -posL+myLightSourcePosition;
    l /= l.norm();
    int intensity = aNormal.dot(l)*std::numeric_limits<typename TImage2D::Value>::max();
    return intensity>0? intensity:0;
  }
  TPoint3D  myLightSourcePosition;
};





// Specular reflectance from Nayar model.
template<typename TImage2D, typename TPoint3D >
struct SpecularNayarShadindFunctor{
  SpecularNayarShadindFunctor(const TPoint3D &lightSourceDirection, const double kld,
                              const double kls, const double sigma ):
    myLightSourceDirection(lightSourceDirection/lightSourceDirection.norm()),
    myKld(kld), myKls(kls), mySigma(sigma){}
  
  inline
  unsigned int operator()(const TPoint3D &aNormal)  const {
    double lambertianIntensity = std::max(aNormal.dot(myLightSourceDirection), 0.0);
    double alpha = acos(((myLightSourceDirection+Z3i::RealPoint(0,0,1.0))/2.0).dot(aNormal/aNormal.norm()));
    double specularIntensity =  exp(-alpha*alpha/(2.0*mySigma));
    double resu = myKld*lambertianIntensity+myKls*specularIntensity;
    
    resu = std::max(resu, 0.0);
    resu = std::min(resu, 1.0);
    return resu*std::numeric_limits<typename TImage2D::Value>::max();
  }

  TPoint3D  myLightSourceDirection;
  double myKld, myKls, mySigma;
};



// Specular reflectance from Nayar model.
template<typename TImage2D, typename TPoint3D >
struct SpecularNayarShadindFunctorAllDirections{
  SpecularNayarShadindFunctorAllDirections(const TPoint3D &lightSourcePosition, const double kld,
                                           const double kls, const double sigma ):
    myLightSourcePosition(lightSourcePosition),
    myKld(kld), myKls(kls), mySigma(sigma){}
  
  inline
  unsigned int operator()(const TPoint3D &aNormal, const Z2i::Point &aPoint, const double h)  const {
    TPoint3D l;
    Z3i::RealPoint posL (aPoint[0], aPoint[1], h);
    l = -posL+myLightSourcePosition;
    l /= l.norm();
  
    double lambertianIntensity = std::max(aNormal.dot(l), 0.0);
    double alpha = acos(((l+Z3i::RealPoint(0,0,1.0))/2.0).dot(aNormal/aNormal.norm()));
    double specularIntensity =  exp(-alpha*alpha/(2.0*mySigma));
    double resu = myKld*lambertianIntensity+myKls*specularIntensity;
    
    resu = std::max(resu, 0.0);
    resu = std::min(resu, 1.0);
    return resu*std::numeric_limits<typename TImage2D::Value>::max();
  }

  TPoint3D  myLightSourcePosition;
  double myKld, myKls, mySigma;
};




// Basic Lambertian reflectance model from one light source positiion emeting in all direction.
template<typename TImage2D, typename TPoint3D >
struct ImageMapReflectance{
  ImageMapReflectance(const std::string  &filename): myImageMap (PPMReader<TImage2D>::importPPM(filename))
  {
   
    myCenterPoint = (myImageMap.domain().upperBound()-myImageMap.domain().lowerBound())/2;
    myImageRadius = min((myImageMap.domain().upperBound()-myImageMap.domain().lowerBound())[1], (myImageMap.domain().upperBound()-myImageMap.domain().lowerBound())[0])/2;
  }
  
  inline
  unsigned int operator()(const TPoint3D &aNormal)  const
  {
    Z2i::Point p(aNormal[0]*myImageRadius,aNormal[1]*myImageRadius ) ;
    p += myCenterPoint;
    if(myImageMap.domain().isInside(p)){
      return myImageMap(p);
    }else{
      return myImageMap(Z2i::Point(0,0,0));
    }
  }
  TImage2D  myImageMap;
  Z2i::Point  myCenterPoint;
  unsigned int myImageRadius;
};


struct IdColor{
  Color operator()( const unsigned int & aValue ) const{
    return DGtal::Color(aValue);
  }
};



int main( int argc, char** argv )
{
  //  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned int> Image2D;
  typedef ImageSelector < Z2i::Domain, unsigned char>::Type Image2D;
  //typedef ImageContainerBySTLVector < Z2i::Domain, char> Image2DChar;
  typedef ImageContainerBySTLVector < Z2i::Domain, Z3i::RealPoint> Image2DNormals;



  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.pgm"};
  std::string normalFileName {""};
  double lx, ly, lz, px, py, pz;
  bool usingAllDirectionLightSource = false;
  std::vector<double> specularModel;
  std::string reflectanceMap;
  std::vector<double> lDir;
  std::vector<double> lPos;
  
  
  app.description("Render a 2D heightfield image into a shading image. You can choose between lambertian model (diffuse reflectance) and specular model (Nayar reflectance model). You can also choose between a single directional light source (using --lightDir option) or use light source which emits in all direction (by specifying the light source position with --lightPos} option). Another rendering mode is given from a bitmap reflectance map which represents the rendering for a normal vector value (mapped according the x/y coordinates).\nExample:\n   ${DGtal}/examples/samples/bunnyHeightField.pgm  shading.pgm --lPos 10.0  -120.0 550.0 --importNormal ${DGtal}/examples/samples/bunnyHeightField_normals.sdp -s 1.0 0.2 0.8 \n");
  app.add_option("-i,--input,1", inputFileName, "input heightfield file (2D image).")
     ->check(CLI::ExistingFile)
     ->required();
  app.add_option("-o,--output,2", outputFileName,"output image.");
  app.add_option("--importNormal", normalFileName, "import normals from file.");
  
  app.add_option("--lightDir,--lDir,--ld", lDir, "light source direction: lx ly lz.")
    ->expected(3);
  app.add_option("--lightPos,--lPos,--lp", lPos, "light source position: px py pz.")
    ->expected(3);
  app.add_option("-s,--specularModel", specularModel, "use specular Nayar model with 3 param Kdiff, Kspec, sigma." )
    ->expected(3);
  app.add_option("-r,--reflectanceMap",reflectanceMap, "specify a image as reflectance map.")
    ->check(CLI::ExistingFile);
   
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------


  
  if(lDir.size() == 3)
    {
      lx = lDir[0];
      ly = lDir[1];
      lz = lDir[2];
    }
  else if(lPos.size() == 3)
    {
      px = lPos[0];
      py = lPos[1];
      pz = lPos[2];  
      usingAllDirectionLightSource = true;
    }
  else if (reflectanceMap == "")
    {
      trace.error() << "You need to specify either the light source direction or position (if you use a all directions model)." << std::endl;
      exit(0);
    }

  LambertianShadindFunctor<Image2D, Z3i::RealPoint> lShade (Z3i::RealPoint(lx,ly,lz));
  LambertianShadindFunctorAllDirections<Image2D, Z3i::RealPoint> lShadePosD (Z3i::RealPoint(px ,py, pz));
  SpecularNayarShadindFunctor<Image2D, Z3i::RealPoint> lSpecular (Z3i::RealPoint(lx,ly,lz), 0, 0, 0);  
  SpecularNayarShadindFunctorAllDirections<Image2D, Z3i::RealPoint> lSpecularPosD (Z3i::RealPoint(px,py,pz), 0, 0, 0);  

  
  bool useSpecular = false;
  if(specularModel.size() == 3){
    useSpecular = true;
    lSpecular.myKld = specularModel[0];
    lSpecular.myKls = specularModel[1];
    lSpecular.mySigma = specularModel[2];
    lSpecularPosD.myKld = specularModel[0];
    lSpecularPosD.myKls = specularModel[1];
    lSpecularPosD.mySigma = specularModel[2];
    if(specularModel[2]==0.0)      
    {
      trace.error()<< "a 0 value for sigma is not possible in the Nayar model, please change it. "<< std::endl;
      exit(1);
    }

  }

  
  trace.info() << "Reading input file " << inputFileName ; 
  Image2D inputImage = DGtal::GenericReader<Image2D>::import(inputFileName);  
  Image2DNormals vectNormals (inputImage.domain());
  Image2D result (inputImage.domain());
  if(normalFileName != ""){
    importNormals(normalFileName, vectNormals);
  }else{
    computerBasicNormalsFromHeightField(inputImage, vectNormals);
  }
  if(reflectanceMap != "")
    {
      ImageMapReflectance<Image2D, Z3i::RealPoint> lMap(reflectanceMap);
      for(typename Image2D::Domain::ConstIterator it = inputImage.domain().begin(); 
          it != inputImage.domain().end(); it++){
        if(reflectanceMap != "")
          {
            result.setValue(*it, lMap(vectNormals(*it)));           
          }       
      }        
      IdColor id;
      PPMWriter<Image2D, IdColor  >::exportPPM(outputFileName, result, id);
    }

  else
    {
      for(typename Image2D::Domain::ConstIterator it = inputImage.domain().begin(); 
          it != inputImage.domain().end(); it++){
        if(usingAllDirectionLightSource)
          {
            result.setValue(*it, useSpecular? lSpecularPosD(vectNormals(*it), *it, inputImage(*it)):
                            lShadePosD(vectNormals(*it), *it, inputImage(*it))); 
          }
        else
          {
            result.setValue(*it, useSpecular? lSpecular(vectNormals(*it)):lShade(vectNormals(*it))); 
          }
        
      }
      result >> outputFileName;
   }
  
  return 0;  
}
