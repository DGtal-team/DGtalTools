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
 * @ingroup Converters
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
 @ingroup convertertools


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
   1 TEXT:FILE REQUIRED                  mesh file (.off)
   2 TEXT=result.pgm                     sequence of discrete point file (.sdp)

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         mesh file (.off)
   -o,--output TEXT=result.pgm           sequence of discrete point file (.sdp)
   -s,--meshScale FLOAT                  change the default mesh scale (each vertex multiplied by the scale)
   -a,--remeshMinArea FLOAT=0.01         ajust the remeshing min triangle are used to avoid empty areas
   --heightFieldMaxScan INT=255          set the maximal scan deep.
   -x,--centerX INT=0                    choose x center of the projected image.
   -y,--centerY INT=0                    choose y center of the projected image.
   -z,--centerZ INT=200                  choose z center of the projected image.
   --nx FLOAT=0                          set the x component of the projection direction.
   --ny FLOAT=0                          set the y component of the projection direction.
   --nz FLOAT=1                          set the z component of the projection direction.
   -v,--invertNormals                    invert normal vector of the mesh
   --width FLOAT=500                     set the width of the area to be extracted as an height field image. (note that the resulting image width also depends of the scale parameter (option --meshScale))
   --height FLOAT=500                    set the height of the area to extracted  as an height field image. (note that the resulting image height also depends of the scale parameter (option --meshScale))
   --orientAutoFrontX                    automatically orients the camera in front according the x axis.
   --orientAutoFrontY                    automatically orients the camera in front according the y axis.
   --orientAutoFrontZ                    automatically orients the camera in front according the z axis.
   --orientBack                          change the camera direction to back instead front as given in  orientAutoFront{X,Y,Z} options.
   --exportNormals                       export mesh normal vectors (given in the image height field basis).
   --backgroundNormalBack                set the normals of background in camera opposite direction (to obtain a black background in rendering).
   --setBackgroundLastDepth              change the default background (black with the last filled intensity).


@endcode

@b Example:
@code
  $ heightfield2shading   ${DGtal}/examples/samples/bunnyHeightField.pgm  shading.pgm --lDir 0.0  1.0 1.0 --importNormal ${DGtal}/examples/samples/bunnyHeightField_normals.sdp -s 1.0 0.2 0.8
@endcode
You will obtain such image:
@image html  resHeightfield2shading.png "Resulting image with a 90° ccw rotation"

@b Other example:
@code
  $ heightfield2shading   ${DGtal}/examples/samples/bunnyHeightField.pgm  shading.ppm  --importNormal ${DGtal}/examples/samples/bunnyHeightField_normals.sdp --hsvShading
@endcode
You will obtain such image:
@image html  resHueHeightfield2shading.png "Resulting image with a 90° ccw rotation (and conversion to png)"
@see heightfield2shading.cpp

*/



template<typename TImage, typename TImageVector>
void 
computerBasicNormalsFromHeightField(const TImage &anHeightMap, TImageVector &vectorField,
                                    bool invertN = false)
{
  for(typename TImage::Domain::ConstIterator it = anHeightMap.domain().begin(); 
      it != anHeightMap.domain().end(); it++){
    if(anHeightMap.domain().isInside(*it+Z2i::Point::diagonal(1))&&
       anHeightMap.domain().isInside(*it-Z2i::Point::diagonal(1))){
      double dx = (anHeightMap(*it-Z2i::Point(1,0))-anHeightMap(*it+Z2i::Point(1,0)))/2.0;
      double dy = (anHeightMap(*it-Z2i::Point(0,1))-anHeightMap(*it+Z2i::Point(0,1)))/2.0;
      Z3i::RealPoint n (dx, dy, 1);
      n /= n.norm();
      vectorField.setValue(*it,invertN? -n : n);
    }
  }
}



template<typename TImageVector>
void 
importNormals(std::string file, TImageVector &vectorField,
              bool invertN = false)
{
  std::vector<Z3i::RealPoint> vp = PointListReader<Z3i::RealPoint>::getPointsFromFile(file);
  trace.info() << "import done: " << vp.size() <<  std::endl;
  for(unsigned int i = 0; i< vp.size()-1; i=i+2){
        Z3i::RealPoint p = vp.at(i);
        Z3i::RealPoint q = vp.at(i+1);
        Z3i::RealPoint n = (q-p)/(p-q).norm();
        vectorField.setValue(Z2i::Point(p[0], p[1]),invertN? -n : n);
  }
  trace.info() <<endl;
}


template<typename TImageVector>
void 
importNormalsOrdDir(std::string file, TImageVector &vectorField,
                    unsigned int width, unsigned int height, bool invertN = false)
{
  std::vector<Z3i::RealPoint> vp = PointListReader<Z3i::RealPoint>::getPointsFromFile(file);

  for(unsigned int i = 0; i< vp.size()-1; i++){
    Z2i::Point p ( i-(width*floor((i/width))), i/width);
    Z3i::RealPoint n = vp[i];
    vectorField.setValue(p,invertN? -n : n);
  }
  trace.info() <<endl;
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



DGtal::Color
colorFromHSB(double h, double saturation, double value){
  double r, g, b;
  DGtal::Color::HSVtoRGB(r, g, b, h,saturation, value);
  return DGtal::Color(r*255.0,g*255.0,b*255.0);
  
}


int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned int> Image2DC;
  typedef ImageSelector < Z2i::Domain, unsigned char>::Type Image2D;
  typedef ImageContainerBySTLVector < Z2i::Domain, Z3i::RealPoint> Image2DNormals;



  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.pgm"};
  std::string normalFileName {""};
  double lx, ly, lz, px, py, pz;
  bool usingAllDirectionLightSource = false;
  bool useOrderedImportNormal = false;
  bool hsvShading = false;
  bool normalMap = false;
  bool invertNormals = false;
  std::vector<double> specularModel;
  std::string reflectanceMap;
  std::vector<double> lDir = {0, 0, 1};
  std::vector<double> lPos;
  std::vector<unsigned int> domain;  
  
  app.description("Render a 2D heightfield image into a shading image. You can choose between lambertian model (diffuse reflectance) and specular model (Nayar reflectance model). You can also choose between a single directional light source (using --lightDir option) or use light source which emits in all direction (by specifying the light source position with --lightPos} option). Another rendering mode is given from a bitmap reflectance map which represents the rendering for a normal vector value (mapped according the x/y coordinates).\nExample:\n heightfield2shading   ${DGtal}/examples/samples/bunnyHeightField.pgm  shading.pgm --lPos 10.0  -120.0 550.0 --importNormal ${DGtal}/examples/samples/bunnyHeightField_normals.sdp -s 1.0 0.2 0.8 \n"
    "Other example: heightfield2shading   ${DGtal}/examples/samples/bunnyHeightField.pgm  shading.ppm  --importNormal ${DGtal}/examples/samples/bunnyHeightField_normals.sdp --hsvShading\n");
  auto opt1 = app.add_option("-i,--input,1", inputFileName, "input heightfield file (2D image).")
    ->check(CLI::ExistingFile);
  auto domOpt =  app.add_option("--domain,-d", domain , "specify the domain (required when normal are imported and if --inout is not given).")
    -> expected(2);
  app.add_option("-o,--output,2", outputFileName,"output image.");
  auto impNOpt = app.add_option("--importNormal", normalFileName, "import normals from file.");
  app.add_flag("--orderedNormalsImport",useOrderedImportNormal, "Use ordered normals." );
  app.add_option("--lightDir,--lDir,--ld", lDir, "light source direction: lx ly lz.")
    ->expected(3);
  app.add_option("--lightPos,--lPos,--lp", lPos, "light source position: px py pz.")
    ->expected(3);
  app.add_option("-s,--specularModel", specularModel, "use specular Nayar model with 3 param Kdiff, Kspec, sigma." )
    ->expected(3);
  app.add_option("-r,--reflectanceMap",reflectanceMap, "specify a image as reflectance map.")
    ->check(CLI::ExistingFile);
  app.add_flag("--hsvShading", hsvShading, "use shading with HSV shading (given from the normal vector)");
  app.add_flag("--normalMap", normalMap, "generates normal map.");
  app.add_flag("--invertNormals,-v", invertNormals, "invert normal orientations.");
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  if(! *opt1 && !(*domOpt && *impNOpt ) ){
    trace.error() << "You need either set input file (--input) or use a domain (--domain) with the --importNormal option." << std::endl;
    exit(0);
  }

  
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
  else if (reflectanceMap == "" && ! hsvShading && !normalMap)
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

  
  Image2D inputImage(Z2i::Domain(Z2i::Point(0,0), Z2i::Point(0,0) ));
  if (inputFileName != "") {
    trace.info() << "Reading input file " << inputFileName ; 
    inputImage = DGtal::GenericReader<Image2D>::import(inputFileName);
    trace.info() << "[done]" << std::endl;
  }
  else {
    inputImage = Image2D(Z2i::Domain(Z2i::Domain(Z2i::Point(0,0),
                                                 Z2i::Point(domain[0],domain[1]) )));
  }
  
  Image2DNormals vectNormals (inputImage.domain());
  Image2D result (inputImage.domain());
  Image2DC resultC (inputImage.domain());

  if(normalFileName != ""){
    trace.info() << "Import normal file " << inputFileName << vectNormals.domain();
    if (useOrderedImportNormal)
    {
      importNormalsOrdDir(normalFileName, vectNormals,
                          inputImage.domain().upperBound()[0]+1,
                          inputImage.domain().upperBound()[1]+1, invertNormals);
    }
    else
    {
      importNormals(normalFileName, vectNormals, invertNormals);
    }
    trace.info() << "[done]" << std::endl;
  }
  else
  {
    computerBasicNormalsFromHeightField(inputImage, vectNormals, invertNormals);
  }
  if (hsvShading)
  {
    for(typename Image2D::Domain::ConstIterator it = inputImage.domain().begin();
        it != inputImage.domain().end(); it++){
      auto n = vectNormals(*it);
      double sat = 1.0*( sin(acos(Z3i::RealPoint(0.0,0.0,1.0).dot(n))));
      double value = 1.0;
      double hue = ((int)(((2.0*M_PI+atan2(Z3i::RealPoint(0.0,1.0,0.0).dot(n),
                     Z3i::RealPoint(1.0,0.0,0.0).dot(n)))/(2.0*M_PI))*360.0+100))%360;
      DGtal::uint32_t colCode = colorFromHSB(hue, sat, value).getRGB();
      resultC.setValue(*it, colCode);
    }
    IdColor id;
    PPMWriter<Image2DC, IdColor  >::exportPPM(outputFileName, resultC, id);
  }
  else if (normalMap)
  {
    DGtal::functors::Rescaling<double, unsigned int> rgRescale (-1.0, 1.0, 0, 255);
    DGtal::functors::Rescaling<double, unsigned int> bRescale (0.0, -1.0, 128, 255);
    for(typename Image2D::Domain::ConstIterator it = inputImage.domain().begin();
        it != inputImage.domain().end(); it++){
        auto n = vectNormals(*it);
        DGtal::Color c (rgRescale(n[0]), rgRescale(n[1]), bRescale(n[2]) );
        resultC.setValue(*it, c.getRGB());
    }
    IdColor id;
    PPMWriter<Image2DC, IdColor  >::exportPPM(outputFileName, resultC, id);
  }
  else if(reflectanceMap != "")
    {
      ImageMapReflectance<Image2D, Z3i::RealPoint> lMap(reflectanceMap);
      for(typename Image2D::Domain::ConstIterator it = inputImage.domain().begin(); 
          it != inputImage.domain().end(); it++){
            result.setValue(*it, lMap(vectNormals(*it)));
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
