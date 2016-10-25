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
 * @file 3dHeightMapViewer.cpp
 * @ingroup Examples
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2014/06/2014
 *
 * An example file named qglViewer.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <climits>

#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include <DGtal/io/colormaps/GradientColorMap.h>
#include <DGtal/io/colormaps/GrayscaleColorMap.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include "DGtal/io/Color.h"

#include <DGtal/base/ConstAlias.h>
#include <DGtal/kernel/BasicPointFunctors.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/io/colormaps/GradientColorMap.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;
using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;


/**
 @page Doc3dHeightMapViewer 3dHeightMapViewer
 
 @brief  Displays 2D image as heightmap by using QGLviewer.

 @b Usage:  3dImageViewer [options] input
 

 @b Allowed @b options @b are :
 
 @code
  -h [ --help ]                  display this message
  -s [ --scale ] arg (=1)        set the scale of the maximal level. (default 
                                 1.0)
  -c [ --colorMap ]              define the heightmap color with a pre-defined 
                                 colormap (GradientColorMap)
  -t [ --colorTextureImage ] arg define the heightmap color from a given color 
                                 image (32 bits image).
  -i [ --input-file ] arg        2d input image representing the height map 
                                 (given as grayscape image cast into 8 bits).

 @endcode


 @b Example: 
 
 @code
 $ visualisation/3dHeightMapViewer -i ${DGtal}/examples/samples/church.pgm -s 0.2
 @endcode 

 You should obtain such a result:

 @image html res3dHeightMapViewer.png "resulting heightmap visualisation."
 

 @see
 @ref 3dHeightMapViewer.cpp

 */

template < typename Space = DGtal::Z3i::Space, typename KSpace = DGtal::Z3i::KSpace>
struct Viewer3DImageSpec: public DGtal::Viewer3D <Space, KSpace>
{
  Viewer3DImageSpec(Z2i::RealPoint pup, Z2i::RealPoint plow): myZScale(1.0), myZposClipping(10.0),
                                                              myPUpper(pup), myPLower(plow){
    DGtal::Viewer3D<>::update();
  };

virtual  void
init(){
   DGtal::Viewer3D<>::init();
   QGLViewer::setKeyDescription ( Qt::Key_Up, "Move Up the cutting plane in the Z axis direction." );
   QGLViewer::setKeyDescription ( Qt::Key_Down, "Move Down the cutting plane in the Z axis direction." );
   QGLViewer::setKeyDescription ( Qt::Key_Shift, "Change the cutting plane move with step 1 (5 by default)" );
   QGLViewer::setKeyDescription ( Qt::Key_Plus, "Increase the openGL scale by 0.25 " );
   QGLViewer::setKeyDescription ( Qt::Key_Minus, "Decrease the openGL scale by 0.25 " );
};


protected:
  double myZScale;
  double myZposClipping;
  Z2i::RealPoint myPUpper, myPLower;

  virtual void draw (  ){
    DGtal::Viewer3D<>::draw();
    glPushMatrix();
    glMultMatrixd ( DGtal::Viewer3D<>::manipulatedFrame()->matrix() );
    glPushMatrix();
    glScalef(1.0, 1.0, myZScale);
    glEnable( GL_LIGHTING );
    glBegin( GL_QUADS );
    glColor4ub( 50, 50, 240,  150 );
    glNormal3f( 0, 0 , 1.0 );
    glVertex3f( myPLower[0], myPLower[1] , myZposClipping );
    glVertex3f( myPUpper[0], myPLower[1] , myZposClipping );
    glVertex3f( myPUpper[0], myPUpper[1], myZposClipping );
    glVertex3f( myPLower[0], myPUpper[1],  myZposClipping );
    glEnd();
    glPopMatrix();
    glPopMatrix();
  }
  virtual void keyPressEvent ( QKeyEvent *e ){
    bool handled = false;
    if( e->key() == Qt::Key_Up){
      if((e->modifiers() & Qt::MetaModifier)){
        myZposClipping+=1;
      }else{
        myZposClipping+=5;
      }

      DGtal::Viewer3D<>::update();
      handled=true;
    }
    if( e->key() == Qt::Key_Down){
      if((e->modifiers() & Qt::MetaModifier)){
        myZposClipping-=1;
      }else{
        myZposClipping-=5;
      }

      DGtal::Viewer3D<>::update();
      handled=true;
    }

    if( e->key() == Qt::Key_Plus){
      myZScale+=0.25;
      DGtal::Viewer3D<Space, KSpace>::setGLScale(1.0, 1.0, myZScale);
      DGtal::Viewer3D<Space, KSpace>::update();
      handled=true;
    }
    if( e->key() == Qt::Key_Minus){
      myZScale-=0.25;
      DGtal::Viewer3D<Space, KSpace>::setGLScale(1.0, 1.0, myZScale);
      DGtal::Viewer3D<>::update();
      handled=true;
    }
    if ( !handled )
      DGtal::Viewer3D<>::keyPressEvent ( e );
  }
};




// Defining a Helper to get the 3D point functor from an 2DImage
template<typename TImage2D, typename TPoint3D >
struct Image3DPredicatFrom2DImage{
  typedef  TPoint3D Point3D;
  /**
   *  Construct the predicat given a 2D Image
   **/
  Image3DPredicatFrom2DImage(DGtal::ConstAlias<TImage2D> anImage, double aScale):myImageRef(anImage),
                                                                                 myScale(aScale){
  }
  inline
  bool operator()(const Point3D &aPoint)  const {
    functors::Projector<SpaceND<2, typename TImage2D::Integer> > projXY;
    return  (*myImageRef)(projXY(aPoint))*myScale >= aPoint[2];
  }
  CountedConstPtrOrConstPtr<TImage2D> myImageRef;
  double myScale;
};








int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("scale,s", po::value<double>()->default_value(1.0), "set the scale of the maximal level. (default 1.0)")
    ("colorMap,c", "define the heightmap color with a pre-defined colormap (GradientColorMap)")
    ("colorTextureImage,t",po::value<std::string>(),  "define the heightmap color from a given color image (32 bits image).")
    ("input-file,i", po::value<std::string>(), "2d input image representing the height map (given as grayscape image cast into 8 bits)." );

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
      std::cout << "Usage: " << argv[0] << " [input-file]\n"
                << "Displays 2D image as heightmap by using QGLviewer. "<< endl
                << general_opt << "\n" <<
                "Exemple of use:  visualisation/3dHeightMapViewer -i ${DGtal}/examples/samples/church.pgm -s 0.2" << std::endl;

      return 0;
    }

  if(! vm.count("input-file"))
    {
      trace.error() << " The file name was defined" << endl;
      return 0;
    }


  string inputFilename = vm["input-file"].as<std::string>();
  double scale = vm["scale"].as<double>();


  typedef DGtal::ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image2DG ;
  typedef DGtal::ImageContainerBySTLVector<Z2i::Domain, unsigned int> Image2DCol ;

  Image2DG image = GenericReader<Image2DG>::import( inputFilename );
  Image2DCol imageTexture(image.domain());
  int maxHeight = (int)(std::numeric_limits<Image2DG::Value>::max()*scale);
  trace.info()<< "Max height from scale:" << maxHeight << std::endl;

  if(vm.count("colorTextureImage")){
    imageTexture =  GenericReader<Image2DCol>::import( vm["colorTextureImage"].as<std::string>() );
  }

  QApplication application(argc,argv);
  Z2i::RealPoint plow (image.domain().lowerBound()[0]-0.5,
                       image.domain().lowerBound()[1]-0.5);
  
  Z2i::RealPoint pup (image.domain().upperBound()[0]+0.5,
                      image.domain().upperBound()[1]+0.5);
  
  Viewer3DImageSpec<> viewer(plow, pup) ;
  viewer.setWindowTitle("Height Map Viewer");
  viewer.show();


  KSpace K;
  K.init(Z3i::Point(0,0,0),Z3i::Point(image.domain().upperBound()[0], image.domain().upperBound()[1], maxHeight+1), true);
  std::set<KSpace::SCell> boundVect;
  Image3DPredicatFrom2DImage<Image2DG, Z3i::Point> image3Dpredicate(image, scale);
  trace.info() << "Constructing boundary... ";
  Surfaces<KSpace>::sMakeBoundary (boundVect, K, image3Dpredicate, Z3i::Point(0,0,0),
                                   Z3i::Point(image.domain().upperBound()[0], image.domain().upperBound()[1], maxHeight+1)  );
  trace.info() << "[done]"<< std::endl;

  viewer << SetMode3D((*(boundVect.begin())).className(), "Basic" );
  GradientColorMap<Image2DG::Value,CMAP_JET>  gradientShade( 0, std::numeric_limits<Image2DG::Value>::max());
  GrayscaleColorMap<Image2DG::Value>  grayShade(0, std::numeric_limits<Image2DG::Value>::max());


  for(std::set<KSpace::SCell>::const_iterator it = boundVect.begin();
      it!= boundVect.end(); it++){
    Z3i::Point pt = K.sCoords(K.sDirectIncident( *it, 2 ));
    functors::Projector<SpaceND<2,int> > proj;
    Image2DG::Value val = image(proj(pt));
    if(vm.count("colorMap")){
      viewer.setFillColor(gradientShade(val));
    }else if (vm.count("colorTextureImage")) {
      viewer.setFillColor(Color(imageTexture(proj(pt))));
    }else{
      viewer.setFillColor(grayShade(val));
    }
  viewer << *it;
  }

  viewer << Viewer3D<>::updateDisplay;
  return application.exec();
}
