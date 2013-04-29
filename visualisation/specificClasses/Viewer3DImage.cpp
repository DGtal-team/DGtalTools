
/**
 * @file Viewer3DImage.cpp
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2011/01/03
 *
 * Implementation of methods defined in Viewer3DImage.h
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////

#include "Viewer3DImage.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/writers/PGMWriter.h"
#include "DGtal/io/readers/PNMReader.h"

#include <QKeyEvent>

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace qglviewer;


///////////////////////////////////////////////////////////////////////////////
// class Viewer3DImage
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :



void
Viewer3DImage::init(){
   Viewer3D::init();
   setKeyDescription ( Qt::Key_X, "Change the current axis to X for the current 2D image slice setting." );
   setKeyDescription ( Qt::Key_Y, "Change the current axis to Y for the current 2D image slice setting." );
   setKeyDescription ( Qt::Key_Z, "Change the current axis to Z for the current 2D image slice setting." );
   setKeyDescription ( Qt::Key_Up, "Move the current 2D image slice to 5 in the positive direction of the current axis." );
   setKeyDescription ( Qt::Key_Down, "Move the current 2D image slice to 5 in the negative direction of the current axis." );
   setKeyDescription ( Qt::Key_Shift, "Change the slice move with step 1 (5 by default)" );


}

void 
Viewer3DImage::setVolImage(Image3D * an3DImage){
  my3dImage = an3DImage;

  // Adding X slice in the viewer.
  Image2D sliceImageX = DGtal::extractLowerDimImage<Image3D, Image2D>(*my3dImage, 0, mySliceXPos);
  (*this) << sliceImageX;
  (*this) << DGtal::UpdateImagePosition(0, DGtal::Display3D::xDirection, mySliceXPos, 0.0, 0.0);


  // Adding Y slice in the viewer.
  Image2D sliceImageY = DGtal::extractLowerDimImage<Image3D, Image2D>(*my3dImage, 1, mySliceYPos);
  (*this) << sliceImageY; 
  (*this) << DGtal::UpdateImagePosition(1, DGtal::Display3D::yDirection, 0.0, mySliceYPos, 0.0);

  // Adding Y slice in the viewer.
  Image2D sliceImageZ = DGtal::extractLowerDimImage<Image3D, Image2D>(*my3dImage, 2, mySliceZPos);
  (*this) << sliceImageZ;
  (*this) << DGtal::UpdateImagePosition(2, DGtal::Display3D::zDirection, 0.0, 0.0, mySliceZPos);

    
  (*this) << Viewer3D::updateDisplay;
}




QString
Viewer3DImage::helpString() const
{
  QString text ( "<h2> Viewer3DImage</h2>" );
  text += "Use the mouse to move the camera around the object. ";
  text += "You can respectively revolve around, zoom and translate with the three mouse buttons. ";
  text += "Left and middle buttons pressed together rotate around the camera view direction axis<br><br>";
  text += "Pressing <b>Alt</b> and one of the function keys (<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
  text += "Simply press the function key again to restore it. Several keyFrames define a ";
  text += "camera path. Paths are saved when you quit the application and restored at next start.<br><br>";
  text += "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
  text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save a snapshot. ";
  text += "See the <b>Keyboard</b> tab in this window for a complete shortcut list.<br><br>";
  text += "Double clicks automates single click actions: A left button double click aligns the closer axis with the camera (if close enough). ";
  text += "A middle button double click fits the zoom of the camera and the right button re-centers the scene.<br><br>";
  text += "A left button double click while holding right button pressed defines the camera <i>Revolve Around Point</i>. ";
  text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
  text += "Press <b>Escape</b> to exit the viewer.";
  return text;
}






void
Viewer3DImage::keyPressEvent ( QKeyEvent *e )
{
  
  bool handled = false;
 
  if( e->key() == Qt::Key_I){
    std::cout << "Image generation" << std::endl;
    handled=true;
  }
  
  if( e->key() == Qt::Key_X){
    std::cout << "Current axis set to X.Image generation" << std::endl;
    myCurrentSliceDim=0;
    handled=true;
  }  
  if( e->key() == Qt::Key_Y){
    std::cout << "Current axis set to Y. Image generation" << std::endl;
    myCurrentSliceDim=1;
    handled=true;
  }
  if( e->key() == Qt::Key_Z){
    std::cout << "Current axis set to Z. Image generation" << std::endl;
    myCurrentSliceDim=2;
    handled=true;
  }
  if( e->key() == Qt::Key_Up ||  e->key() == Qt::Key_Down){
    int dirStep = (e->key() == Qt::Key_Up)? 5: -5; 
    if((e->modifiers() & Qt::ShiftModifier)){
      dirStep/=5;
    }
    int aSliceNum=0;
    int aSliceMax=0;
    bool stoped=false;
    if(myCurrentSliceDim==0){
      aSliceMax=my3dImage->extent()[0];
      if(mySliceXPos+dirStep<aSliceMax&&mySliceXPos+dirStep>=0){
	 mySliceXPos+=dirStep;
      }else{
	stoped=true;
      }
      aSliceNum=mySliceXPos;
    }else if(myCurrentSliceDim==1){
      aSliceMax=my3dImage->extent()[1];
      if(mySliceYPos+dirStep<aSliceMax&&mySliceYPos+dirStep>=0){
	 mySliceYPos+=dirStep;
      }else{
	stoped=true;
      }
      aSliceNum=mySliceYPos;
    }else if(myCurrentSliceDim==2){
       aSliceMax=my3dImage->extent()[2];
       if(mySliceZPos+dirStep<aSliceMax&&mySliceZPos+dirStep>=0){
	 mySliceZPos+=dirStep;
       }else{
	 stoped=true;
       }
       aSliceNum=mySliceZPos;
    }
    if(!stoped){
      std::cerr << "slice num =" << aSliceNum << std::endl;
      
	Image2D sliceImage = DGtal::extractLowerDimImage<Image3D, Image2D>(*my3dImage, myCurrentSliceDim, aSliceNum);
	(*this) << DGtal::UpdateImageData<Image2D>(myCurrentSliceDim, sliceImage, 
						   (myCurrentSliceDim==0)? dirStep: 0.0, 
						   (myCurrentSliceDim==1)? dirStep: 0.0,
						   (myCurrentSliceDim==2)? dirStep: 0.0);
	(*this).updateList(false);
	(*this).update();
      
    }
    handled=true;
   }

   
  if ( !handled )
    DGtal::Viewer3D::keyPressEvent ( e );
  
}



