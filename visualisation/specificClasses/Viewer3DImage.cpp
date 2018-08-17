
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
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/images/ImageSelector.h"
#include <sstream>      
#include <QKeyEvent>

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace qglviewer;


///////////////////////////////////////////////////////////////////////////////
// class Viewer3DImage
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :


template < typename Space, typename KSpace>
void
Viewer3DImage< Space, KSpace>::init(){
   DGtal::Viewer3D<>::init();

   QGLViewer::setKeyDescription ( Qt::Key_X, "Change the current axis to X for the current 2D image slice setting." );
   QGLViewer::setKeyDescription ( Qt::Key_Y, "Change the current axis to Y for the current 2D image slice setting." );
   QGLViewer::setKeyDescription ( Qt::Key_Z, "Change the current axis to Z for the current 2D image slice setting." );
   QGLViewer::setKeyDescription ( Qt::Key_Up, "Move the current 2D image slice to 5 in the positive direction of the current axis." );
   QGLViewer::setKeyDescription ( Qt::Key_Down, "Move the current 2D image slice to 5 in the negative direction of the current axis." );
   QGLViewer::setKeyDescription ( Qt::Key_Shift, "Change the slice move with step 1 (5 by default)" );
   QGLViewer::setKeyDescription ( Qt::Key_M, "Hide/Display message informations (slice numbers, selected axis)" );   
}


template < typename Space, typename KSpace>
void 
Viewer3DImage< Space, KSpace>::setVolImage(Image3D * an3DImage){
  my3dImage = an3DImage;

  switch (myMode) {
  case BoundingBox:
    (*this) << DGtal::SetMode3D("ConstImageAdapter", "BoundingBox");
    break;
  case InterGrid:
    (*this) << DGtal::SetMode3D("ConstImageAdapter", "InterGrid");
    break;
  case Grid:
    (*this) << DGtal::SetMode3D("ConstImageAdapter", "Grid");
    break;
  case Empty:
    (*this) << DGtal::SetMode3D("ConstImageAdapter", "");
    break;
  }

  myImageOrigin = my3dImage->domain().lowerBound();
  (*this).updateList(false);
  (*this).update();
  
  mySliceXPos=myImageOrigin[0];
  mySliceYPos=myImageOrigin[1];
  mySliceZPos=myImageOrigin[2];

 // Adding X slice in the viewer.
  DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctorX; invFunctorX.initRemoveOneDim(0);
  DGtal::Z2i::Domain domain2DX(invFunctorX(my3dImage->domain().lowerBound()), 
			       invFunctorX(my3dImage->domain().upperBound()));
  
    
  DGtal::functors::SliceRotator2D<DGtal::Z3i::Domain> aSliceFunctorX(0, my3dImage->domain(), mySliceXPos,2, myAngleRotation );
  const DGtal::functors::Identity identityFunctor{};
  MyRotatorSliceImageAdapter sliceImageX( *my3dImage, domain2DX, aSliceFunctorX, identityFunctor ); 
    
  std::cout << "image:" << sliceImageX.className();
  (*this) << sliceImageX;
  (*this) << DGtal::UpdateImagePosition< Space, KSpace >(0, DGtal::Viewer3D<>::xDirection, mySliceXPos,myImageOrigin[1], myImageOrigin[2]);


  // Adding Y slice in the viewer.
  DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctorY; invFunctorY.initRemoveOneDim(1);
  DGtal::Z2i::Domain domain2DY(invFunctorY(my3dImage->domain().lowerBound()), 
			       invFunctorY(my3dImage->domain().upperBound()));
 
  DGtal::functors::Projector<DGtal::Z3i::Space> aSliceFunctorY(mySliceYPos); aSliceFunctorY.initAddOneDim(1);
  SliceImageAdapter sliceImageY(*my3dImage, domain2DY, aSliceFunctorY, identityFunctor ); 
  (*this) << sliceImageY;
  (*this) << DGtal::UpdateImagePosition< Space, KSpace >(1, DGtal::Viewer3D<>::yDirection, myImageOrigin[0], mySliceYPos, myImageOrigin[2]);



 
 
   // Adding Z slice in the viewer.
  DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctorZ; invFunctorZ.initRemoveOneDim(2);
  DGtal::Z2i::Domain domain2DZ(invFunctorZ(my3dImage->domain().lowerBound()), 
			       invFunctorZ(my3dImage->domain().upperBound()));
 
  DGtal::Z3i::Point centerZ((my3dImage->domain().upperBound())[0]/2, (my3dImage->domain().upperBound())[1]/2, mySliceZPos);
    
  DGtal::functors::SliceRotator2D<DGtal::Z3i::Domain> aSliceFunctorZ(2, my3dImage->domain(), mySliceZPos, 2, myAngleRotation, centerZ );
  MyRotatorSliceImageAdapter sliceImageZ( *my3dImage, domain2DZ, aSliceFunctorZ, identityFunctor ); 
  (*this) << sliceImageZ;
  (*this) << DGtal::UpdateImagePosition< Space, KSpace > (2, DGtal::Viewer3D<>::zDirection, myImageOrigin[0], myImageOrigin[1], mySliceZPos);
  (*this) << DGtal::Viewer3D<>::updateDisplay;    
  
}


template < typename Space, typename KSpace>
QString
Viewer3DImage< Space, KSpace>::helpString() const
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





template < typename Space, typename KSpace>
void
Viewer3DImage< Space, KSpace>::keyPressEvent ( QKeyEvent *e )
{
  
  bool handled = false;
  if( e->key() == Qt::Key_M){
    myDisplayingInfo = !myDisplayingInfo;
    handled=true;
  } 
  if( e->key() == Qt::Key_I){
    std::cout << "Image generation" << std::endl;
    handled=true;
  }
  
  if( e->key() == Qt::Key_X){
    std::cout << "Current axis set to X" << std::endl;
    (*this).displayMessage(QString("Current axis set to X"), 100000); 
    myCurrentSliceDim=0;
    handled=true;
  }  
  if( e->key() == Qt::Key_Y){
    (*this).displayMessage(QString("Current axis set to Y"), 100000); 
    myCurrentSliceDim=1;
    handled=true;
  }
  if( e->key() == Qt::Key_Z){
    (*this).displayMessage(QString("Current axis set to Z"), 100000); 
    myCurrentSliceDim=2;
    handled=true;
  }
  if( e->key() == Qt::Key_Up ||  e->key() == Qt::Key_Down){
    int dirStep = (e->key() == Qt::Key_Up)? 5: -5; 
    if((e->modifiers() & Qt::ShiftModifier)){
      dirStep/=5;
    }
    if((e->modifiers() & Qt::MetaModifier)){
      myAngleRotation =(e->key() == Qt::Key_Down)? 0.1: -0.1;
      dirStep=0;
    }else{
      myAngleRotation=0.0;
    }
    int aSliceNum=0;
    int aSliceMax=0;
    bool stoped=false;
    if(myCurrentSliceDim==0){
      aSliceMax=my3dImage->domain().upperBound()[0]+1;
      if(mySliceXPos+dirStep<aSliceMax&&mySliceXPos+dirStep>=myImageOrigin[0]){
	 mySliceXPos+=dirStep;
      }else{
	stoped=true;
      }

      aSliceNum=mySliceXPos;
    }else if(myCurrentSliceDim==1){
      aSliceMax=my3dImage->domain().upperBound()[1]+1;
      if(mySliceYPos+dirStep<aSliceMax&&mySliceYPos+dirStep>=myImageOrigin[1]){
	 mySliceYPos+=dirStep;
      }else{
	stoped=true;
      }
      aSliceNum=mySliceYPos;
    }else if(myCurrentSliceDim==2){
      aSliceMax=my3dImage->domain().upperBound()[2]+1;
       if(mySliceZPos+dirStep<aSliceMax&&mySliceZPos+dirStep>=myImageOrigin[2]){
	 mySliceZPos+=dirStep;
       }else{
	 stoped=true;
       }
       aSliceNum=mySliceZPos;
    }

    if(!stoped){
      double myTotalAngleRotation=0.0;
      if(myCurrentSliceDim==0){
	myTotalAngleRotationX+=myAngleRotation;
	myTotalAngleRotation = myTotalAngleRotationX;
      }else if (myCurrentSliceDim==1){
	myTotalAngleRotationY+=myAngleRotation;
	myTotalAngleRotation = myTotalAngleRotationY;
      }else if (myCurrentSliceDim==2){
	myTotalAngleRotationZ+=myAngleRotation;
	myTotalAngleRotation = myTotalAngleRotationZ;
      }
      
      DGtal::functors::SliceRotator2D<DGtal::Z3i::Domain> aSliceFunctor(myCurrentSliceDim, my3dImage->domain(), 
							      aSliceNum, 2, myTotalAngleRotation);
      // Adding X slice in the viewer.
      DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(myCurrentSliceDim);
      DGtal::Z2i::Domain domain2D(invFunctor(my3dImage->domain().lowerBound()), 
				  invFunctor(my3dImage->domain().upperBound()));
     
      const DGtal::functors::Identity identityFunctor{};
      MyRotatorSliceImageAdapter sliceImage( *my3dImage, domain2D, aSliceFunctor, identityFunctor ); 
      
      (*this) << DGtal::UpdateImageData<MyRotatorSliceImageAdapter>(myCurrentSliceDim, sliceImage, 
                                                                   (myCurrentSliceDim==0)? dirStep: 0.0, 
                                                                   (myCurrentSliceDim==1)? dirStep: 0.0,
                                                                   (myCurrentSliceDim==2)? dirStep: 0.0, 
							   myAngleRotation, DGtal::Viewer3D<>::zDirection);
      (*this).updateList(false);
      (*this).update();      
    }
    if(myDisplayingInfo){
      std::stringstream sstring;
      sstring << "slice X: " << mySliceXPos ;
      sstring << " slice Y: " << mySliceYPos ;
      sstring << " slice Z: " << mySliceZPos ;
      (*this).displayMessage(QString(sstring.str().c_str()), 100000);
    }



    
    handled=true;
  }
   
  if ( !handled )
    DGtal::Viewer3D<>::keyPressEvent ( e );
  
}



