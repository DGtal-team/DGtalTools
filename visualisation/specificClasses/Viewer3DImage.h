
#if defined(Viewer3DImage_RECURSES)
#error Recursive header files inclusion detected in Viewer3DImage.h
#else // defined(Viewer3DImage_RECURSES)
/** Prevents recursive inclusion of headers. */
#define Viewer3DImage_RECURSES


#if !defined Viewer3DImage_h
#define Viewer3DImage_h
/** Prevents repeated inclusion of headers. */



#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/Display3D.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ConstImageAdapter.h"

template < typename Space = DGtal::Z3i::Space, typename KSpace = DGtal::Z3i::KSpace>
class Viewer3DImage: public DGtal::Viewer3D <Space, KSpace>
{

  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned char> Image3D;
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image2D;
 
  typedef DGtal::ConstImageAdapter<Image3D, Image2D::Domain, DGtal::functors::Projector< DGtal::Z3i::Space>,
                                   Image3D::Value,  DGtal::functors::Identity >  SliceImageAdapter;

  typedef DGtal::ConstImageAdapter<Image3D, DGtal::Z2i::Domain, DGtal::functors::SliceRotator2D< DGtal::Z3i::Domain >,
				   Image3D::Value,  DGtal::functors::Identity >  MyRotatorSliceImageAdapter;



public:
  enum ModeVisu { BoundingBox, InterGrid, Grid, Empty};
  
  Viewer3DImage(ModeVisu aMode=BoundingBox) : DGtal::Viewer3D<Space, KSpace>(), myImageOrigin(DGtal::Z3i::Point(0,0,0)){
    mySliceXPos=0;
    mySliceYPos=0;
    mySliceZPos=0;
    myTotalAngleRotationX=0.0;
    myTotalAngleRotationY=0.0;
    myTotalAngleRotationZ=0.0;
    myAngleRotation=0.0;
    myCurrentSliceDim=0;
    myMode=aMode;
    myDisplayingInfo=true;
  }
  
  
  Viewer3DImage(Image3D * an3DImage){
    my3dImage = an3DImage; 
    Viewer3DImage();
  }
  
  void setVolImage(Image3D * an3DImage);
  
  

  
protected:
  virtual QString helpString() const;  
  virtual void keyPressEvent ( QKeyEvent *e );
  virtual void init();
  Image3D *my3dImage;
  float myScaleX;
  float myScaleY;
  float myScaleZ;
  int mySliceXPos;
  int mySliceYPos;
  int mySliceZPos;
  int myCurrentSliceDim ;
  double myAngleRotation;
  double myTotalAngleRotationX;
  double myTotalAngleRotationY;
  double myTotalAngleRotationZ;
  bool myDisplayingInfo;
  DGtal::Z3i::Point myImageOrigin; 
  ModeVisu myMode;
};


#endif // undefined viewer3dimage

#undef Viewer3DImage_RECURSES
#endif // else defined(Viewer3DImage_RECURSES)

