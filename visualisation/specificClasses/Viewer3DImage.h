#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/Display3D.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/ImageContainerBySTLVector.h"


class Viewer3DImage: public DGtal::Viewer3D 
{

  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned char> Image3D;
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image2D;
 





public:
  
  Viewer3DImage(){
    mySliceXPos=0;
    mySliceYPos=0;
    mySliceZPos=0;
    myCurrentSliceDim=0;
    Viewer3D();
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

  int mySliceXPos;
  int mySliceYPos;
  int mySliceZPos;

  int myCurrentSliceDim ;

};




