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
 * @file sliceViewer.cpp
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2014/07/04
 *
 * An example file named qglViewer.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <QtGui/qapplication.h>

#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#ifdef WITH_ITK
#include "DGtal/io/readers/DicomReader.h"
#endif
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/readers/PointListReader.h"

#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "sliceViewer.h"
#include "ui_sliceViewer.h"


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;
using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;





template <typename TImage>
static QImage 
getImage(const TImage &anImage ){
  unsigned int height = anImage.domain().upperBound()[1]+1;
  unsigned int width = anImage.domain().upperBound()[0]+1;
  uchar * data = new uchar [height*width*4];
  for(unsigned int i=0; i<height; i++){
    for(unsigned int j=0; j<width; j++){
      data[(j+width*i)*4]=anImage(Z2i::Point(j,i));
      data[(j+width*i)*4+1]=anImage(Z2i::Point(j,i));
      data[(j+width*i)*4+2]=anImage(Z2i::Point(j,i));
      data[(j+width*i)*4+3]=anImage(Z2i::Point(j,i));
    }
  }
   QImage result(  data, width,  height, QImage::Format_RGB32 );
  return result;
}


MainWindow::MainWindow(DGtal::Viewer3D<> *aViewer,       
                       DGtal::ImageContainerBySTLVector < DGtal::Z3i::Domain, unsigned char > *anImage,
                       QWidget *parent, Qt::WindowFlags flags) :
    myViewer(aViewer),
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    myImage3D(anImage)
{

    ui->setupUi(this);
    ui->verticalLayout_5->addWidget(aViewer);

    QObject::connect(ui->_horizontalSliderX, SIGNAL(valueChanged(int)), this, SLOT(updateSliceImageX()));
    QObject::connect(ui->_horizontalSliderY, SIGNAL(valueChanged(int)), this, SLOT(updateSliceImageY()));
    QObject::connect(ui->_horizontalSliderZ, SIGNAL(valueChanged(int)), this, SLOT(updateSliceImageZ()));

    ui->_horizontalSliderZ->setMinimum(0);
    ui->_horizontalSliderZ->setMaximum(anImage->domain().upperBound()[2]);

    ui->_horizontalSliderY->setMinimum(0);
    ui->_horizontalSliderY->setMaximum(anImage->domain().upperBound()[1]);

    ui->_horizontalSliderX->setMinimum(0);
    ui->_horizontalSliderX->setMaximum(anImage->domain().upperBound()[0]);
    
}

MainWindow::~MainWindow()
{
     delete ui;
}

void MainWindow::setImageProjX(const QPixmap &aPixMap){
  ui->ImageProjX->setPixmap(aPixMap);  
}
void MainWindow::setImageProjY(const QPixmap &aPixMap){
  ui->ImageProjY->setPixmap(aPixMap);  
}
void MainWindow::setImageProjZ(const QPixmap &aPixMap){
  ui->ImageProjZ->setPixmap(aPixMap);  
}


void MainWindow::updateSliceImageX(){
  updateSliceImageX(ui->_horizontalSliderX->value(), false);
}

void MainWindow::updateSliceImageY(){
  updateSliceImageY(ui->_horizontalSliderY->value(), false);
}

void MainWindow::updateSliceImageZ(){
  updateSliceImageZ(ui->_horizontalSliderZ->value(), false);
}


void MainWindow::updateSliceImageX(unsigned int sliceNumber, bool init){
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char > Image2D;

  typedef DGtal::ConstImageAdapter<Image3D, Image2D::Domain, DGtal::functors::Projector< DGtal::Z3i::Space>,
				   Image3D::Value,  functors::Identity >  SliceImageAdapter;

  DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(0);
  DGtal::Z2i::Domain domain2D(invFunctor(myImage3D->domain().lowerBound()),
                              invFunctor(myImage3D->domain().upperBound()));
  DGtal::functors::Projector<DGtal::Z3i::Space> aSliceFunctor(sliceNumber); aSliceFunctor.initAddOneDim(0);
  SliceImageAdapter sliceImage(*myImage3D, domain2D, aSliceFunctor, functors::Identity());
  QImage anImage = getImage(sliceImage); 
  setImageProjX(QPixmap::fromImage(anImage));
  if(init){
    (*myViewer) << sliceImage;
    (*myViewer) << DGtal::UpdateImagePosition< Space, KSpace >(0, DGtal::Viewer3D<>::xDirection, sliceNumber, 0.0, 0.0);
    (*myViewer) << Viewer3D<>::updateDisplay;
  }else{
    (*myViewer) << DGtal::UpdateImageData< SliceImageAdapter > (0, sliceImage, 0,0, 0,0,  DGtal::Viewer3D<>::xDirection);
    (*myViewer) << DGtal::UpdateImagePosition< Space, KSpace >(0, DGtal::Viewer3D<>::xDirection, sliceNumber, 0.0, 0.0);
    (*myViewer).updateList(init);
    (*myViewer).update();      
  }
  

}


void MainWindow::updateSliceImageY(unsigned int sliceNumber, bool init){
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char > Image2D;

  typedef DGtal::ConstImageAdapter<Image3D, Image2D::Domain, DGtal::functors::Projector< DGtal::Z3i::Space>,
				   Image3D::Value,  functors::Identity >  SliceImageAdapter;

  DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(1);
  DGtal::Z2i::Domain domain2D(invFunctor(myImage3D->domain().lowerBound()),
                              invFunctor(myImage3D->domain().upperBound()));
  DGtal::functors::Projector<DGtal::Z3i::Space> aSliceFunctor(sliceNumber); aSliceFunctor.initAddOneDim(1);
  SliceImageAdapter sliceImage(*myImage3D, domain2D, aSliceFunctor, functors::Identity());
  QImage anImage = getImage(sliceImage); 
  setImageProjY(QPixmap::fromImage(anImage));
  if(init){
    (*myViewer) << sliceImage;
    (*myViewer) << DGtal::UpdateImagePosition< Space, KSpace >(1, DGtal::Viewer3D<>::yDirection, 0.0, sliceNumber, 0.0);
    (*myViewer) << Viewer3D<>::updateDisplay;
  }else{
    (*myViewer) << DGtal::UpdateImageData< SliceImageAdapter > (1, sliceImage, 0,0, 0, 0,  DGtal::Viewer3D<>::yDirection);
    (*myViewer) << DGtal::UpdateImagePosition< Space, KSpace >(1, DGtal::Viewer3D<>::yDirection, 0.0, sliceNumber, 0.0);
    (*myViewer).updateList(init);
    (*myViewer).update();      
  }



}


void MainWindow::updateSliceImageZ(unsigned int sliceNumber, bool init){
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char > Image2D;

  typedef DGtal::ConstImageAdapter<Image3D, Image2D::Domain, DGtal::functors::Projector< DGtal::Z3i::Space>,
				   Image3D::Value,  functors::Identity >  SliceImageAdapter;

  DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(2);
  DGtal::Z2i::Domain domain2D(invFunctor(myImage3D->domain().lowerBound()),
                              invFunctor(myImage3D->domain().upperBound()));
  DGtal::functors::Projector<DGtal::Z3i::Space> aSliceFunctor(sliceNumber); aSliceFunctor.initAddOneDim(2);
  SliceImageAdapter sliceImage(*myImage3D, domain2D, aSliceFunctor, functors::Identity());
  QImage anImage = getImage(sliceImage); 
  setImageProjZ(QPixmap::fromImage(anImage));
 if(init){
    (*myViewer) << sliceImage;
    (*myViewer) << DGtal::UpdateImagePosition< Space, KSpace >(2, DGtal::Viewer3D<>::zDirection, 0.0, 0.0, sliceNumber);
    (*myViewer) << Viewer3D<>::updateDisplay;
  }else{
    (*myViewer) << DGtal::UpdateImageData< SliceImageAdapter > (2, sliceImage, 0,0, 0, 0,  DGtal::Viewer3D<>::zDirection);
    (*myViewer) << DGtal::UpdateImagePosition< Space, KSpace >(2, DGtal::Viewer3D<>::zDirection, 0.0,  0.0, sliceNumber);
    (*myViewer).updateList(init);
    (*myViewer).update();      
 }
  
}



int main( int argc, char** argv )
{

 po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input-file,i", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
#ifdef WITH_ITK
    ("dicomMin", po::value<int>()->default_value(-1000), "set minimum density threshold on Hounsfield scale")
    ("dicomMax", po::value<int>()->default_value(3000), "set maximum density threshold on Hounsfield scale")
#endif    
    ;

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
                << "Display volume file as a voxel set by using QGLviewer"<< endl
                << general_opt << "\n";
      return 0;
    }
  
  if(! vm.count("input-file"))
    {
      trace.error() << " The file name was defined" << endl;      
      return 0;
    }
  string inputFilename = vm["input-file"].as<std::string>();
  

  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char > Image2D;
  
    string extension = inputFilename.substr(inputFilename.find_last_of(".") + 1);
  if(extension!="vol" && extension != "p3d" && extension != "pgm3D" && extension != "pgm3d" && extension != "sdp" && extension != "pgm" 
 #ifdef WITH_ITK
    && extension !="dcm"
#endif
){
    trace.info() << "File extension not recognized: "<< extension << std::endl;
    return 0;
  }
  
  if(extension=="vol" || extension=="pgm3d" || extension=="pgm3D"
#ifdef WITH_ITK
    || extension =="dcm"
#endif
){
    unsigned int numDisplayed=0;
    
#ifdef WITH_ITK
   int dicomMin = vm["dicomMin"].as<int>();
   int dicomMax = vm["dicomMax"].as<int>();
   typedef functors::Rescaling<int ,unsigned char > RescalFCT;
   Image3D image = extension == "dcm" ? DicomReader< Image3D,  RescalFCT  >::importDicom( inputFilename, 
                                                                                      RescalFCT(dicomMin,
                                                                                                dicomMax,
												    0, 255) ) : 
     GenericReader<Image3D>::import( inputFilename );
   trace.info() << "Imported ITK..."<< std::endl;
#else
   Image3D image = GenericReader<Image3D>::import (inputFilename );
   trace.info() << "Imported..."<< std::endl;
#endif

  

  QApplication application(argc,argv);
  Viewer3D<> *viewer = new Viewer3D<>();
  MainWindow w(viewer, &image, 0,0);
  w.setWindowTitle ( QString("sliceViewer"));
  w.updateSliceImageX(0, true);
  w.updateSliceImageY(0, true);
  w.updateSliceImageZ(0, true);


  w.show();
  application.exec();
  delete viewer;
  }
}


