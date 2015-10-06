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
#ifndef Q_MOC_RUN
#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/Color.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#endif
#ifdef WITH_ITK
#include "DGtal/io/readers/DicomReader.h"
#endif

#include "sliceViewer.h"
#include "ui_sliceViewer.h"


#ifndef Q_MOC_RUN
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/viewers/DrawWithViewer3DModifier.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/images/ConstImageAdapter.h"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#endif

using namespace std;
using namespace DGtal;
using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;


// Set to define slider int value and grid size

static const int MIN_ZOOM_FACTOR = 10.0;
static const int MAX_ZOOM_FACTOR = 40.0;
static const int INIT_SCALE1_ZOOM_FACTOR = 20.0;



template <typename TImage>
static QImage
getImage(const TImage &anImage, double gridSize, const MainWindow::ColorMapFunctor &colFunctor ){
  typedef ConstImageAdapter<TImage, typename TImage::Domain,
                            functors::BasicDomainSubSampler<typename TImage::Domain, int, double>,
                            unsigned int,
                            MainWindow::ColorMapFunctor> ConstImageAdapterForSubSampling;

  std::vector<double> scales;
  scales.push_back(gridSize);
  scales.push_back(gridSize);  
  functors::BasicDomainSubSampler<typename TImage::Domain,  int, double> subSampler (anImage.domain(),
                                                                                             scales, Z2i::Point(0,0));
  typename TImage::Domain newDomain = subSampler.getSubSampledDomain();
  ConstImageAdapterForSubSampling  scaledImage (anImage, newDomain, subSampler, colFunctor );
  unsigned int height = scaledImage.domain().upperBound()[1]-scaledImage.domain().lowerBound()[1];
  unsigned int width = scaledImage.domain().upperBound()[0]-scaledImage.domain().lowerBound()[0];
  QImage res (width, height,QImage::Format_RGB32 );
  for(unsigned int i=0; i<height; i++){
    for(unsigned int j=0; j<width; j++){
      res.setPixel(j, height-i-1, scaledImage(Z2i::Point(j,i)+scaledImage.domain().lowerBound()));
    }
  }
  return  res;
}


MainWindow::MainWindow(DGtal::Viewer3D<> *aViewer,
                       DGtal::ImageContainerBySTLVector < DGtal::Z3i::Domain, unsigned char > *anImage,
                       const ColorMapFunctor &aFunctor, QWidget *parent, Qt::WindowFlags flags) :
  QMainWindow(parent),
  ui(new Ui::MainWindow),
  myViewer(aViewer),
  myImage3D(anImage),
  myColorMap(aFunctor)
{

  ui->setupUi(this);
  ui->verticalLayout_5->addWidget(aViewer);

  QObject::connect(ui->_horizontalSliderX, SIGNAL(valueChanged(int)), this, SLOT(updateSliceImageX()));
  QObject::connect(ui->_horizontalSliderY, SIGNAL(valueChanged(int)), this, SLOT(updateSliceImageY()));
  QObject::connect(ui->_horizontalSliderZ, SIGNAL(valueChanged(int)), this, SLOT(updateSliceImageZ()));
  QObject::connect(ui->_zoomXSlider, SIGNAL(valueChanged(int)), this, SLOT(updateZoomImageX()));
  QObject::connect(ui->_zoomYSlider, SIGNAL(valueChanged(int)), this, SLOT(updateZoomImageY()));
  QObject::connect(ui->_zoomZSlider, SIGNAL(valueChanged(int)), this, SLOT(updateZoomImageZ()));

  QObject::connect(ui->_scale1ButtonX, SIGNAL(clicked()), this, SLOT(setScale1_1_ImageX()));
  QObject::connect(ui->_scale1ButtonY, SIGNAL(clicked()), this, SLOT(setScale1_1_ImageY()));
  QObject::connect(ui->_scale1ButtonZ, SIGNAL(clicked()), this, SLOT(setScale1_1_ImageZ()));

  QObject::connect(ui->_CoolButton, SIGNAL(clicked()), this, SLOT(changeCoolColorMap()));
  QObject::connect(ui->_HotButton, SIGNAL(clicked()), this, SLOT(changeHotColorMap()));
  QObject::connect(ui->_HueButton, SIGNAL(clicked()), this, SLOT(changeHueColorMap()));
  QObject::connect(ui->_NormalButton, SIGNAL(clicked()), this, SLOT(changeNormalColorMap()));

  ui->_horizontalSliderZ->setMinimum(anImage->domain().lowerBound()[2]);
  ui->_horizontalSliderZ->setMaximum(anImage->domain().upperBound()[2]);

  ui->_horizontalSliderY->setMinimum(anImage->domain().lowerBound()[1]);
  ui->_horizontalSliderY->setMaximum(anImage->domain().upperBound()[1]);

  ui->_horizontalSliderX->setMinimum(anImage->domain().lowerBound()[0]);
  ui->_horizontalSliderX->setMaximum(anImage->domain().upperBound()[0]);

  ui->_zoomXSlider->setMinimum( MIN_ZOOM_FACTOR);
  ui->_zoomXSlider->setMaximum( MAX_ZOOM_FACTOR);
  ui->_zoomXSlider->setValue(INIT_SCALE1_ZOOM_FACTOR);

  ui->_zoomYSlider->setMinimum(MIN_ZOOM_FACTOR);
  ui->_zoomYSlider->setMaximum(MAX_ZOOM_FACTOR);
  ui->_zoomYSlider->setValue(INIT_SCALE1_ZOOM_FACTOR);

  ui->_zoomZSlider->setMinimum(MIN_ZOOM_FACTOR);
  ui->_zoomZSlider->setMaximum(MAX_ZOOM_FACTOR);
  ui->_zoomZSlider->setValue(INIT_SCALE1_ZOOM_FACTOR);

  

}

MainWindow::~MainWindow()
{
  delete ui;
}

void MainWindow::updateAllDisplayedImages(){
  updateSliceImageX();
  updateSliceImageY();
  updateSliceImageZ();
}  
void MainWindow::changeNormalColorMap(){
  myColorMap = MainWindow::Id ; 
  updateAllDisplayedImages();
}

void MainWindow::changeCoolColorMap(){
  myColorMap = MainWindow::GradientMapCool; 
  updateAllDisplayedImages();
}
void MainWindow::changeHueColorMap(){
  myColorMap = MainWindow::HueshadeCM;
  updateAllDisplayedImages();
}
void MainWindow::changeHotColorMap(){
  myColorMap = MainWindow::GradientMapHot;
  updateAllDisplayedImages();
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


void MainWindow::setScale1_1_ImageX(){
  ui->_zoomXSlider->setValue(INIT_SCALE1_ZOOM_FACTOR);
  updateZoomImageX();
}

void MainWindow::setScale1_1_ImageY(){
  ui->_zoomYSlider->setValue(INIT_SCALE1_ZOOM_FACTOR);
  updateZoomImageY();
}

void MainWindow::setScale1_1_ImageZ(){
  ui->_zoomZSlider->setValue(INIT_SCALE1_ZOOM_FACTOR);
  updateZoomImageZ();
}



void MainWindow::updateZoomImageX(){
  double gridSize = (double)INIT_SCALE1_ZOOM_FACTOR/ui->_zoomXSlider->value();
  updateZoomImageX(ui->_horizontalSliderX->value(), gridSize );
  QString gridStr = QString::number(gridSize, 'f', 3);
  QString scaleStr = QString::number(1.0/gridSize, 'f', 3);
  ui->_groupBoxX->setTitle(QString("Slice View X: sampling grid size: ").append(gridStr).
                           append(QString(" (zoom x "). append(scaleStr).append(QString(")") )));
}
void MainWindow::updateZoomImageY(){
  double gridSize = (double)INIT_SCALE1_ZOOM_FACTOR/ui->_zoomYSlider->value();
  updateZoomImageY(ui->_horizontalSliderY->value(), gridSize );
  QString gridStr = QString::number(gridSize, 'f', 3);
  QString scaleStr = QString::number(1.0/gridSize, 'f', 3);
  ui->_groupBoxY->setTitle(QString("Slice View Y: sampling grid size: ").append(gridStr).
                           append(QString(" (zoom x "). append(scaleStr).append(QString(")") )));

}
void MainWindow::updateZoomImageZ(){
  double gridSize = (double)INIT_SCALE1_ZOOM_FACTOR/ui->_zoomZSlider->value();
  updateZoomImageZ(ui->_horizontalSliderZ->value(), gridSize );
  QString gridStr = QString::number(gridSize, 'f', 3);
  QString scaleStr = QString::number(1.0/gridSize, 'f', 3);
  ui->_groupBoxZ->setTitle(QString("Slice View Z: sampling grid size: ").append(gridStr).
                           append(QString(" (zoom x "). append(scaleStr).append(QString(")") )));

}


void MainWindow::updateZoomImageX(unsigned int sliceNumber, double gridSize){
  DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(0);
  DGtal::Z2i::Domain domain2D(invFunctor(myImage3D->domain().lowerBound()),
                              invFunctor(myImage3D->domain().upperBound()));
  DGtal::functors::Projector<DGtal::Z3i::Space> aSliceFunctor(sliceNumber); aSliceFunctor.initAddOneDim(0);
  SliceImageAdapter sliceImage(*myImage3D, domain2D, aSliceFunctor, functors::Identity());
  QImage anImage = getImage(sliceImage, gridSize, myColorMap);
  setImageProjX(QPixmap::fromImage(anImage));
}

void MainWindow::updateZoomImageY(unsigned int sliceNumber, double gridSize){
  DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(1);
  DGtal::Z2i::Domain domain2D(invFunctor(myImage3D->domain().lowerBound()),
                              invFunctor(myImage3D->domain().upperBound()));
  DGtal::functors::Projector<DGtal::Z3i::Space> aSliceFunctor(sliceNumber); aSliceFunctor.initAddOneDim(1);
  SliceImageAdapter sliceImage(*myImage3D, domain2D, aSliceFunctor, functors::Identity());

  QImage anImage = getImage(sliceImage, gridSize, myColorMap);
  setImageProjY(QPixmap::fromImage(anImage));
}


void MainWindow::updateZoomImageZ(unsigned int sliceNumber, double gridSize){
  DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(2);
  DGtal::Z2i::Domain domain2D(invFunctor(myImage3D->domain().lowerBound()),
                              invFunctor(myImage3D->domain().upperBound()));
  DGtal::functors::Projector<DGtal::Z3i::Space> aSliceFunctor(sliceNumber); aSliceFunctor.initAddOneDim(2);
  SliceImageAdapter sliceImage(*myImage3D, domain2D, aSliceFunctor, functors::Identity());
  QImage anImage = getImage(sliceImage, gridSize, myColorMap );
  setImageProjZ(QPixmap::fromImage(anImage));
}


void MainWindow::updateSliceImageX(int sliceNumber, bool init){
  DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(0);
  DGtal::Z2i::Domain domain2D(invFunctor(myImage3D->domain().lowerBound()),
                              invFunctor(myImage3D->domain().upperBound()));
  DGtal::functors::Projector<DGtal::Z3i::Space> aSliceFunctor(sliceNumber); aSliceFunctor.initAddOneDim(0);
  SliceImageAdapter sliceImage (*myImage3D, domain2D, aSliceFunctor, functors::Identity());
  
  double gridSize = ((double)INIT_SCALE1_ZOOM_FACTOR)/ui->_zoomXSlider->value();
  QImage anImage = getImage(sliceImage, gridSize,myColorMap);
  setImageProjX(QPixmap::fromImage(anImage));
  Z3i::Point imageOrigin = myImage3D->domain().lowerBound();
  if(init){
    (*myViewer) << DGtal::AddTextureImage2DWithFunctor<SliceImageAdapter, ColorMapFunctor, Z3i::Space, Z3i::KSpace>(sliceImage, myColorMap, DGtal::Viewer3D<Z3i::Space, Z3i::KSpace>::RGBMode);
    (*myViewer) << DGtal::UpdateImagePosition< Space, KSpace >(0, DGtal::Viewer3D<>::xDirection, sliceNumber,
                                                               imageOrigin[1], imageOrigin[2]);
    (*myViewer) << Viewer3D<>::updateDisplay;
  }else{
    (*myViewer) << DGtal::UpdateImageData< SliceImageAdapter, ColorMapFunctor > (0, sliceImage, 0, 0, 0 ,0,  DGtal::Viewer3D<>::xDirection, myColorMap);
    (*myViewer) << DGtal::UpdateImagePosition< Space, KSpace >(0, DGtal::Viewer3D<>::xDirection, sliceNumber, imageOrigin[1],
                                                               imageOrigin[2]);
    
    (*myViewer).updateList(init);
    (*myViewer).update();
  }
  
  
}


void MainWindow::updateSliceImageY( int sliceNumber, bool init){
  
  DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(1);
  DGtal::Z2i::Domain domain2D(invFunctor(myImage3D->domain().lowerBound()),
                              invFunctor(myImage3D->domain().upperBound()));
  DGtal::functors::Projector<DGtal::Z3i::Space> aSliceFunctor(sliceNumber); aSliceFunctor.initAddOneDim(1);
  SliceImageAdapter sliceImage(*myImage3D, domain2D, aSliceFunctor, functors::Identity());
  
  double gridSize = ((double)INIT_SCALE1_ZOOM_FACTOR)/ui->_zoomYSlider->value();
  QImage anImage = getImage(sliceImage, gridSize, myColorMap);
  setImageProjY(QPixmap::fromImage(anImage));
  Z3i::Point imageOrigin = myImage3D->domain().lowerBound();
  if(init){
    (*myViewer) << DGtal::AddTextureImage2DWithFunctor<SliceImageAdapter, ColorMapFunctor, Z3i::Space, Z3i::KSpace>(sliceImage, myColorMap, DGtal::Viewer3D<Z3i::Space, Z3i::KSpace>::RGBMode);
    (*myViewer) << DGtal::UpdateImagePosition< Space, KSpace >(1, DGtal::Viewer3D<>::yDirection, imageOrigin[0],
                                                               sliceNumber, imageOrigin[2]);
    (*myViewer) << Viewer3D<>::updateDisplay;
  }else{
    (*myViewer) << DGtal::UpdateImageData< SliceImageAdapter, ColorMapFunctor > (1, sliceImage, 0,0, 0, 0,  DGtal::Viewer3D<>::yDirection, myColorMap);
    (*myViewer) << DGtal::UpdateImagePosition< Space, KSpace >(1, DGtal::Viewer3D<>::yDirection, imageOrigin[0],
                                                               sliceNumber, imageOrigin[2]);
    (*myViewer).updateList(init);
    (*myViewer).update();
  }
  
}


void MainWindow::updateSliceImageZ(int sliceNumber, bool init){
  
  DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(2);
  DGtal::Z2i::Domain domain2D(invFunctor(myImage3D->domain().lowerBound()),
                              invFunctor(myImage3D->domain().upperBound()));
  DGtal::functors::Projector<DGtal::Z3i::Space> aSliceFunctor(sliceNumber); aSliceFunctor.initAddOneDim(2);
  SliceImageAdapter sliceImage(*myImage3D, domain2D, aSliceFunctor, functors::Identity());
  double gridSize = (double)INIT_SCALE1_ZOOM_FACTOR/ui->_zoomZSlider->value();
  QImage anImage = getImage(sliceImage, gridSize, myColorMap);
  setImageProjZ(QPixmap::fromImage(anImage));
  Z3i::Point imageOrigin = myImage3D->domain().lowerBound();
  if(init){
    (*myViewer) << DGtal::AddTextureImage2DWithFunctor<SliceImageAdapter, ColorMapFunctor, Z3i::Space, Z3i::KSpace>(sliceImage, myColorMap, DGtal::Viewer3D<Z3i::Space, Z3i::KSpace>::RGBMode);
    (*myViewer) << DGtal::UpdateImagePosition< Space, KSpace >(2, DGtal::Viewer3D<>::zDirection, imageOrigin[0],
                                                               imageOrigin[1], sliceNumber);
  
    (*myViewer) << Viewer3D<>::updateDisplay;
  }else{
    (*myViewer) << DGtal::UpdateImageData< SliceImageAdapter,ColorMapFunctor > (2, sliceImage, 0,0, 0, 0,  DGtal::Viewer3D<>::zDirection,
                                                                                myColorMap);
    (*myViewer) << DGtal::UpdateImagePosition< Space, KSpace >(2, DGtal::Viewer3D<>::zDirection, imageOrigin[0],
                                                               imageOrigin[1], sliceNumber);
    (*myViewer).updateList(init);
    (*myViewer).update();
  }

      
}



int main( int argc, char** argv )
{

  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
    ("hueColorMap", "use hue color map to display images." )
    ("gradHotColorMap", "use hot gradient color map to display images." )
    ("gradCoolColorMap", "use cool gradient color map to display images." )

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
      std::cout << "Usage: " << argv[0] << " [input]\n"
                << "Display volume file as a voxel set by using QGLviewer"<< endl
                << general_opt << "\n";
      return 0;
    }

  if(! vm.count("input"))
    {
      trace.error() << " The file name was defined" << endl;
      return 0;
    }
  string inputFilename = vm["input"].as<std::string>();


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
    bool usehm = vm.count("hueColorMap");
    bool usegh = vm.count("gradHotColorMap");
    bool usegc = vm.count("gradCoolColorMap");
    
    MainWindow w(viewer, &image, MainWindow::ColorMapFunctor(usehm? MainWindow::HueshadeCM:
                                                             usegh? MainWindow::GradientMapHot:
                                                             usegc? MainWindow::GradientMapCool:
                                                             MainWindow::Id), 0,0);
    w.setWindowTitle ( QString("sliceViewer"));
    w.updateSliceImageX(0, true);
    w.updateSliceImageY(0, true);
    w.updateSliceImageZ(0, true);
    w.show();
    Z3i::Point size = image.domain().upperBound() - image.domain().lowerBound();
    Z3i::Point center = image.domain().lowerBound()+size/2;
    unsigned int maxDist = std::max(std::max(size[2], size[1]), size[2]);
    viewer->camera()->setPosition(qglviewer::Vec(center[0],center[1], 
                                                 2.0*maxDist));
    viewer->camera()->setSceneCenter(qglviewer::Vec(center[0],center[1],center[2]));
    application.exec();
    delete viewer;
  }
}


