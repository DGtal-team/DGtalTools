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
 * @file 3dCompSurfelData.cpp
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2014/08/01
 *
 * A tool to compare generic surfel data informations given from two data files.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <QGLViewer/qglviewer.h>
#include <stdio.h>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/colormaps/GradientColorMap.h"

#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/CanonicCellEmbedder.h"
#include "DGtal/math/Statistic.h"

#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"

#include "CLI11.hpp"


#include <limits>

using namespace std;
using namespace DGtal;
using namespace Z3i;



/**
 @page CompSurfelData 3DCompSurfelData

 @brief  Computes generic scalar surfel data comparisons (squared error) (given from an input data file and from a reference one). 
 
 

 @b Usage: 3dCompSurfelData [input] [reference]
 

 @b Allowed @b options @b are:
 
 @code

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         input file: sdpa (sequence of discrete points with attribute)
   -r,--reference TEXT:FILE REQUIRED     input file: sdpa (sequence of discrete points with attribute)
   -l,--compAccordingLabels              apply the comparisos only on points with same labels (by default fifth colomn)
   -a,--drawSurfelAssociations           Draw the surfel association.
   -o,--fileMeasureOutput TEXT           specify the output file to store (append) the error stats else the result is given to std output.
   -n,--noWindows                        Don't display Viewer windows.
   -d,--doSnapShotAndExit TEXT           save display snapshot into file. Notes that  the camera setting is set by default according the last saved configuration (use SHIFT+Key_M to save current camera setting in the Viewer3D).
   --fixMaxColorValue FLOAT              fix the maximal color value for the scale error display (else the scale is set from the maximal value)
   --labelIndex UINT                     set the index of the label (by default set to 4)
   --SDPindex UINT x 4                   specify the sdp index (by default 0,1,2,3).


 @endcode

 @b Example: 

 To use this tools  we need first to have two differents surfel set with an attribute to be compared. We can use for instance the curvature from the sample file of DGtal/examples/samples/cat10.vol: 
 @code
 # generating another input vol file using tutorial example (eroded.vol):
 $ $DGtal/build/examples/tutorial-examples/FMMErosion
 # Estimate curvature using other DGtalTools program:
 $ 3dCurvatureViewer -i $DGtal/examples/samples/cat10.vol -r 3 --exportOnly -d curvatureCat10R3.dat
 $ 3dCurvatureViewer -i eroded.vol -r 3 --exportOnly -d curvatureCat10ErodedR3.dat
 @endcode
 Now we compare the different curvature values from the two shapes:
 @code
  ./visualisation/3dCompSurfelData -i curvatureCat10ErodedR3.dat -r curvatureCat10R3.dat --fixMaxColorValue 1.0  -o  statMeasures.dat                       
 @endcode

 You should obtain such a visualisation:
 @image html res3dCompSurfelData.png "resulting visualisation of curvature comparison."
 

 @see
 @ref 3dCompSurfelData.cpp
 @ref 3dCurvatureViewer 
 */



template < typename Space = DGtal::Z3i::Space, typename KSpace = DGtal::Z3i::KSpace>
struct ViewerSnap: DGtal::Viewer3D <Space, KSpace>
{

  ViewerSnap(const KSpace &KSEmb, bool saveSnap): Viewer3D<Space, KSpace>(KSEmb), mySaveSnap(saveSnap){
  };

  virtual  void
  init(){
    DGtal::Viewer3D<>::init();
    if(mySaveSnap){
      QObject::connect(this, SIGNAL(drawFinished(bool)), this, SLOT(saveSnapshot(bool)));
    }
  };
  bool mySaveSnap;
};


template < typename Point>
void
getBoundingUpperAndLowerPoint(const std::vector<Point> &vectorPt, Point &ptLower, Point &ptUpper){
  for(unsigned int i =1; i<vectorPt.size(); i++){
    if(vectorPt.at(i)[0] < ptLower[0]){
      ptLower[0] = vectorPt.at(i)[0];
    }
    if(vectorPt.at(i)[1] < ptLower[1]){
      ptLower[1] = vectorPt.at(i)[1];
    }
   if(vectorPt.at(i)[2] < ptLower[2]){
      ptLower[2] =vectorPt.at(i)[2];
    }
   if(vectorPt.at(i)[0] < ptLower[0]){
      ptLower[0] = vectorPt.at(i)[0];
   }
   if(vectorPt.at(i)[1] < ptLower[1]){
     ptLower[1] = vectorPt.at(i)[1];
    }
   if(vectorPt.at(i)[2] < ptLower[2]){
      ptLower[2] =vectorPt.at(i)[2];
    }
  }
}


int main( int argc, char** argv )
{
  typedef PointVector<4, double> Point4D;
  typedef PointVector<1, int> Point1D;
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string referenceFileName;
  std::string fileMeasureOutput;
  std::string snapShotName;
  bool useLabels {false};
  bool drawSurfelAssociations {false};
  bool noWindows {false};
  double maxVal;
  unsigned int labelIndex {4};
  std::vector<unsigned int> vectSDPindex;
  
  app.description("Computes generic scalar surfel data comparisons (squared error) (given from an input data file and from a reference one).\n Each surfels are associated to the nearest one of the reference surfels (computed by a 'brut force' search). This association can also be limited to surfel of same label (if available in the data and by using the --compAccordingLabels option). The comparison and surfel association can be displayed and result statistics are saved on output file (--fileMeasureOutput). You can also remove the interactive 3d view by just doing a snapshot and exit with option --doSnapShotAndExit.\n \t  Example of use: \n \t \t 3dCompSurfelData -i surfelCurvatureInput.dat -r surfelCurvatureRef.dat  --fixMaxColorValue 0.8 -d visuSEcurvature.png -o  statMeasures.dat \n \n  \t \t => From the two compared files you should obtain the result of the comparison (statMeasures.dat) with the associated visualisation (visuSEcurvature.png).\n");
  
  app.add_option("-i,--input,1", inputFileName, "input file: sdpa (sequence of discrete points with attribute)" )
  ->required()
  ->check(CLI::ExistingFile);
  
  app.add_option("-r,--reference,2", referenceFileName, "input file: sdpa (sequence of discrete points with attribute)" )
  ->required()
  ->check(CLI::ExistingFile);
  app.add_flag("-l,--compAccordingLabels", useLabels, "apply the comparisos only on points with same labels (by default fifth colomn)" );
  app.add_flag("-a,--drawSurfelAssociations", drawSurfelAssociations, "Draw the surfel association.");
  app.add_option("-o,--fileMeasureOutput",fileMeasureOutput,  "specify the output file to store (append) the error stats else the result is given to std output. ");
  
  app.add_flag("-n,--noWindows","Don't display Viewer windows." );
    
  app.add_option("-d,--doSnapShotAndExit", snapShotName, "save display snapshot into file. Notes that  the camera setting is set by default according the last saved configuration (use SHIFT+Key_M to save current camera setting in the Viewer3D).");
  auto maxValSpe = app.add_option("--fixMaxColorValue", maxVal, "fix the maximal color value for the scale error display (else the scale is set from the maximal value)");
  auto labOptio = app.add_option("--labelIndex", labelIndex,"set the index of the label (by default set to 4)  " );
  app.add_option("--SDPindex",vectSDPindex, "specify the sdp index (by default 0,1,2,3).")
   ->expected(4);
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------




  Z3i::KSpace K;

  
  std::vector<Point4D> surfelAndScalarInput;
  std::vector<Point4D> surfelAndScalarReference;
  std::vector<Point1D> vectLabelsInput;
  std::vector<Point1D> vectLabelsReference;

  if(useLabels){
    if(labOptio->count() > 0){
      std::vector<unsigned int > vectIndex;
      vectIndex.push_back(labelIndex);
      vectLabelsInput = PointListReader<Point1D>::getPointsFromFile(inputFileName, vectIndex);
      vectLabelsReference = PointListReader<Point1D>::getPointsFromFile(referenceFileName, vectIndex);
    }else{
      vectLabelsInput = PointListReader<Point1D>::getPointsFromFile(inputFileName);
      vectLabelsReference = PointListReader<Point1D>::getPointsFromFile(referenceFileName);
    }
  }


  if(vectSDPindex.size() == 4) {
    surfelAndScalarInput = PointListReader<Point4D>::getPointsFromFile(inputFileName, vectSDPindex);
    surfelAndScalarReference = PointListReader<Point4D>::getPointsFromFile(referenceFileName, vectSDPindex);
  }else{
    surfelAndScalarInput = PointListReader<Point4D>::getPointsFromFile(inputFileName);
    surfelAndScalarReference = PointListReader<Point4D>::getPointsFromFile(referenceFileName);
  }


  Point4D ptLower = surfelAndScalarReference.at(0);
  Point4D ptUpper = surfelAndScalarReference.at(0);
  getBoundingUpperAndLowerPoint(surfelAndScalarReference,  ptLower, ptUpper);
  getBoundingUpperAndLowerPoint(surfelAndScalarInput,  ptLower, ptUpper);

  K.init(Z3i::Point(2*ptLower[0]+1, 2*ptLower[1]+1, 2*ptLower[2]+1),
         Z3i::Point(2*ptUpper[0]+1, 2*ptUpper[1]+1, 2*ptUpper[2]+1), true);


  std::vector<Cell> vectSurfelsReference;
  std::vector<Cell> vectSurfelsInput;

  // Construction of the set of surfels
  for(unsigned int i =0; i<surfelAndScalarReference.size(); i++){
    Point4D pt4d = surfelAndScalarReference.at(i);
    Cell c = K.uCell(Z3i::Point(pt4d[0], pt4d[1], pt4d[2]));
    vectSurfelsReference.push_back(c);
  }
  // Construction of the set of surfels
  for(unsigned int i =0; i<surfelAndScalarInput.size(); i++){
    Point4D pt4d = surfelAndScalarInput.at(i);
    Cell c = K.uCell(Z3i::Point(pt4d[0], pt4d[1], pt4d[2]));
    vectSurfelsInput.push_back(c);
  }


  // Brut force association of each surfel of the input to the reference with the minimal distance.
  // For each surfel of the input we associate an index to the nearest surfel of the reference.
  // if the option --compAccordingLabels is selected then we search only surfel of same label (more precise comparisons in some cases)


  CanonicCellEmbedder<KSpace> embeder(K);
  std::vector<unsigned int> vectIndexMinToReference;

  trace.info() << "Associating each input surfel to reference:" << std::endl;
  trace.progressBar(0, surfelAndScalarInput.size());

  for(unsigned int i=0;i <vectSurfelsInput.size(); i++){
    trace.progressBar(i, vectSurfelsInput.size());
    Z3i::RealPoint ptCenterInput = embeder(vectSurfelsInput.at(i));
    unsigned int indexDistanceMin = 0;
    double distanceMin = std::numeric_limits<int>::max() ;
    for(unsigned int j=0; j <vectSurfelsReference.size(); j++){
      if(useLabels && vectLabelsReference.at(j) != vectLabelsInput.at(i)){
        continue;
      }
      Z3i::RealPoint ptCenterRef = embeder(vectSurfelsReference.at(j));
      double distance = (ptCenterRef - ptCenterInput).norm();
      if(distance < distanceMin){
        distanceMin = distance;
        indexDistanceMin = j;
      }
    }
    vectIndexMinToReference.push_back(indexDistanceMin);

  }
  trace.progressBar(surfelAndScalarInput.size(), surfelAndScalarInput.size());
  trace.info() << std::endl;



  //-------------------------
  // Displaying input with error and computing statistics


  QApplication application(argc,argv);
  typedef ViewerSnap<> Viewer;


  Viewer viewer(K, snapShotName!="");
  if(snapShotName!=""){
    viewer.setSnapshotFileName(snapShotName.c_str());
  }
  viewer.setWindowTitle("3dCompSurfel Viewer");
  viewer.show();

  Statistic<double> statErrors(true);
  std::ofstream outputStatStream;
  
  if(fileMeasureOutput != ""){
    outputStatStream.open(fileMeasureOutput.c_str(), ios::app );
  }



  double maxSqError=0;
  for(unsigned int i=0;i <surfelAndScalarInput.size(); i++){
    double scalarInput = surfelAndScalarInput.at(i)[3];
    double scalarRef = surfelAndScalarReference.at(vectIndexMinToReference.at(i))[3];
    double sqError = (scalarRef-scalarInput)*(scalarRef-scalarInput);
    statErrors.addValue(sqError);
    if(sqError> maxSqError){
      maxSqError =sqError;
    }
  }
  if (maxValSpe->count() == 0 )
  {
    maxVal = maxSqError;
  }
  

  GradientColorMap<double> gradientColorMap( 0, maxVal );
  gradientColorMap.addColor( Color(255,255,255,100 ));
  gradientColorMap.addColor( Color(255,0,0,100 ) );
  gradientColorMap.addColor( Color(0,0,255,100 ) );


  //trace.info() << "Maximal error:" << maxSqError << std::endl;
  // Hack waiting issue #899 if maxSqError =0, don't use gradientColorMap
  //bool useGrad = maxSqError!=0.0;

  viewer << SetMode3D(vectSurfelsInput.at(0).className(), "Basic");
  for(unsigned int i=0; i <surfelAndScalarInput.size(); i++){
    double scalarInput = surfelAndScalarInput.at(i)[3];
    double scalarRef = surfelAndScalarReference.at(vectIndexMinToReference.at(i))[3];
    double sqError = (scalarRef-scalarInput)*(scalarRef-scalarInput);
    viewer.setFillColor(gradientColorMap(sqError));
    
    viewer << vectSurfelsInput.at(i);
    if(drawSurfelAssociations){
      viewer.addLine(embeder(vectSurfelsInput.at(i)),embeder(vectSurfelsReference.at(vectIndexMinToReference.at(i))));
    }
  }

  statErrors.terminate();
  if(fileMeasureOutput != ""){
    outputStatStream << "input= " <<  inputFileName << " reference=" << referenceFileName << " " ;
    outputStatStream << statErrors << std::endl;
  }else{
    trace.info()  << statErrors;
  }
  viewer << Viewer::updateDisplay;
  
  if(snapShotName != ""){
    // Appy cleaning just save the last snap
    viewer.restoreStateFromFile();
    std::string extension = snapShotName.substr(snapShotName.find_last_of(".") + 1);
    std::string basename = snapShotName.substr(0, snapShotName.find_last_of("."));
    for(int i=0; i< viewer.snapshotCounter()-1; i++){
      std::stringstream s;
      s << basename << "-"<< setfill('0') << setw(4)<<  i << "." << extension;
      trace.info() << "erase temp file: " << s.str() << std::endl;
      remove(s.str().c_str());
    }
    std::stringstream s;
    s << basename << "-"<< setfill('0') << setw(4)<<  viewer.snapshotCounter()-1 << "." << extension;
    rename(s.str().c_str(), snapShotName.c_str());
    return 0;
  }

  if(noWindows){
    return 0;
  }else{
    return application.exec();
  }
}
