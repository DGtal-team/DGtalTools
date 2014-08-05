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
 * @file 3dSDPViewer.cpp
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2014/04/01
 *
 * An simple viewer to display 3d SDP files (sequence of discrete points).
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <QtGui/qapplication.h>

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

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <limits> 

using namespace std;
using namespace DGtal;
using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

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
    
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "input file: sdpa (sequence of discrete points with attribute)" )
    ("reference,r", po::value<std::string>(), "input reference file: sdpa (sequence of discrete points with attribute)" )
    ("SDPindex", po::value<std::vector <unsigned int> >()->multitoken(), "specify the sdp index (by default 0,1,2,3).");


  bool parseOK=true;
  bool cannotStart= false;
  Z3i::KSpace K;

  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.error()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);    
  if(parseOK && ! vm.count("input"))
    {
      trace.error() << " The input file name was not defined" << endl;      
      cannotStart = true;
    }
   
  
  if( !parseOK || cannotStart ||  vm.count("help")||argc<=1)
    {
      trace.info() << "Usage: " << argv[0] << " [input]\n"
		<< "Display sequence of 3d discrete points by using QGLviewer."
		<< general_opt << "\n";
      return 0;
    }

  string inputFilename = vm["input"].as<std::string>();
  string referenceFilename = vm["reference"].as<std::string>();

  std::vector<Point4D> surfelAndScalarInput;
  std::vector<Point4D> surfelAndScalarReference;

  
  if(vm.count("SDPindex")) {
    std::vector<unsigned int > vectIndex = vm["SDPindex"].as<std::vector<unsigned int > >();
    if(vectIndex.size()!=3){
      trace.error() << "you need to specify the three indexes of vertex." << std::endl; 
      return 0;
    }
    surfelAndScalarInput = PointListReader<Point4D>::getPointsFromFile(inputFilename, vectIndex);
    surfelAndScalarReference = PointListReader<Point4D>::getPointsFromFile(referenceFilename, vectIndex);
  }else{
    surfelAndScalarInput = PointListReader<Point4D>::getPointsFromFile(inputFilename);
    surfelAndScalarReference = PointListReader<Point4D>::getPointsFromFile(referenceFilename);
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
      Z3i::RealPoint ptCenterRef = embeder(vectSurfelsReference.at(j));
      double distance = (ptCenterRef - ptCenterInput).norm();
      if(distance < distanceMin){
        distanceMin = distance;
        indexDistanceMin = j;
      }
    }
    vectIndexMinToReference.push_back(indexDistanceMin);
    
  }
  
  
  
  

  //-------------------------
  // Displaying input
  
  QApplication application(argc,argv);
  typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;
   

  Viewer viewer( K );
  viewer.setWindowTitle("3dCompSurfel Viewer");
  viewer.show();
  

  double maxSqError=0;
  for(unsigned int i=0;i <surfelAndScalarInput.size(); i++){
    double curvatureInput = surfelAndScalarInput.at(i)[3];
    double curvatureRef = surfelAndScalarReference.at(vectIndexMinToReference.at(i))[3];
    double sqError = (curvatureRef-curvatureInput)*(curvatureRef-curvatureInput);
    if(sqError> maxSqError){
      maxSqError =sqError;
    }
  }
  GradientColorMap<double> gradientColorMap( 0, maxSqError );
  gradientColorMap.addColor( Color(255,255,255,100 ));
  gradientColorMap.addColor( Color(255,0,0,100 ) );
  gradientColorMap.addColor( Color(0,0,255,100 ) );


  
  viewer << SetMode3D(vectSurfelsInput.at(0).className(), "Basic");
  for(unsigned int i=0;i <surfelAndScalarInput.size(); i++){

    double curvatureInput = surfelAndScalarInput.at(i)[3];
    double curvatureRef = surfelAndScalarReference.at(vectIndexMinToReference.at(i))[3];
    double sqError = (curvatureRef-curvatureInput)*(curvatureRef-curvatureInput);

    viewer.setFillColor(gradientColorMap(sqError));
    
    viewer << vectSurfelsInput.at(i);
    
    //viewer.addLine(embeder(vectSurfelsInput.at(i)),embeder(vectSurfelsReference.at(vectIndexMinToReference.at(i))));
    
  }
  

  
  // vector<Z3i::Cell> vectSurfelInput;
  // if(vm.count("SDPindex")) {
  //   std::vector<unsigned int > vectIndex = vm["SDPindex"].as<std::vector<unsigned int > >();
  //   if(vectIndex.size()!=3){
  //     trace.error() << "you need to specify the three indexes of vertex." << std::endl; 
  //     return 0;
  //   }
  //   vectVoxels = PointListReader<Z3i::RealPoint>::getPointsFromFile(inputFilename, vectIndex);
  // }else{
  //   vectVoxels = PointListReader<Z3i::RealPoint>::getPointsFromFile(inputFilename);
  // }
  
  // if(!vm.count("noPointDisplay")){
  //   double percent = vm["filter"].as<double>();
  //   int step = max(1, (int) (100/percent));
  //   for(int i=0;i< vectVoxels.size(); i=i+step){
  //     if(typePrimitive=="voxel"){
  //       viewer << Z3i::Point((int)vectVoxels.at(i)[0],
  //                            (int)vectVoxels.at(i)[1],
  //                            (int)vectVoxels.at(i)[2]); 
  //     }else{
  //       viewer.addBall(vectVoxels.at(i), sphereRadius);    
  //     }
  //   }  
  // }
  
  // viewer << CustomColors3D(lineColor, lineColor);
  // if(vm.count("drawLines")){
  //   for(int i=1;i< vectVoxels.size(); i++){
  //     viewer.addLine(vectVoxels.at(i-1), vectVoxels.at(i), lineSize); 
  //   }  
  // }

  
  // if(vm.count("drawVectors")){
  //   std::string vectorsFileName = vm["drawVectors"].as<std::string>(); 
  //   std::vector<Z3i::RealPoint> vectorsPt = PointListReader<Z3i::RealPoint>::getPointsFromFile(vectorsFileName);
  //   if (vectorsPt.size()%2==1){
  //     trace.info()<<"Warning the two set of points doesn't contains the same number of points, some vectors will be skipped." << std::endl;
  //   }
  //   for(unsigned int i =0; i<vectorsPt.size()-1; i=i+2){
  //     viewer.addLine(vectorsPt.at(i),vectorsPt.at(i+1));      
  //   } 
 
    
  // }
  
  
  viewer << Viewer3D<>::updateDisplay;
  return application.exec();

}
