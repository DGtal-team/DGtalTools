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
 * @file volSegment.cpp
 * @ingroup volumetric
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2014/07/20
 *
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include <DGtal/topology/SetOfSurfels.h>
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;


/**
 @page volSegment volSegment
 
 @brief Segments volumetric file from a simple threshold which can be set automatically from the otsu estimation.


 @b Usage: volSegment [input] [output]


 @b Allowed @b options @b are : 
 @code

 Positionals:
   1 TEXT:FILE REQUIRED                  volumetric input file (.vol, .pgm, .pgm3d, .longvol).

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         volumetric input file (.vol, .pgm, .pgm3d, .longvol).
   -o,--output TEXT=result.vol           volumetric output file (.vol, .pgm, .pgm3d, .longvol)
   --labelBackground                     option to define a label to regions associated to object background.
   -m,--thresholdMin INT=0               min threshold (if not given the max threshold is computed with Otsu algorithm).
   -M,--thresholdMax INT=255             max threshold
   @endcode

 @b Example: 
 
 You can test the segmentation in the lobster volume file:

 @code
 $ volSegment -i ${DGtal}/examples/samples/lobster.vol  -o segmentation.vol  -m 70 -M 255
 @endcode

You will obtain a volumetric file representing for each voxel a label associated to a connected component. You can display this segmentation results by extracting it in SDP format with the @ref vol2sdp tool (with option -e to export also the image labels): 
@code
 $ vol2sdp -i segmentation.vol -o segmentation.sdp  -e -m 1 -M 255
@endcode

and display them with @ref Doc3DSDPViewer :

 @code
 $ 3dSDPViewer -i segmentation.sdp  --importColorLabels
 @endcode


 You should obtain such a result:
 @image html resVolSegment.png "Segmentation result displayed with colors representing the segmentation labels."
 
  @see
  @ref volSegment.cpp
 */

typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;

std::vector<unsigned int> getHistoFromImage(const Image3D &image){
  const Image3D::Domain &imgDom = image.domain();
  std::vector<unsigned int> vectHisto(UCHAR_MAX);
  for(Image3D::Domain::ConstIterator it=imgDom.begin(); it!= imgDom.end(); ++it){
    vectHisto[image(*it)]++;
  }
  return vectHisto;
}


unsigned int 
getOtsuThreshold(const Image3D &image){
  std::vector<unsigned int> histo = getHistoFromImage(image);
  unsigned int imageSize = image.domain().size();
  unsigned int sumA = 0;
  unsigned int sumB = imageSize;
  unsigned int muA=0;
  unsigned int muB=0;
  unsigned int sumMuAll= 0;
  for( unsigned int t=0; t< histo.size();t++){
    sumMuAll+=histo[t]*t;
  }
  
  unsigned int thresholdRes=0;
  double valMax=0.0;
  for( unsigned int t=0; t< histo.size(); t++){
    sumA+=histo[t];
    if(sumA==0)
      continue; 
    sumB=imageSize-sumA;
    if(sumB==0){
      break;
    }
    
    muA+=histo[t]*t;
    muB=sumMuAll-muA;
    double muAr=muA/(double)sumA;
    double muBr=muB/(double)sumB;
    double sigma=  (double)sumA*(double)sumB*(muAr-muBr)*(muAr-muBr);
    if(valMax<=sigma){
      valMax=sigma;
      thresholdRes=t;
    }
  }
  return thresholdRes;
}



int main( int argc, char** argv )
{
  
  typedef Z3i::KSpace::SurfelSet SurfelSet;
  typedef SetOfSurfels< Z3i::KSpace, SurfelSet > MySetOfSurfels;
  
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string outputFileName {"result.vol"};
  bool labelBackground {false};
  int minTh {0};
  int maxTh {255};
  
  app.description("Segment volumetric file from a simple threshold which can be set automatically from the otsu estimation.\n The segmentation result is given by an integer label given in the resulting image. Example:\n volSegment -i ${DGtal}/examples/samples/lobster.vol -o segmentation.vol \n");
  
  
  app.add_option("-i,--input,1", inputFileName, "volumetric input file (.vol, .pgm, .pgm3d, .longvol)." )
  ->required()
  ->check(CLI::ExistingFile);

  app.add_option("-o,--output", outputFileName, "volumetric output file (.vol, .pgm, .pgm3d, .longvol)", true);
  auto labelOpt = app.add_flag("--labelBackground",labelBackground, "option to define a label to regions associated to object background.");
  app.add_option("-m,--thresholdMin",minTh, "min threshold (if not given the max threshold is computed with Otsu algorithm).", true );
  auto maxThOpt = app.add_option("-M,--thresholdMax", maxTh, "max threshold", true );

  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  trace.info() << "Reading input file " << inputFileName ;
  Image3D inputImage = DGtal::GenericReader<Image3D>::import(inputFileName);
  Image3D imageResuSegmentation(inputImage.domain());
  
  trace.info() << " [done] " << std::endl ; 
  std::ofstream outStream;
  outStream.open(outputFileName.c_str());
  if(maxThOpt->count()==0){
    maxTh = getOtsuThreshold(inputImage);
    trace.info() << "maximal threshold value not specified, using Otsu value: "  << maxTh << std::endl;
  }
  trace.info() << "Processing image to output file " << outputFileName << std::endl;

  functors::IntervalForegroundPredicate<Image3D> simplePredicate ( inputImage, minTh, maxTh );
  SurfelAdjacency< Z3i::KSpace::dimension > SAdj ( true );
  Z3i::KSpace K;
  bool space_ok = K.init( inputImage.domain().lowerBound(), inputImage.domain().upperBound(), false );
  if(!space_ok){
     trace.error() << "problem initializing 3d space" << endl;
  }

  std::vector< std::vector<Z3i::SCell > > vectConnectedSCell;
  Surfaces<Z3i::KSpace>::extractAllConnectedSCell(vectConnectedSCell,K, SAdj, simplePredicate, false);
  trace.progressBar(0, vectConnectedSCell.size());
  for(unsigned int i = 0; i<vectConnectedSCell.size(); i++)
    {
      trace.progressBar(i, vectConnectedSCell.size());
      MySetOfSurfels  aSet(K, SAdj);
      Z3i::Point lowerPoint, upperPoint;
      Z3i::Point p1;
      Z3i::Point p2;
      for(std::vector<Z3i::SCell>::const_iterator it= vectConnectedSCell.at(i).begin(); it != vectConnectedSCell.at(i).end(); ++it)
        {
          aSet.surfelSet().insert(aSet.surfelSet().begin(),  *it);
          unsigned int orth_dir = K.sOrthDir( *it );                                                     
          p1 =  K.sCoords( K.sIncident( *it, orth_dir, true ) );               
          p2 =  K.sCoords( K.sIncident( *it, orth_dir, false ) );
          if(p1[0] < lowerPoint[0]) lowerPoint[0]= p1[0];
          if(p1[1] < lowerPoint[1]) lowerPoint[1]= p1[1];
          if(p1[2] < lowerPoint[2]) lowerPoint[2]= p1[2];

          if(p1[0] > upperPoint[0]) upperPoint[0]= p1[0];
          if(p1[1] > upperPoint[1]) upperPoint[1]= p1[1];
          if(p1[2] > upperPoint[2]) upperPoint[2]= p1[2];

          if(p2[0] < lowerPoint[0]) lowerPoint[0]= p2[0];
          if(p2[1] < lowerPoint[1]) lowerPoint[1]= p2[1];
          if(p2[2] < lowerPoint[2]) lowerPoint[2]= p2[2];

          if(p2[0] > upperPoint[0]) upperPoint[0]= p2[0];
          if(p2[1] > upperPoint[1]) upperPoint[1]= p2[1];
          if(p2[2] > upperPoint[2]) upperPoint[2]= p2[2];
          
        }    
       
       Z3i::KSpace kRestr ;
       kRestr.init( lowerPoint, upperPoint, false );
       if(simplePredicate(p2)){
         DGtal::Surfaces<Z3i::KSpace>::uFillInterior( kRestr,  aSet.surfelPredicate(), 
                                                      imageResuSegmentation,
                                                      i, false, false);
       }else if (labelOpt->count() > 0){
         DGtal::Surfaces<Z3i::KSpace>::uFillExterior( kRestr,  aSet.surfelPredicate(), 
                                                      imageResuSegmentation,
                                                      i+1, false, false);
       }
    }
  trace.progressBar(vectConnectedSCell.size(), vectConnectedSCell.size());
  trace.info() << std::endl;
  GenericWriter<Image3D>::exportFile(outputFileName, imageResuSegmentation);
  return 0;
}

