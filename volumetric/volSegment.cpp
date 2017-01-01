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


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;


/**
 @page volSegment volSegment
 
 @brief Segments volumetric file from a simple threshold which can be set automatically from the otsu estimation.


 @b Usage: volSegment [input] [output]


 @b Allowed @b options @b are : 
 @code
  -h [ --help ]                  display this message
  -i [ --input ] arg             volumetric input file (.vol, .pgm, .pgm3d, 
                                 .longvol) 
  -o [ --output ] arg            volumetric output file (.vol, .pgm, .pgm3d, 
                                 .longvol) 
  --labelBackground              option to define a label to regions associated
                                 to object background. 
  -m [ --thresholdMin ] arg (=0) min threshold (if not given the max threshold 
                                 is computed with Otsu algorithm).
  -M [ --thresholdMax ] arg      max threshold (default 255)
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
  
  // parse command line ----------------------------------------------
  po::options_description general_opt(" \n Allowed options are ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "volumetric input file (.vol, .pgm, .pgm3d, .longvol) " )
    ("output,o", po::value<std::string>(), "volumetric output file (.vol, .pgm, .pgm3d, .longvol) " )
    ("labelBackground", "option to define a label to regions associated to object background. ")
    ("thresholdMin,m", po::value<int>()->default_value(0), "min threshold (if not given the max threshold is computed with Otsu algorithm)." )
    ("thresholdMax,M", po::value<int>(), "max threshold (default 255)" );
  
  
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
      std::cout << "Usage: " << argv[0] << " [input] [output]\n"
		<< "Segment volumetric file from a simple threshold which can be set automatically from the otsu estimation.\n"
                << "The segmentation result is given by an integer label given in the resulting image."
		<< general_opt << "\n";
      std::cout << "Example:\n"
		<< "volSegment -i ${DGtal}/examples/samples/lobster.vol -o segmentation.vol \n";
      return 0;
    }
  
  
  if(! vm.count("input") ||! vm.count("output"))
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;      
      return 0;
    }

  
  string inputFilename = vm["input"].as<std::string>();
  string outputFilename = vm["output"].as<std::string>();
  
  trace.info() << "Reading input file " << inputFilename ; 
  Image3D inputImage = DGtal::GenericReader<Image3D>::import(inputFilename);
  Image3D imageResuSegmentation(inputImage.domain());
  
  trace.info() << " [done] " << std::endl ; 
  std::ofstream outStream;
  outStream.open(outputFilename.c_str());
  int minTh = vm["thresholdMin"].as<int>();
  int maxTh = 128;
  if(!vm.count("thresholdMax")){
    maxTh = getOtsuThreshold(inputImage);
    trace.info() << "maximal threshold value not specified, using Otsu value: "  << maxTh << std::endl;
  }else{
   maxTh =  vm["thresholdMax"].as<int>();
  }  
  trace.info() << "Processing image to output file " << outputFilename << std::endl; 

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
       }else if (vm.count("labelBackground")){
         DGtal::Surfaces<Z3i::KSpace>::uFillExterior( kRestr,  aSet.surfelPredicate(), 
                                                      imageResuSegmentation,
                                                      i+1, false, false);
       }
    }
  trace.progressBar(vectConnectedSCell.size(), vectConnectedSCell.size());
  trace.info() << std::endl;
  GenericWriter<Image3D>::exportFile(outputFilename, imageResuSegmentation);   
  return 0;
}

