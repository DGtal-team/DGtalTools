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
 * @file displayContours.cpp
 * @ingroup visualization
 * @author Bertrand Kerautret (\c kerautre@loria.fr)
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2011/27/04
 *
 * Display discrete contours
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

//STL
#include <vector>
#include <string>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/topology/helpers/Surfaces.h"

//image
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/readers/TableReader.h"
#include "DGtal/io/Color.h"

#include "DGtal/io/readers/GenericReader.h"

//contour
#include "DGtal/geometry/curves/FreemanChain.h"

//processing
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/GreedySegmentation.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/FP.h"
#include "DGtal/geometry/curves/StabbingCircleComputer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/SegmentComputerUtils.h"

#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/boards/CDrawableWithBoard2D.h"

#include "CLI11.hpp"

using namespace DGtal;



/**
 @page displayContours displayContours
 
 @brief Displays discrete contours.
 
 
 @b Usage:   	 displayContours [options] -i  <fileName>
 
 
 @b Allowed @b options @b are :
 
 @code

 Positionals:
   1 TEXT:FILE                           input FreemanChain file name

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE                  input FreemanChain file name
   -o,--outputFile TEXT                  save output file automatically according the file format extension.
   --SDP TEXT:FILE                       Import a contour as a Sequence of Discrete Points (SDP format)
   --SFP TEXT:FILE                       mport a contour as a Sequence of Floating Points (SFP format)
   --drawContourPoint FLOAT              <size> display contour points as disk of radius <size>
   --fillContour                         fill the contours with default color (gray)
   --lineWidth FLOAT                     Define the linewidth of the contour (SDP format)
   -f,--drawPointOfIndex UINT             Draw the contour point of index.
   --pointSize FLOAT                     <size> Set the display point size of the point displayed by drawPointofIndex option (default 2.0)
   --noXFIGHeader BOOLEAN                to exclude xfig header in the resulting output stream (no effect with option -outputFile).
   --withProcessing TEXT:{MS,FP,MLP}     Processing (used only when the input is a Freeman chain (--input)):
                                            DSS segmentation {DSS}
                                             Maximal segments {MS}
                                            Faithful Polygon {FP}
                                            Minimum Length Polygon {MLP}
   -v,--displayVectorField TEXT          Add the display of a vector field represented by two floating coordinates. Each vector is displayed starting from the corresponding contour point coordinates.
   --scaleVectorField FLOAT=1            set the scale of the vector field (default 1) (used with --displayVectorField).
   --vectorFieldIndex UINT=[0,1] x 2     specify the vector field index (by default 0,1) (used with --displayVectorField).
   --vectorFromAngle UINT                specify that the vectors are defined from an angle value represented at the given index  (by default 0) (used with --displayVectorField).
   --rotateVectorField                   apply a CCW rotation of 90° (used with --displayVectorField).
   --outputStreamEPS                      specify eps for output stream format.
   --outputStreamSVG                      specify svg for output stream format.
   --outputStreamFIG                      specify fig for output stream format.
   --invertYaxis Needs: --SDP             invertYaxis invert the Y axis for display contours (used only with --SDP)
   --backgroundImage TEXT:FILE           backgroundImage <filename> : display image as background
   --alphaBG FLOAT=1                     alphaBG <value> 0-1.0 to display the background image in transparency (default 1.0), (transparency works only if cairo is available)
   --scale FLOAT=1                       scale <value> 1: normal; >1 : larger ; <1 lower resolutions)
   
 @endcode
 
 
 @b Example:
 
 In this example we show how to display of a set of contours extracted
 in a single image. The first step is to extract a set contours by
 using the tool
 
 @code
 $ img2freeman -i $DGtal/examples/samples/church.pgm  -R 0 20 255 -s 200  > church.fc
 @endcode
 
 Then, we display the set of contours with the background images (you need to have compiled DGTal with the options (-DWITH_MAGICK=true and -DWITH_CAIRO=true) :
 
 @code
 $  displayContours -i church.fc   --backgroundImage $DGtal/examples/samples/church.png --alphaBG 0.75 --outputFile church.pdf
 @endcode
 
 
 You should obtain such a result:
 
 @image html resDisplayContours.png "Resulting visualization."
 
 
 @see
 @ref displayContours.cpp
 
 */


int main( int argc, char** argv )
{
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  app.description("Display discrete contours. \n Basic example: \t displayContours [options] --input  <fileName>");
  std::string inputFileName;
  std::string outputFileName {"result.svg"};
  std::string inputSDFileName;
  std::string inputSFPFileName;
  std::string processingName;
  unsigned int indexPoint;
  double pointSize {0.0};
  double scaleVectorField {1.0};
  bool fillContour {false};
  bool rotateVectorField {false};
  bool outputStreamEPS {false};
  bool outputStreamSVG {false};
  bool outputStreamFIG {false};
  bool noXFIGHeader {false};
  bool invertYaxis {false};
  double alphaBG {1.0};
  double lineWidth {0.0};
  double scale {1.0};
  unsigned int vectorFromAngle {0};
  std::string displayVectorField = "";
  std::vector <unsigned int> vectorFieldIndex = {0,1};
  std::string backgroundImage;
  
  app.add_option("-i,--input,1", inputFileName, "input FreemanChain file name" )
  ->check(CLI::ExistingFile);
  app.add_option("--outputFile,-o", outputFileName, "save output file automatically according the file format extension.");
  
  CLI::Option* optSDP = app.add_option("--SDP",inputSDFileName, "Import a contour as a Sequence of Discrete Points (SDP format)")
  ->check(CLI::ExistingFile);
  app.add_option("--SFP",inputSFPFileName, "mport a contour as a Sequence of Floating Points (SFP format)")
  ->check(CLI::ExistingFile);
  app.add_option("--drawContourPoint", pointSize, "<size> display contour points as disk of radius <size>");
  app.add_flag("--fillContour", fillContour, "fill the contours with default color (gray)");
  app.add_option("--lineWidth", lineWidth, "Define the linewidth of the contour (SDP format)" );
  auto indexPointOpt =  app.add_option("--drawPointOfIndex,-f", indexPoint, " Draw the contour point of index." );
  app.add_option("--pointSize", pointSize, "<size> Set the display point size of the point displayed by drawPointofIndex option (default 2.0) " );
  app.add_option("--noXFIGHeader", noXFIGHeader, "to exclude xfig header in the resulting output stream (no effect with option -outputFile).");
  
  app.add_option("--withProcessing",processingName, "Processing (used only when the input is a Freeman chain (--input)):\n\t DSS segmentation {DSS}\n\t  Maximal segments {MS}\n\t Faithful Polygon {FP}\n\t Minimum Length Polygon {MLP}" )
  -> check(CLI::IsMember({"MS", "FP", "MLP"}));
  
  app.add_option("--displayVectorField,-v", displayVectorField, "Add the display of a vector field represented by two floating coordinates. Each vector is displayed starting from the corresponding contour point coordinates.");
  app.add_option("--scaleVectorField",scaleVectorField, "set the scale of the vector field (default 1) (used with --displayVectorField).", true);
  app.add_option("--vectorFieldIndex", vectorFieldIndex ,"specify the vector field index (by default 0,1) (used with --displayVectorField).", true)
  ->expected(2);
  
  auto optVFAngke = app.add_option("--vectorFromAngle", vectorFromAngle, "specify that the vectors are defined from an angle value represented at the given index  (by default 0) (used with --displayVectorField).");
  app.add_flag("--rotateVectorField", rotateVectorField, "apply a CCW rotation of 90° (used with --displayVectorField). ");
  app.add_flag("--outputStreamEPS"," specify eps for output stream format.");
  app.add_flag("--outputStreamSVG"," specify svg for output stream format.");
  app.add_flag("--outputStreamFIG"," specify fig for output stream format.");
  app.add_flag("--invertYaxis", invertYaxis, " invertYaxis invert the Y axis for display contours (used only with --SDP)")
  ->needs(optSDP);
  app.add_option("--backgroundImage", backgroundImage, "backgroundImage <filename> : display image as background ")
  ->check(CLI::ExistingFile);
  app.add_option("--alphaBG", alphaBG, "alphaBG <value> 0-1.0 to display the background image in transparency (default 1.0), (transparency works only if cairo is available)", true);
  app.add_option("--scale", scale, "scale <value> 1: normal; >1 : larger ; <1 lower resolutions)", true);
  
  
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  if (argc == 1 )
  {
    trace.info() << app.get_description() << std::endl;
    trace.error() << "You need at least add one input file using -i, --SDP, --SFP (see --help for more details)" << std::endl;
    return EXIT_FAILURE;
  }
  
  Board2D aBoard;
  aBoard.setUnit (0.05*scale, LibBoard::Board::UCentimeter);
  
  if(backgroundImage != "")
  {
    typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image;
    Image img = DGtal::GenericReader<Image>::import( backgroundImage );
    Z2i::Point ptInf = img.domain().lowerBound();
    Z2i::Point ptSup = img.domain().upperBound();
    unsigned int width = abs(ptSup[0]-ptInf[0]+1);
    unsigned int height = abs(ptSup[1]-ptInf[1]+1);
    aBoard.drawImage(backgroundImage, 0-0.5,height-0.5, width, height, -1, alphaBG);
  }
  
  if(inputFileName != ""){
    std::vector< FreemanChain<int> > vectFc =  PointListReader< Z2i::Point>:: getFreemanChainsFromFile<int> (inputFileName);
    aBoard << CustomStyle( vectFc.at(0).className(),
                          new CustomColors( Color::Red  , fillContour?  Color::Gray: Color::None  ) );
    aBoard.setLineWidth (lineWidth);
    for(unsigned int i=0; i<vectFc.size(); i++){
      aBoard <<  vectFc.at(i) ;
      if( indexPointOpt->count() != 0 ){
        aBoard.setPenColor(Color::Blue);
        aBoard.fillCircle((double)(vectFc.at(i).getPoint(indexPoint)[0]),
                          (double)(vectFc.at(i).getPoint(indexPoint)[1]), pointSize);
      }
      
      if(processingName != ""){
        std::vector<Z2i::Point> vPts(vectFc.at(i).size()+1);
        copy ( vectFc.at(i).begin(), vectFc.at(i).end(), vPts.begin() );
        bool isClosed;
        if ( vPts.at(0) == vPts.at(vPts.size()-1) )
        {
          isClosed = true;
          vPts.pop_back();
        }
        else
          isClosed = false;
        
        if (processingName == "DSS")
        {
          typedef ArithmeticalDSSComputer<std::vector<Z2i::Point>::iterator,int,4> DSS4;
          typedef GreedySegmentation<DSS4> Decomposition4;
          
          DSS4 computer;
          Decomposition4 theDecomposition( vPts.begin(),vPts.end(),computer );
          
          //for each segment
          std::string className;
          for ( Decomposition4::SegmentComputerIterator it = theDecomposition.begin();
               it != theDecomposition.end(); ++it )
          {
            DSS4::Primitive segment(it->primitive());
            
            aBoard << SetMode( segment.className(), "BoundingBox" );
            className = segment.className() + "/BoundingBox";
            aBoard << CustomStyle( className,
                                  new CustomPenColor( DGtal::Color::Gray ) );
            aBoard << segment; // draw each segment
          }
          
        } else if (processingName == "MS") {
          
          typedef ArithmeticalDSSComputer<std::vector<Z2i::Point>::iterator,int,4> DSS4;
          typedef SaturatedSegmentation<DSS4> Decomposition4;
          
          //Segmentation
          DSS4 computer;
          Decomposition4 theDecomposition( vPts.begin(),vPts.end(),computer );
          
          //for each segment
          std::string className;
          for ( Decomposition4::SegmentComputerIterator it = theDecomposition.begin();
               it != theDecomposition.end(); ++it )
          {
            DSS4::Primitive segment(it->primitive());
            aBoard << SetMode( segment.className(), "BoundingBox" );
            className = segment.className() + "/BoundingBox";
            aBoard << CustomStyle( className,
                                  new CustomPenColor( DGtal::Color::Gray ) );
            aBoard << segment; // draw each segment
          }
          
        } else if (processingName == "FP")
        {
          typedef FP<std::vector<Z2i::Point>::iterator,int,4> FP;
          FP theFP( vPts.begin(),vPts.end() );
          aBoard << CustomStyle( theFP.className(),
                                new CustomPenColor( DGtal::Color::Black ) );
          aBoard << theFP;
          
        } else if (processingName == "MLP")
        {
          typedef FP<std::vector<Z2i::Point>::iterator,int,4> FP;
          FP theFP( vPts.begin(),vPts.end() );
          
          std::vector<FP::RealPoint> v( theFP.size() );
          theFP.copyMLP( v.begin() );
          
          //polyline to draw
          std::vector<LibBoard::Point> polyline;
          std::vector<FP::RealPoint>::const_iterator it = v.begin();
          for ( ;it != v.end();++it)
          {
            FP::RealPoint p = (*it);
            polyline.push_back(LibBoard::Point(p[0],p[1]));
          }
          if (isClosed)
          {
            FP::RealPoint p = (*v.begin());
            polyline.push_back(LibBoard::Point(p[0],p[1]));
          }
          aBoard.setPenColor(DGtal::Color::Black);
          aBoard.drawPolyline(polyline);
          
        } else if (processingName == "MDCA")
        {
          typedef KhalimskySpaceND<2,int> KSpace;
          typedef GridCurve<KSpace> Curve;
          Curve curve; //grid curve
          curve.initFromPointsVector( vPts );
          typedef Curve::IncidentPointsRange Range; //range
          Range r = curve.getIncidentPointsRange(); //range
          typedef Range::ConstCirculator ConstCirculator; //iterator
          typedef StabbingCircleComputer<ConstCirculator> SegmentComputer; //segment computer
          //typedef GeometricalDCA<ConstIterator> SegmentComputer; //segment computer
          typedef SaturatedSegmentation<SegmentComputer> Segmentation;
          //Segmentation theSegmentation( r.begin(), r.end(), SegmentComputer() );
          Segmentation theSegmentation( r.c(), r.c(), SegmentComputer() );
          theSegmentation.setMode("Last");
          // board << curve;
          Segmentation::SegmentComputerIterator it = theSegmentation.begin();
          Segmentation::SegmentComputerIterator itEnd = theSegmentation.end();
          Board2D otherBoard;
          otherBoard.setPenColor(DGtal::Color::Black);
          otherBoard << curve;
          for ( ; it != itEnd; ++it )
          {
            aBoard << SetMode(SegmentComputer().className(), "") << (*it);
            otherBoard << SetMode(SegmentComputer().className(), "") << (*it);
          }
          otherBoard.saveSVG("mdca.svg", Board2D::BoundingBox, 5000 );
        }
      }
      
    }
  }
  

  if( inputSDFileName != ""  || inputSFPFileName != "" )
  {
    std::vector<LibBoard::Point> contourPt;
    if( inputSDFileName != "" )
    {
      std::vector< Z2i::Point  >  contour = PointListReader< Z2i::Point >::getPointsFromFile(inputSDFileName);
      for(unsigned int j=0; j<contour.size(); j++)
      {
        LibBoard::Point pt((double)(contour.at(j)[0]),
                           (invertYaxis? (double)(-contour.at(j)[1]+contour.at(0)[1]):(double)(contour.at(j)[1])));
        contourPt.push_back(pt);
        if(pointSize != 0.0)
        {
          aBoard.fillCircle(pt.x, pt.y, pointSize);
        }
      }
    }
    
    if( inputSFPFileName != "" )
    {
      std::vector< PointVector<2,double>  >  contour =
      PointListReader<  PointVector<2,double>  >::getPointsFromFile(inputSFPFileName);
      for(unsigned int j=0; j<contour.size(); j++)
      {
        LibBoard::Point pt((double)(contour.at(j)[0]),
                           (invertYaxis? (double)(-contour.at(j)[1]+contour.at(0)[1]):(double)(contour.at(j)[1])));
        contourPt.push_back(pt);
        if(pointSize != 0.0)
        {
          aBoard.fillCircle(pt.x, pt.y, pointSize);
        }
      }
    }
    
    aBoard.setPenColor(Color::Red);
    aBoard.setFillColor(Color::Gray);
    aBoard.setLineStyle (LibBoard::Shape::SolidStyle );
    aBoard.setLineWidth (lineWidth);
    if(!fillContour)
    {
      aBoard.drawPolyline(contourPt);
    }else
    {
      aBoard.fillPolyline(contourPt);
    }
    if( indexPointOpt->count() != 0 )
    {
      aBoard.fillCircle((double)(contourPt.at(indexPoint).x), (double)(contourPt.at(indexPoint).y), pointSize);
    }
    
    
    // display vector field
    if(displayVectorField != "")
    {
      std::vector< PointVector<2,double>  >  vField;
      if(optVFAngke->count() != 0)
      {
        
        std::vector<double> vAngles  = TableReader<double>::getColumnElementsFromFile(displayVectorField, vectorFromAngle);
        for(unsigned int i = 0; i < vAngles.size(); i++)
        {
          vField.push_back(Z2i::RealPoint(cos(vAngles[i]),sin(vAngles[i])));
        }
      }
      else
      {
        vField = PointListReader<  PointVector<2,double>  >::getPointsFromFile(displayVectorField, vectorFieldIndex);
      }
      for(unsigned int i = 0; i< contourPt.size(); i++)
      {
        vField[i] = vField[i].getNormalized();
        auto p = contourPt[i];
        if(!rotateVectorField)
        {
          aBoard.drawArrow(p.x, p.y, p.x+vField[i][0]*scaleVectorField, p.y+vField[i][1]*scaleVectorField  );
        }
        else
        {
          aBoard.drawArrow(p.x, p.y, p.x-vField[i][1]*scaleVectorField, p.y+vField[i][0]*scaleVectorField  );
        }
      }
    }
  }
  
  
  
  
  if( outputFileName != "" )
  {
    std::string extension = outputFileName.substr(outputFileName.find_last_of(".") + 1);
    if(extension=="svg")
    {
      aBoard.saveSVG(outputFileName.c_str());
    }
#ifdef WITH_CAIRO
    else
      if (extension=="eps")
      {
        aBoard.saveCairo(outputFileName.c_str(),Board2D::CairoEPS );
      }
      else if (extension=="pdf")
      {
          aBoard.saveCairo(outputFileName.c_str(),Board2D::CairoPDF );
      }
      else if (extension=="png")
      {
            aBoard.saveCairo(outputFileName.c_str(),Board2D::CairoPNG );
      }
#endif
      else if(extension=="eps")
      {
         aBoard.saveEPS(outputFileName.c_str());
      }
      else if(extension=="fig")
      {
         aBoard.saveFIG(outputFileName.c_str(),LibBoard::Board::BoundingBox, 10.0, !noXFIGHeader );
      }
  }
  
  if (outputStreamSVG)
  {
    aBoard.saveSVG(std::cout);
  }
  else if (outputStreamFIG)
  {
      aBoard.saveFIG(std::cout, LibBoard::Board::BoundingBox, 10.0,  !noXFIGHeader);
  } else if (outputStreamEPS)
  {
        aBoard.saveEPS(std::cout);
  }
  
}

