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
 * DGtal convert grey scales image to fremann contour. 
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

//boost
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

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

using namespace DGtal;



/**
 @page displayContours displayContours
 
 @brief Displays discrete contours. 


 @b Usage:   	 displayContours [options] -i  <fileName>  


 @b Allowed @b options @b are :
 
 @code
  -h [ --help ]           display this message
  -i [ --input ] arg      input Freeman chain file name
  --SDP arg               Import a contour as a Sequence of Discrete Points 
                          (SDP format)
  --SFP arg               Import a contour as a Sequence of Floating Points 
                          (SFP format)
  --drawContourPoint arg  <size> display contour points as disk of radius 
                          <size>
  --fillContour           fill the contours with default color (gray)
  --lineWidth arg (=1)    Define the linewidth of the contour (SDP format)
  --drawPointOfIndex arg  <index> Draw the contour point of index <index> 
                          (default 0) 
  --pointSize arg (=2)    <size> Set the display point size of the point 
                          displayed by drawPointofIndex option (default 2.0) 
  --noXFIGHeader           to exclude xfig header in the resulting output 
                          stream (no effect with option -outputFile).
  --withProcessing arg    Processing (used only when the input is a Freeman chain (--input)):
                           DSS segmentation {DSS}
                            Maximal segments {MS}
                           Faithful Polygon {FP}
                           Minimum Length Polygon {MLP}
  -o [ --outputFile ] arg  <filename> save output file automatically according 
                          the file format extension.
  -v [ --displayVectorField ] arg    Add the display of a vector field 
                                     represented by two floating coordinates. 
                                     Each vector is displayed starting from the
                                     corresponding contour point coordinates.
  -v [ --scaleVectorField ] arg (=1) set the scale of the vector field (default
                                     1) (used with --displayVectorField).
  --vectorFieldIndex arg             specify the vector field index (by default
                                     0,1) (used with --displayVectorField).
  --vectorFromAngle arg              specify that the vectors are defined from 
                                     an angle value represented at the given 
                                     index  (by default 0) (used with 
                                     --displayVectorField).

  --rotateVectorField                apply a CCW rotation of 90° (used with 
                                     --displayVectorField).  

  --outputStreamEPS        specify eps for output stream format.
  --outputStreamSVG        specify svg for output stream format.
  --outputStreamFIG        specify fig for output stream format.
  --invertYaxis            invertYaxis invert the Y axis for display contours 
                          (used only with --SDP)
  --backgroundImage arg   backgroundImage <filename> : display image as 
                          background 
  --alphaBG arg           alphaBG <value> 0-1.0 to display the background image
                          in transparency (default 1.0), (transparency works 
                          only if cairo is available)
  --scale arg             scale <value> 1: normal; >1 : larger ; <1 lower 
                          resolutions  )
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



///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "input FreemanChain file name")
    ("SDP", po::value<std::string>(), "Import a contour as a Sequence of Discrete Points (SDP format)")
    ("SFP", po::value<std::string>(), "Import a contour as a Sequence of Floating Points (SFP format)")
    ("drawContourPoint", po::value<double>(), "<size> display contour points as disk of radius <size>")    
    ("fillContour", "fill the contours with default color (gray)")
    ("lineWidth", po::value<double>()->default_value(1.0), "Define the linewidth of the contour (SDP format)") 
    ("drawPointOfIndex", po::value<int>(), "<index> Draw the contour point of index <index> (default 0) ") 
    ("pointSize", po::value<double>()->default_value(2.0), "<size> Set the display point size of the point displayed by drawPointofIndex option (default 2.0) ") 
    ("noXFIGHeader", " to exclude xfig header in the resulting output stream (no effect with option -outputFile).")
    ("withProcessing", po::value<std::string>(), "Processing (used only when the input is a Freeman chain (--input)):\n\t DSS segmentation {DSS}\n\t  Maximal segments {MS}\n\t Faithful Polygon {FP}\n\t Minimum Length Polygon {MLP}")   
    ("outputFile,o", po::value<std::string>(), " <filename> save output file automatically according the file format extension.")
    ("displayVectorField,v", po::value<std::string>(), "Add the display of a vector field represented by two floating coordinates. Each vector is displayed starting from the corresponding contour point coordinates.")
    ("scaleVectorField,v", po::value<double>()->default_value(1.0), "set the scale of the vector field (default 1) (used with --displayVectorField).")
    ("vectorFieldIndex", po::value<std::vector <unsigned int> >()->multitoken(), "specify the vector field index (by default 0,1) (used with --displayVectorField)." )
    ("vectorFromAngle", po::value<unsigned int>(), "specify that the vectors are defined from an angle value represented at the given index  (by default 0) (used with --displayVectorField)." )
    ("rotateVectorField", "apply a CCW rotation of 90° (used with --displayVectorField).  ") 
    ("outputStreamEPS", " specify eps for output stream format.")
    ("outputStreamSVG", " specify svg for output stream format.")
    ("outputStreamFIG", " specify fig for output stream format.")
    ("invertYaxis", " invertYaxis invert the Y axis for display contours (used only with --SDP)")

    ("backgroundImage", po::value<std::string>(), "backgroundImage <filename> : display image as background ")
    ("alphaBG", po::value<double>(), "alphaBG <value> 0-1.0 to display the background image in transparency (default 1.0), (transparency works only if cairo is available)")

    ("scale", po::value<double>(), "scale <value> 1: normal; >1 : larger ; <1 lower resolutions  )");
  
  
 
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }

  po::notify(vm);    
  if(!parseOK||vm.count("help")||argc<=1 || (!(vm.count("input")) && !(vm.count("SDP")) && !(vm.count("SFP"))&&
					     !(vm.count("backgroundImage")) ) )
    {
      trace.info()<< "Display discrete contours. " <<std::endl << "Basic usage: "<<std::endl
		  << "\t displayContours [options] --input  <fileName>  "<<std::endl
		  << general_opt << "\n";
      return 0;
    }
  
  
  
  double lineWidth=  vm["lineWidth"].as<double>();
  bool filled = vm.count("fillContour");
  double scale=1.0;
  if(vm.count("scale")){
    scale = vm["scale"].as<double>();
  }
  
  Board2D aBoard;
  aBoard.setUnit (0.05*scale, LibBoard::Board::UCentimeter);
  




  double alpha=1.0;
  if(vm.count("alphaBG")){
    alpha = vm["alphaBG"].as<double>(); 
  }
  
  if(vm.count("backgroundImage")){
    std::string imageName = vm["backgroundImage"].as<std::string>();
    typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image;
    Image img = DGtal::GenericReader<Image>::import( imageName );
    Z2i::Point ptInf = img.domain().lowerBound(); 
    Z2i::Point ptSup = img.domain().upperBound(); 
    unsigned int width = abs(ptSup[0]-ptInf[0]+1);
    unsigned int height = abs(ptSup[1]-ptInf[1]+1);
    
    aBoard.drawImage(imageName, 0-0.5,height-0.5, width, height, -1, alpha );
  }

 

 
  if(vm.count("input")){
    std::string fileName = vm["input"].as<std::string>();
    std::vector< FreemanChain<int> > vectFc =  PointListReader< Z2i::Point>:: getFreemanChainsFromFile<int> (fileName); 
    aBoard << CustomStyle( vectFc.at(0).className(), 
			   new CustomColors( Color::Red  ,  filled?  Color::Gray: Color::None  ) );    
    aBoard.setLineWidth (lineWidth);
    for(unsigned int i=0; i<vectFc.size(); i++){  
      aBoard <<  vectFc.at(i) ;
      if(vm.count("drawPointOfIndex")){
	int index = vm["drawPointOfIndex"].as<int>();
	double size = vm["pointSize"].as<double>();
	aBoard.setPenColor(Color::Blue);
	     
	aBoard.fillCircle((double)(vectFc.at(i).getPoint(index)[0]), (double)(vectFc.at(i).getPoint(index)[1]), size);
      }

      if(vm.count("withProcessing")){
	std::string processingName = vm["withProcessing"].as<std::string>();

	std::vector<Z2i::Point> vPts(vectFc.at(i).size()+1); 
	copy ( vectFc.at(i).begin(), vectFc.at(i).end(), vPts.begin() ); 
	bool isClosed;
	if ( vPts.at(0) == vPts.at(vPts.size()-1) ) { 
          isClosed = true;
          vPts.pop_back(); 
	} else isClosed = false;

	if (processingName == "DSS") {

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

	} else if (processingName == "FP") {

	  typedef FP<std::vector<Z2i::Point>::iterator,int,4> FP;
	  FP theFP( vPts.begin(),vPts.end() );
          aBoard << CustomStyle( theFP.className(), 
				 new CustomPenColor( DGtal::Color::Black ) ); 
          aBoard << theFP;


	} else if (processingName == "MLP") {

	  typedef FP<std::vector<Z2i::Point>::iterator,int,4> FP;
	  FP theFP( vPts.begin(),vPts.end() );

          std::vector<FP::RealPoint> v( theFP.size() );
          theFP.copyMLP( v.begin() );

          //polyline to draw
	  std::vector<LibBoard::Point> polyline;
	  std::vector<FP::RealPoint>::const_iterator it = v.begin();
	  for ( ;it != v.end();++it) {
	    FP::RealPoint p = (*it);
	    polyline.push_back(LibBoard::Point(p[0],p[1]));
	  }
          if (isClosed) {
	    FP::RealPoint p = (*v.begin());
	    polyline.push_back(LibBoard::Point(p[0],p[1]));
          }
          aBoard.setPenColor(DGtal::Color::Black);
	  aBoard.drawPolyline(polyline);
	  	  
	} else if (processingName == "MDCA") {
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
	  for ( ; it != itEnd; ++it ) {
	    aBoard << SetMode(SegmentComputer().className(), "") << (*it); 
	    otherBoard << SetMode(SegmentComputer().className(), "") << (*it); 
	  }
	  otherBoard.saveSVG("mdca.svg", Board2D::BoundingBox, 5000 ); 
	}
      }

    }



  }
 
 

  if(vm.count("SDP") || vm.count("SFP")){
    bool drawPoints= vm.count("drawContourPoint");
    bool invertYaxis = vm.count("invertYaxis");
    double pointSize=1.0;
    if(drawPoints){
      pointSize = vm["drawContourPoint"].as<double>();
    }
    std::vector<LibBoard::Point> contourPt;
    if(vm.count("SDP")){
      std::string fileName = vm["SDP"].as<std::string>();
      std::vector< Z2i::Point  >  contour = PointListReader< Z2i::Point >::getPointsFromFile(fileName); 
      for(unsigned int j=0; j<contour.size(); j++){
	LibBoard::Point pt((double)(contour.at(j)[0]),
			   (invertYaxis? (double)(-contour.at(j)[1]+contour.at(0)[1]):(double)(contour.at(j)[1])));
	contourPt.push_back(pt);
	if(drawPoints){
	  aBoard.fillCircle(pt.x, pt.y, pointSize);
	}
      }
    }
 
    if(vm.count("SFP")){
      std::string fileName = vm["SFP"].as<std::string>();
      std::vector< PointVector<2,double>  >  contour = 
	PointListReader<  PointVector<2,double>  >::getPointsFromFile(fileName); 
      for(unsigned int j=0; j<contour.size(); j++){
	LibBoard::Point pt((double)(contour.at(j)[0]),
			   (invertYaxis? (double)(-contour.at(j)[1]+contour.at(0)[1]):(double)(contour.at(j)[1])));
	contourPt.push_back(pt);
	if(drawPoints){
	  aBoard.fillCircle(pt.x, pt.y, pointSize);
	}
      }
      
    }
    
    
    aBoard.setPenColor(Color::Red);
    aBoard.setFillColor(Color::Gray);
    aBoard.setLineStyle (LibBoard::Shape::SolidStyle );
    aBoard.setLineWidth (lineWidth);
    if(!filled){
      aBoard.drawPolyline(contourPt);
    }else{
      aBoard.fillPolyline(contourPt);
    }
    if(vm.count("drawPointOfIndex")){
      int index = vm["drawPointOfIndex"].as<int>();
      double size = vm["pointSize"].as<double>();
      aBoard.fillCircle((double)(contourPt.at(index).x), (double)(contourPt.at(index).y), size);
    }   

    
    // display vector field
    if(vm.count("displayVectorField"))
      {
        bool rotate = vm.count("rotateVectorField");
        double sv = vm["scaleVectorField"].as<double>();
        std::vector<unsigned int> vIndex = {0,1};
        if(vm.count("vectorFieldIndex"))
          {
            vIndex = vm["vectorFieldIndex"].as<std::vector<unsigned int>>();
          }
        std::string vname = vm["displayVectorField"].as<std::string>();
        std::vector< PointVector<2,double>  >  vField;
        if(vm.count("vectorFromAngle"))
          {
            unsigned int aIndex = vm["vectorFromAngle"].as<unsigned int>();
            std::vector<double> vAngles  = TableReader<double>::getColumnElementsFromFile(vname, aIndex); 
            for(unsigned int i = 0; i < vAngles.size(); i++)
              {
                vField.push_back(Z2i::RealPoint(cos(vAngles[i]),sin(vAngles[i])));
              }
          }
        else
          {
            vField = PointListReader<  PointVector<2,double>  >::getPointsFromFile(vname, vIndex);
          }
        for(unsigned int i = 0; i< contourPt.size(); i++)
          {
            vField[i] = vField[i].getNormalized();
            auto p = contourPt[i];
            if(!rotate)
              {
                aBoard.drawArrow(p.x, p.y, p.x+vField[i][0]*sv, p.y+vField[i][1]*sv  );
              }
            else
              {
                aBoard.drawArrow(p.x, p.y, p.x-vField[i][1]*sv, p.y+vField[i][0]*sv  );
              }
              
          }
        
      }
  }



  
  if(vm.count("outputFile")){
    std::string outputFileName= vm["outputFile"].as<std::string>();
    std::string extension = outputFileName.substr(outputFileName.find_last_of(".") + 1);

    if(extension=="svg"){
      aBoard.saveSVG(outputFileName.c_str());
    }
    #ifdef WITH_CAIRO
    else
      if (extension=="eps"){
	aBoard.saveCairo(outputFileName.c_str(),Board2D::CairoEPS );
      } else 
	if (extension=="pdf"){
	  aBoard.saveCairo(outputFileName.c_str(),Board2D::CairoPDF );
	} else 
	  if (extension=="png"){
	    aBoard.saveCairo(outputFileName.c_str(),Board2D::CairoPNG );
	  }
    #endif
    else if(extension=="eps"){
      aBoard.saveEPS(outputFileName.c_str());
    }else if(extension=="fig"){
      aBoard.saveFIG(outputFileName.c_str(),LibBoard::Board::BoundingBox, 10.0, !vm.count("noXFIGHeader") );
    }
  }
    
    if (vm.count("outputStreamSVG")){
    aBoard.saveSVG(std::cout);
  } else   
      if (vm.count("outputStreamFIG")){
    aBoard.saveFIG(std::cout, LibBoard::Board::BoundingBox, 10.0,  !vm.count("noXFIGHeader"));
  } else
	if (vm.count("outputStreamEPS")){
    aBoard.saveEPS(std::cout);
  } 
    
  }

