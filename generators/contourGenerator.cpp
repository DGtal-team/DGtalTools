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
 * @file contourGenerator.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr)
 * LIRIS (CNRS, UMR 5205),
 *
 * @date 2011/07/04
 *
 * DGtal shape generator
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "CLI11.hpp"

#include "DGtal/base/Common.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"

#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"


#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/boards/Board2D.h"

#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/geometry/curves/GridCurve.h"
#include "DGtal/geometry/curves/estimation/TrueLocalEstimatorOnPoints.h"
#include "DGtal/geometry/curves/estimation/TrueGlobalEstimatorOnPoints.h"
#include "DGtal/geometry/curves/estimation/ParametricShapeCurvatureFunctor.h"
#include "DGtal/geometry/curves/estimation/ParametricShapeTangentFunctor.h"
#include "DGtal/geometry/curves/estimation/ParametricShapeArcLengthFunctor.h"

#include "DGtal/topology/helpers/Surfaces.h"


using namespace DGtal;

/**
 *
 @page contourGenerator contourGenerator
 
 @brief  Generates multigrid contours of 2d digital shapes using DGtal library.
 
 The associated geometric information like contour, curvature can also be displayed. 

 @b Usage: 	contourGenerator --shape <shapeName> [requiredParam] [otherOptions]
 
 @b Allowed @b options @b are:
 
 @code
  -h,--help                             Print this help message and exit
  -l,--list                             List all available shapes
  -s,--shape TEXT                       Shape name
  -R,--radius FLOAT                     Radius of the shape
  -A,--axis1 FLOAT                      Half big axis of the shape (ellipse)
  -a,--axis2 FLOAT                      Half small axis of the shape (ellipse)
  -r,--smallradius FLOAT=5              Small radius of the shape (default 5)
  -v,--varsmallradius FLOAT=5           Variable small radius of the shape (default 5)
  -k UINT=3                             Number of branches or corners the shape (default 3)
  --phi FLOAT=0                         Phase of the shape (in radian, default 0.0)
  -w,--width FLOAT=10                   Width of the shape (default 10.0)
  -p,--power FLOAT=2                    Power of the metric (default 2.0)
  -x,--center_x FLOAT=0                 x-coordinate of the shape center (default 0.0)
  -y,--center_y FLOAT=0                 y-coordinate of the shape center (default 0.0)
  -g,--gridstep FLOAT=1                 Gridstep for the digitization (default 1.0)
  -f,--format TEXT=pts                  Output format:
                                            List of pointel coordinates {pts}
                                            Freeman chaincode Vector {fc} (default pts)
  -o,--outputGeometry TEXT              Base name of the file containing the shape geometry (points, tangents, curvature)
 @endcode

 You can also list all possible shapes:
 @code 
 contourGenerator --list
 2D Shapes:
	ball	Ball for the Euclidean metric.
		Required parameter(s): --radius [-R]   
	square	square (no signature).
		Required parameter(s): --width [-w]   
	lpball	Ball for the l_power metric (no signature).
		Required parameter(s): --radius [-R], --power [-p]  
	flower	Flower with k petals with radius ranging from R+/-v.
		Required parameter(s): --radius [-R], --varsmallradius [-v], --k [-k], --phi
	ngon	Regular k-gon.
		Required parameter(s): --radius [-R], --k [-k], --phi 
	accflower	Accelerated Flower with k petals.
		Required parameter(s): --radius [-R], --varsmallradius [-v], --k [-k], --phi
	ellipse	Ellipse.
		Required parameter(s): --axis1 [-A], --axis2 [-a], --phi 
 @endcode

 
 @b Example:
 @code
 # generate a flower contour with 5 petals of maximal radius 5 and small radius 3 with grid size = 0.1:
 ./shapeGenerator/contourGenerator -s flower -R 5 -v 3  -k 5 --outputGeometry result   -g 0.1
 # Display the results using gnuplot:
 $gnuplot
 [gnuplot>  plot [][-5:2]'res.geom'  using 1:6 w l  t "curvature"
 [gnuplot>  plot 'res.geom'  using 2:3 w l  t "contour flower"
@endcode

 
 You should obtain such a visualization:
 @image html resContourGenerator.png "resulting visualisation of generated contour (a) and curvature (b)."
 
 @see 
 @ref contourGenerator.cpp
 @ref shapeGenerator 
 */


/**
 * Global vectors to describe the available shapes and their
 * parameters.
 */
std::vector<std::string> shapes2D;
std::vector<std::string> shapesDesc;
std::vector<std::string> shapesParam1;
std::vector<std::string> shapesParam2;
std::vector<std::string> shapesParam3;
std::vector<std::string> shapesParam4;


/**
 * Create the static list of shapes.
 *
 */
void createList()
{
  shapes2D.push_back("ball");
  shapesDesc.push_back("Ball for the Euclidean metric.");
  shapesParam1.push_back("--radius [-R]");
  shapesParam2.push_back("");
  shapesParam3.push_back("");
  shapesParam4.push_back("");
  
  shapes2D.push_back("square");
  shapesDesc.push_back("square (no signature).");
  shapesParam1.push_back("--width [-w]");
  shapesParam2.push_back("");
  shapesParam3.push_back("");
  shapesParam4.push_back("");
  
  shapes2D.push_back("lpball");
  shapesDesc.push_back("Ball for the l_power metric (no signature).");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--power [-p]");
  shapesParam3.push_back("");
  shapesParam4.push_back("");
  
  shapes2D.push_back("flower");
  shapesDesc.push_back("Flower with k petals with radius ranging from R+/-v.");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--varsmallradius [-v],");
  shapesParam3.push_back("--k [-k],");
  shapesParam4.push_back("--phi");
  
  shapes2D.push_back("ngon");
  shapesDesc.push_back("Regular k-gon.");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--k [-k],");
  shapesParam3.push_back("--phi");
  shapesParam4.push_back("");
  
  shapes2D.push_back("accflower");
  shapesDesc.push_back("Accelerated Flower with k petals.");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--varsmallradius [-v],");
  shapesParam3.push_back("--k [-k],");
  shapesParam4.push_back("--phi");
  
  shapes2D.push_back("ellipse");
  shapesDesc.push_back("Ellipse.");
  shapesParam1.push_back("--axis1 [-A],");
  shapesParam2.push_back("--axis2 [-a],");
  shapesParam3.push_back("--phi");
  shapesParam4.push_back("");
  
  
}

/**
 * Display the shape list with parameters.
 *
 */
void displayList()
{
  trace.emphase()<<"2D Shapes:"<<std::endl;
  for(unsigned int i=0; i<shapes2D.size(); ++i)
    trace.info()<<"\t"<<shapes2D[i]<<"\t"
    <<shapesDesc[i]<<std::endl
    <<"\t\tRequired parameter(s): "
    << shapesParam1[i]<<" "
    << shapesParam2[i]<<" "
    << shapesParam3[i]<<" "
    << shapesParam4[i]<<std::endl;
  
}


/**
 * Check if a given shape is available. If not, we exist with an error.
 * If it is, we return the corresponding index in the global vectors.
 *
 * @param shapeName name of the shape to search.
 *
 * @return index of the shape in the shape vectors.
 */
unsigned int checkAndReturnIndex(const std::string &shapeName)
{
  unsigned int pos=0;
  
  while ((pos < shapes2D.size()) && (shapes2D[pos] != shapeName))
    pos++;
  
  if (pos == shapes2D.size())
  {
    trace.error() << "The specified shape has not found.";
    trace.info()<<std::endl;
    exit(1);
  }
  
  return pos;
}


/**
 * Provides ground truth for tangent and curvature estimations.
 *
 * @param s implicit shape
 * @param h grid step
 * @param r range of points at which the estimation is computed
 * @param points returned vector of points
 * @param tangents returned vector of (normalized) tangent vectors
 * @param curvature returned vector of curvature values
 *
 */
template <typename Shape, typename Range, typename Point, typename Quantity>
void
estimateGeometry(Shape& s,
                 const double& h,
                 const Range& r,
                 std::vector<Point>& points,
                 std::vector<Point>& tangents,
                 std::vector<Quantity>& curvatures) {
  
  typedef typename Range::ConstIterator ConstIterator;
  for (ConstIterator i = r.begin(); i != r.end(); ++i) {
    Point p( *i );
    p *= h;
    points.push_back(p);
  }
  
  typedef typename Range::ConstCirculator ConstCirculator;
  
  typedef ParametricShapeTangentFunctor< Shape > TangentFunctor;
  TrueLocalEstimatorOnPoints< ConstCirculator, Shape, TangentFunctor >
  trueTangentEstimator;
  trueTangentEstimator.attach(&s);
  trueTangentEstimator.init( h, r.c(), r.c());
  trueTangentEstimator.eval(r.c(), r.c(), std::back_inserter(tangents) );
  
  typedef ParametricShapeCurvatureFunctor< Shape > CurvatureFunctor;
  TrueLocalEstimatorOnPoints< ConstCirculator, Shape, CurvatureFunctor >
  trueCurvatureEstimator;
  trueCurvatureEstimator.attach(&s);
  trueCurvatureEstimator.init( h, r.c(), r.c());
  trueCurvatureEstimator.eval(r.c(), r.c(), std::back_inserter(curvatures) );
  
}

template <typename Space, typename Shape>
bool
generateContour(
                Shape & aShape,
                double h,
                const std::string & outputFormat,
                bool withGeom,
                const std::string & outputFileName  )
{
  // Types
  typedef typename Space::Point Point;
  typedef typename Space::Vector Vector;
  typedef typename Space::RealPoint RealPoint;
  typedef typename Space::Integer Integer;
  typedef HyperRectDomain<Space> Domain;
  typedef KhalimskySpaceND<Space::dimension,Integer> KSpace;
  typedef typename KSpace::SCell SCell;
  typedef typename GridCurve<KSpace>::PointsRange Range;
  typedef typename Range::ConstIterator ConstIteratorOnPoints;
  typedef typename GridCurve<KSpace>::MidPointsRange MidPointsRange;
  
  // Digitizer
  GaussDigitizer<Space,Shape> dig;
  dig.attach( aShape ); // attaches the shape.
  Vector vlow(-1,-1); Vector vup(1,1);
  dig.init( aShape.getLowerBound()+vlow, aShape.getUpperBound()+vup, h );
  Domain domain = dig.getDomain();
  // Create cellular space
  KSpace K;
  bool ok = K.init( dig.getLowerBound(), dig.getUpperBound(), true );
  if ( ! ok )
  {
    std::cerr << "[generateContour]"
    << " error in creating KSpace." << std::endl;
    return false;
  }
  try {
    // Extracts shape boundary
    SurfelAdjacency<KSpace::dimension> SAdj( true );
    SCell bel = Surfaces<KSpace>::findABel( K, dig, 10000 );
    // Getting the consecutive surfels of the 2D boundary
    std::vector<Point> points;
    Surfaces<KSpace>::track2DBoundaryPoints( points, K, SAdj, dig, bel );
    // Create GridCurve
    GridCurve<KSpace> gridcurve;
    gridcurve.initFromVector( points );
    // gridcurve contains the digital boundary to analyze.
    Range r = gridcurve.getPointsRange(); //building range
    
    if ( outputFormat == "pts" )
    {
      
      for ( ConstIteratorOnPoints it = r.begin(), it_end = r.end();
           it != it_end; ++it )
      {
        Point p = *it;
        std::cout << p[ 0 ] << " " << p[ 1 ] << std::endl;
      }
    }
    else if ( outputFormat == "fc" )
    {
      ConstIteratorOnPoints it = r.begin();
      Point p = *it++;
      std::cout << p[ 0 ] << " " << p[ 1 ] << " ";
      for ( ConstIteratorOnPoints it_end = r.end(); it != it_end; ++it )
      {
        Point p2 = *it;
        Vector v = p2 - p;
        if ( v[0 ]== 1 ) std::cout << '0';
        if ( v[ 1 ] == 1 ) std::cout << '1';
        if ( v[ 0 ] == -1 ) std::cout << '2';
        if ( v[ 1 ] == -1 ) std::cout << '3';
        p = p2;
      }
      // close freemanchain if necessary.
      Point p2= *(r.begin());
      Vector v = p2 - p;
      if ( v.norm1() == 1 )
      {
        if ( v[ 0 ] == 1 ) std::cout << '0';
        if ( v[ 1 ] == 1 ) std::cout << '1';
        if ( v[ 0 ] == -1 ) std::cout << '2';
        if ( v[ 1 ] == -1 ) std::cout << '3';
      }
      std::cout << std::endl;
    }
    
    if (withGeom)
    {
      // write geometry of the shape
      std::stringstream s;
      s << outputFileName << ".geom";
      std::ofstream outstream(s.str().c_str()); //output stream
      if (!outstream.is_open()) return false;
      else {
        outstream << "# " << outputFileName << std::endl;
        outstream << "# Pointel (x,y), Midpoint of the following linel (x',y')" << std::endl;
        outstream << "# id x y tangentx tangenty curvaturexy"
        << " x' y' tangentx' tangenty' curvaturex'y'" << std::endl;
        
        std::vector<RealPoint> truePoints, truePoints2;
        std::vector<RealPoint> trueTangents, trueTangents2;
        std::vector<double> trueCurvatures, trueCurvatures2;
        
        estimateGeometry<Shape, Range, RealPoint, double>
        (aShape, h, r, truePoints, trueTangents, trueCurvatures);
        
        estimateGeometry<Shape, MidPointsRange, RealPoint, double>
        (aShape, h, gridcurve.getMidPointsRange(), truePoints2, trueTangents2, trueCurvatures2);
        
        
        unsigned int n = (unsigned int)r.size();
        for (unsigned int i = 0; i < n; ++i ) {
          outstream << std::setprecision( 15 ) << i
          << " " << truePoints[ i ][ 0 ]
          << " " << truePoints[ i ][ 1 ]
          << " " << trueTangents[ i ][ 0 ]
          << " " << trueTangents[ i ][ 1 ]
          << " " << trueCurvatures[ i ]
          << " " << truePoints2[ i ][ 0 ]
          << " " << truePoints2[ i ][ 1 ]
          << " " << trueTangents2[ i ][ 0 ]
          << " " << trueTangents2[ i ][ 1 ]
          << " " << trueCurvatures2[ i ]
          << std::endl;
        }
        
      }
      outstream.close();
    }
    
    /////////////////
    
  }
  catch ( InputException e )
  {
    std::cerr << "[generateContour]"
    << " error in finding a bel." << std::endl;
    return false;
  }
  return true;
}

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam(std::string param)
{
  trace.error() <<" Parameter: "<<param<<" is required..";
  trace.info()<<std::endl;
  exit(1);
}

///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
  // parse command line CLI ----------------------------------------------
  CLI::App app;
  std::string shapeName;
  std::string outputFileName;
  std::string outputFormat {"pts"};
  double radius;
  double power {2.0};
  double smallradius {5};
  double varsmallradius {5};
  double cx {0.0}, cy {0.0};
  double h {1.0};
  unsigned int k {3};
  double phi {0.0};
  double width {10.0};
  double axis1, axis2;

  app.description("Generates multigrid contours of 2d digital shapes using DGtal library.\n Typical use example:\n \t contourGenerator --shape <shapeName> [requiredParam] [otherOptions]\n");
  auto listOpt = app.add_flag("--list,-l","List all available shapes");
  auto shapeNameOpt = app.add_option("--shape,-s", shapeName, "Shape name");
  auto radiusOpt = app.add_option("--radius,-R", radius, "Radius of the shape" );
  auto axis1Opt = app.add_option("--axis1,-A", axis1, "Half big axis of the shape (ellipse)" );
  auto axis2Opt = app.add_option("--axis2,-a", axis2, "Half small axis of the shape (ellipse)" );
  auto smallradiusOpt = app.add_option("--smallradius,-r", smallradius, "Small radius of the shape (default 5)", true);
  auto varsmallradiusOpt = app.add_option("--varsmallradius,-v", varsmallradius, "Variable small radius of the shape (default 5)", true );
  auto kOpt = app.add_option("-k", k, "Number of branches or corners the shape (default 3)", true );
  auto phiOpt = app.add_option("--phi", phi, "Phase of the shape (in radian, default 0.0)", true );
  auto widthOpt = app.add_option("--width,-w", width, "Width of the shape (default 10.0)", true );
  auto powerOpt = app.add_option("--power,-p", power, "Power of the metric (default 2.0)", true );
  app.add_option("--center_x,-x", cx, "x-coordinate of the shape center (default 0.0)", true );
  app.add_option("--center_y,-y", cy, "y-coordinate of the shape center (default 0.0)", true );
  app.add_option("--gridstep,-g", h, "Gridstep for the digitization (default 1.0)", true );
  auto outputFormatOpt = app.add_option("--format,-f", outputFormat, "Output format:\n\t  List of pointel coordinates {pts}\n\t  Freeman chaincode Vector {fc} (default pts)", true );
  auto outputFileNameOpt = app.add_option("--outputGeometry,-o", outputFileName, "Base name of the file containing the shape geometry (points, tangents, curvature)" );
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  //List creation
  createList();

  if ( listOpt->count() > 0 )
  {
    displayList();
    return 0;
  }
  
  if(shapeNameOpt->count()==0) missingParam("--shape");
  bool withGeom = true;
  if (outputFileNameOpt->count()==0) withGeom = false;
  
  //We check that the shape is known
  unsigned int id = checkAndReturnIndex(shapeName);
  
  // standard types
  typedef Z2i::Space Space;
  typedef Space::RealPoint RealPoint;
  
  RealPoint center( cx, cy );

  if (id ==0)
  {
    if (radiusOpt->count()==0) missingParam("--radius");
    Ball2D<Space> ball(Z2i::Point(0,0), radius);
    generateContour<Space>( ball, h, outputFormat, withGeom, outputFileName );
  }
  else if (id ==1)
  {
    //if (widthOpt->count()==0) missingParam("--width");
    ImplicitHyperCube<Space> object(Z2i::Point(0,0), width/2);
    trace.error()<< "Not available.";
    trace.info()<<std::endl;
  }
  else if (id ==2)
  {
    //if (powerOpt->count()==0) missingParam("--power");
    if (radiusOpt->count()==0) missingParam("--radius");
    ImplicitRoundedHyperCube<Space> ball(Z2i::Point(0,0), radius, power);
    trace.error()<< "Not available.";
    trace.info()<<std::endl;
  }
  else if (id ==3)
  {
    //if (varsmallradiusOpt->count()==0) missingParam("--varsmallradius");
    if (radiusOpt->count()==0) missingParam("--radius");
    //if (kOpt->count()==0) missingParam("--k");
    //if (phiOpt->count()==0) missingParam("--phi");
    Flower2D<Space> flower( center, radius, varsmallradius, k, phi );
    generateContour<Space>( flower, h, outputFormat, withGeom, outputFileName  );
  }
  else if (id ==4)
  {
    if (radiusOpt->count()==0) missingParam("--radius");
    //if (kOpt->count()==0) missingParam("--k");
    //if (phiOpt->count()==0) missingParam("--phi");
    NGon2D<Space> object( center, radius, k, phi );
    generateContour<Space>( object, h, outputFormat, withGeom, outputFileName  );
  }
  else if (id ==5)
  {
    //if (varsmallradiusOpt->count()==0) missingParam("--varsmallradius");
    if (radiusOpt->count()==0) missingParam("--radius");
    //if (kOpt->count()==0) missingParam("--k");
    //if (phiOpt->count()==0) missingParam("--phi");
    AccFlower2D<Space> accflower( center, radius, varsmallradius, k, phi );
    generateContour<Space>( accflower, h, outputFormat, withGeom, outputFileName  );
  }
  else if (id ==6)
  {
    if (axis1Opt->count()==0) missingParam("--axis1");
    if (axis2Opt->count()==0) missingParam("--axis2");
    //if (phiOpt->count()==0) missingParam("--phi");
    Ellipse2D<Space> ellipse( center, axis1, axis2, phi );
    generateContour<Space>( ellipse, h, outputFormat, withGeom, outputFileName  ); 
  } 
}
