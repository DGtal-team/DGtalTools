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
 * @file LengthEstimator.cpp
 * @ingroup Tools
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr ) 
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr ) 
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS,
 * France
 *
 *
 * @date 2011/07/07
 *
 * DGtal tool for length estimations on implicit shapes
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
#include "DGtal/base/Clock.h"

//space / domain
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/topology/KhalimskySpaceND.h"

//shape and digitizer
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/geometry/curves/GridCurve.h"

//estimators
#include "DGtal/geometry/curves/estimation/TrueLocalEstimatorOnPoints.h"
#include "DGtal/geometry/curves/estimation/TrueGlobalEstimatorOnPoints.h"
#include "DGtal/geometry/curves/estimation/ParametricShapeArcLengthFunctor.h"

#include "DGtal/geometry/curves/estimation/L1LengthEstimator.h"
#include "DGtal/geometry/curves/estimation/BLUELocalLengthEstimator.h"
#include "DGtal/geometry/curves/estimation/RosenProffittLocalLengthEstimator.h"
#include "DGtal/geometry/curves/estimation/MLPLengthEstimator.h"
#include "DGtal/geometry/curves/estimation/FPLengthEstimator.h"
#include "DGtal/geometry/curves/estimation/DSSLengthEstimator.h"

using namespace DGtal;



/**
 @page lengthEstimators lengthEstimators
 
 @brief Generates multigrid length estimations of paramteric shapes using DGtal library.

It will output length estimations (and timings) using several algorithms for decreasing grid steps.

@b Usage: 	LengthEstimators [options] --shape <shapeName>



 @b Allowed @b options @b are : 
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
  --hMin FLOAT=0.0001                   Minimum value for the grid step h (double, default 0.0001)
  --steps INT=32                        Number of multigrid steps between 1 and hMin (integer, default 32)

 @endcode

 The list of 2D shape are :
  -ball	Ball for the Euclidean metric.
       	Required parameter(s): --radius [-R]
  - square	square (no signature).
    Required parameter(s): --width [-w]
  - lpball	Ball for the l_power metric (no signature).
    Required parameter(s): --radius [-R], --power [-p]
  - flower	Flower with k petals.
    Required parameter(s): --radius [-R], --varsmallradius [-v], --k [-k], --phi
  - ngon	Regular k-gon.
    Required parameter(s): --radius [-R], --k [-k], --phi 
  - accflower	Accelerated Flower with k petals.
    Required parameter(s): --radius [-R], --varsmallradius [-v], --k [-k], --phi
  - ellipse	Ellipse.
    Required parameter(s): --axis1 [-A], --axis2 [-a], --phi 

 @b Example: 

Application of the multigrid length estimation on a flower shape with 5 petals and maximal radius 20 and min 5: 
@code
$ lengthEstimators -s flower -k 5 -R 20 -r 5 --steps 256 > length.dat 
@endcode

You can display using gnuplot:

@code
$ gnuplot 
gnuplot> plot [0:][-0.5: 90] 'length.dat' using 1:(($7-$3)*($7-$3)) w l title "squared error length estimation using DSS", 'length.dat' using 1:(($8-$3)*($8-$3)) w l title "squared error length estimation using MLP", 'length.dat' using 1:(($9-$3)*($9-$3)) w l title "squared error length estimation using FP" linewidth 2
@endcode



 You should obtain such a result:
 @image html resLengthEstimators.png "Resulting visualization."
 
 @see
 @ref lengthEstimators.cpp

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
  shapesDesc.push_back("Flower with k petals.");
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
    << shapesDesc[i]<<std::endl
    << "\t\tRequired parameter(s): "
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
unsigned int checkAndRetrunIndex(const std::string &shapeName)
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

template <typename Shape, typename Space>
bool
lengthEstimators( const std::string & /*name*/,
      Shape & aShape, 
      double h )
{
  // Types
  typedef typename Space::Point Point;
  typedef typename Space::Vector Vector;
  typedef typename Space::Integer Integer;
  typedef HyperRectDomain<Space> Domain;
  typedef KhalimskySpaceND<Space::dimension,Integer> KSpace;
  typedef typename KSpace::SCell SCell;
  typedef typename GridCurve<KSpace>::PointsRange PointsRange;
  typedef typename GridCurve<KSpace>::ArrowsRange ArrowsRange;

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
    std::cerr << "[lengthEstimators]"
    << " error in creating KSpace." << std::endl;
    return false;
  }
  try
  {
    // Extracts shape boundary
    SurfelAdjacency<KSpace::dimension> SAdj( true );
    SCell bel = Surfaces<KSpace>::findABel( K, dig, 10000 );
    // Getting the consecutive surfels of the 2D boundary
    std::vector<Point> points;
    Surfaces<KSpace>::track2DBoundaryPoints( points, K, SAdj, dig, bel );
    // Create GridCurve
    GridCurve<KSpace> gridcurve;
    gridcurve.initFromPointsVector( points );
    // Ranges
    ArrowsRange ra = gridcurve.getArrowsRange(); 
    PointsRange rp = gridcurve.getPointsRange(); 

    // Estimations
    typedef typename PointsRange::ConstIterator ConstIteratorOnPoints; 
    typedef ParametricShapeArcLengthFunctor< Shape > Length;
    TrueGlobalEstimatorOnPoints< ConstIteratorOnPoints, Shape, Length  >  trueLengthEstimator;
    trueLengthEstimator.attach( aShape );

    L1LengthEstimator< typename PointsRange::ConstCirculator > l1length;
    DSSLengthEstimator< typename PointsRange::ConstCirculator > DSSlength;
    MLPLengthEstimator< typename PointsRange::ConstCirculator > MLPlength;
    FPLengthEstimator< typename PointsRange::ConstCirculator > FPlength;
    BLUELocalLengthEstimator< typename ArrowsRange::ConstCirculator > BLUElength;
    RosenProffittLocalLengthEstimator< typename ArrowsRange::ConstCirculator > RosenProffittlength;

    // Output
    double trueValue = trueLengthEstimator.eval();
    double l1, blue, rosen,dss,mlp,fp;
    double Tl1, Tblue, Trosen,Tdss,Tmlp,Tfp;
    
    Clock c;

    //Length evaluation & timing
    c.startClock();
    l1 = l1length.eval( rp.c(), rp.c(), h );
    Tl1 = c.stopClock();

    c.startClock();
    blue = BLUElength.eval( ra.c(), ra.c(), h );
    Tblue = c.stopClock();

    c.startClock();
    rosen = RosenProffittlength.eval( ra.c(), ra.c(), h );
    Trosen = c.stopClock();

    c.startClock();
    dss = DSSlength.eval( rp.c(), rp.c(), h );
    Tdss = c.stopClock();

    c.startClock();
    mlp = MLPlength.eval( rp.c(), rp.c(), h );
    Tmlp = c.stopClock();

    c.startClock();
    fp = FPlength.eval( rp.c(), rp.c(), h );
    Tfp = c.stopClock();

    std::cout << std::setprecision( 15 ) << h << " " << rp.size() << " " << trueValue 
    << " " << l1
    << " " << blue
    << " " << rosen
    << " " << dss
    << " " << mlp
    << " " << fp
    << " " << Tl1
    << " " << Tblue
    << " " << Trosen
    << " " << Tdss
    << " " << Tmlp
    << " " << Tfp
    << std::endl;
    return true;
  }
  catch ( InputException e )
  {
    std::cerr << "[lengthEstimators]"
     << " error in finding a bel." << std::endl;
    return false;
  }
}
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
  // parse command line CLI ----------------------------------------------
  CLI::App app;
  std::string shapeName;
  double hMin {0.0001};
  int nbSteps {32};
  double radius;
  double power {2.0};
  double smallradius {5};
  double varsmallradius {5};
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
  app.add_option("--hMin", hMin, "Minimum value for the grid step h (double, default 0.0001)", true );
  app.add_option("--steps", nbSteps, "Number of multigrid steps between 1 and hMin (integer, default 32)", true );
  
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

  //Parse options
  if(shapeNameOpt->count()==0) missingParam("--shape");
    
  //We check that the shape is known
  unsigned int id = checkAndRetrunIndex(shapeName);


///////////////////////////////////
  std::cout << "#h nbPoints true-length L1 BLUE RosenProffit "
       << "DSS MLP FP Time-L1 Time-BLUE Time-RosenProffitt "
       << "Time-DSS Time-MLP Time-FP" << std::endl;
  std::cout << "# timings are given in msec." << std::endl;

  double h = 1; 
  double step = exp( log(hMin) / (double)nbSteps);
  while (h > hMin) {

    if (id ==0) ///ball
      {
        if (radiusOpt->count()==0) missingParam("--radius");
        
        Ball2D<Z2i::Space> ball(Z2i::Point(0,0), radius);
        
        lengthEstimators<Ball2D<Z2i::Space>,Z2i::Space>("ball",ball,h); 
      }
    else
      if (id ==1)
        {
          //if (widthOpt->count()==0) missingParam("--width");
          
          ImplicitHyperCube<Z2i::Space> object(Z2i::Point(0,0), width/2);
        
                trace.error()<< "Not available.";
                trace.info()<<std::endl;
        }
      else
        if (id ==2)
          {
            //if (powerOpt->count()==0) missingParam("--power");
            if (radiusOpt->count()==0) missingParam("--radius");
            
            ImplicitRoundedHyperCube<Z2i::Space> ball(Z2i::Point(0,0), radius, power);

                trace.error()<< "Not available.";
                trace.info()<<std::endl;
          }
        else
          if (id ==3)
            {
              //if (varsmallradiusOpt->count()==0) missingParam("--varsmallradius");
              if (radiusOpt->count()==0) missingParam("--radius");
              //if (kOpt->count()==0) missingParam("--k");
              //if (phiOpt->count()==0) missingParam("--phi");
              
              Flower2D<Z2i::Space> flower(Z2i::Point(0,0), radius, varsmallradius,k,phi);

              lengthEstimators<Flower2D<Z2i::Space>,Z2i::Space>("flower",flower,h); 
            }
    else
      if (id ==4)
        {
          if (radiusOpt->count()==0) missingParam("--radius");
          //if (kOpt->count()==0) missingParam("--k");
          //if (phiOpt->count()==0) missingParam("--phi");
          
          NGon2D<Z2i::Space> object(Z2i::Point(0,0), radius,k,phi);

          lengthEstimators<NGon2D<Z2i::Space>,Z2i::Space>("NGon",object,h); 

        }
      else
        if (id ==5)
          {
            //if (varsmallradiusOpt->count()==0) missingParam("--varsmallradius");
            if (radiusOpt->count()==0) missingParam("--radius");
            //if (kOpt->count()==0) missingParam("--k");
            //if (phiOpt->count()==0) missingParam("--phi");
                
            AccFlower2D<Z2i::Space> flower(Z2i::Point(0,0), radius, varsmallradius,k,phi);
            lengthEstimators<AccFlower2D<Z2i::Space>,Z2i::Space>("accFlower",flower,h); 

          } 
        else
          //if (id ==6)
          {
            if (axis1Opt->count()==0) missingParam("--axis1");
            if (axis2Opt->count()==0) missingParam("--axis2");
            //if (phiOpt->count()==0) missingParam("--phi");
            
            Ellipse2D<Z2i::Space> ell(Z2i::Point(0,0), axis1, axis2,phi);

      lengthEstimators<Ellipse2D<Z2i::Space>,Z2i::Space>("Ellipse",ell,h); 
    } 

    h = h * step;
  }
  return 0;
}
