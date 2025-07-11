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
 * @file 3dCurveViewer.cpp
 * @ingroup visualisationTools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 *
 * @date 2013/05/09
 *
 * This file is part of the DGtal library
 */


/**
 @page Doc3dCurveViewer 3dCurveViewer
 
 @brief  Displays a 3D curve given as the input filename (with possibly projections and/or tangent information) by using PolyscopeViewer.
 

 @b Usage:  3dCurveViewer [options] input
 

 @b Allowed @b options @b are :
 
 @code

 Positionals:
   1 TEXT:FILE REQUIRED                  the name of the text file containing the list of 3D points (x y z per line).

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         the name of the text file containing the list of 3D points (x y z per line).
   -b,--box INT=0                        specifies the the tightness of the bounding box around the curve with a given integer displacement <arg> to enlarge it (0 is tight)
   -v,--viewBox TEXT:{WIRED,COLORED}=WIRED
                                         displays the bounding box, <arg>=WIRED means that only edges are displayed, <arg>=COLORED adds colors for planes (XY is red, XZ green, YZ, blue).
   -C,--curve3d                          displays the 3D curve.
   -c,--curve2d                          displays the 2D projections of the 3D curve on the bounding box.
   -3,--cover3d                          displays the 3D tangential cover of the curve.
   -2,--cover2d                          displays the 2D projections of the 3D tangential cover of the curve
   -n,--nbColors INT=3                   sets the number of successive colors used for displaying 2d and 3d maximal segments (default is 3: red, green, blue)
   -t,--tangent                          displays the tangents to the curve
   
 @endcode


 @b Example: 
 
 @code
 $  3dCurveViewer -C -b 1 -3 -2 -c ${DGtal}/examples/samples/sinus.dat
 @endcode

 You should obtain such a visualisation:
 @image html res3dCurveViewer.png "resulting visualisation of 3d curve with tangential cover."
 

 @see
 @ref 3dCurveViewer.cpp

 */


#include <iostream>
#include <iterator>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>

#include "CLI11.hpp"

#include "DGtal/base/Common.h"
#include "DGtal/base/Exceptions.h"
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/geometry/curves/GridCurve.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"

#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include "DGtal/io/viewers/PolyscopeViewer.h"


using namespace DGtal;
using namespace std;

const Color  AXIS_COLOR( 0, 0, 0, 255 );
const double AXIS_LINESIZE = 0.1;
const Color  XY_COLOR( 0, 0, 255, 50 );
const Color  XZ_COLOR( 0, 255, 0, 50 );
const Color  YZ_COLOR( 255, 0, 0, 50 );
const Color  CURVE3D_COLOR( 100, 100, 140, 128 );
const Color  CURVE2D_COLOR( 200, 200, 200, 100 );
const double MS3D_LINESIZE = 0.05;

///////////////////////////////////////////////////////////////////////////////
// Functions for displaying the tangential cover of a 3D curve.
///////////////////////////////////////////////////////////////////////////////
template <typename Point, typename RealPoint, typename space, typename kspace >
void displayAxes( PolyscopeViewer<space, kspace> & viewer,
                  const Point & lowerBound, const Point & upperBound,
      const std::string & mode )
{
  RealPoint p0( (double)lowerBound[ 0 ]-0.5,
                (double)lowerBound[ 1 ]-0.5,
                (double)lowerBound[ 2 ]-0.5 );
  RealPoint p1( (double)upperBound[ 0 ]-0.5,
                (double)upperBound[ 1 ]-0.5,
                (double)upperBound[ 2 ]-0.5 );
  if ( ( mode == "WIRED" ) || ( mode == "COLORED" ) )
    {
      viewer.drawColor(AXIS_COLOR);
      viewer.drawLine( DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p0[ 2 ]),
          DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p0[ 2 ]));
      viewer.drawLine( DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p0[ 2 ]),
          DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p0[ 2 ]));
      viewer.drawLine( DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p0[ 2 ]),
          DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p1[ 2 ]));
      viewer.drawLine( DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p0[ 2 ]),
          DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p0[ 2 ]));
      viewer.drawLine( DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p0[ 2 ]),
          DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p1[ 2 ]));
      viewer.drawLine( DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p0[ 2 ]),
          DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p0[ 2 ]));
      viewer.drawLine( DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p0[ 2 ]),
          DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p1[ 2 ]));
      viewer.drawLine( DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p1[ 2 ]),
          DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p1[ 2 ]));
      viewer.drawLine( DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p1[ 2 ]),
          DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p1[ 2 ]));
      viewer.drawLine( DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p0[ 2 ]),
          DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p1[ 2 ]));
      viewer.drawLine( DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p1[ 2 ]),
          DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p1[ 2 ]));
      viewer.drawLine( DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p1[ 2 ]),
          DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p1[ 2 ]));
    }
  if ( mode == "COLORED" )
    {
      viewer.drawColor(XY_COLOR);
      viewer.drawQuad(DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p1[ 2 ]),
         DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p1[ 2 ]),
         DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p1[ 2 ]),
         DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p1[ 2 ]) );
      viewer.drawColor(XZ_COLOR);
      viewer.drawQuad(DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p1[ 2 ]),
         DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p1[ 2 ]),
         DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p0[ 2 ]),
         DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p0[ 2 ]));
      viewer.drawColor(YZ_COLOR);
      viewer.drawQuad(DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p1[ 2 ]),
         DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p1[ 2 ]),
         DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p0[ 2 ]),
         DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p0[ 2 ]));
    }
}

template <typename KSpace, typename StandardDSS6Computer, typename space, typename kspace >
void displayDSS3d( PolyscopeViewer<space, kspace> & viewer,
       const KSpace & ks, const StandardDSS6Computer & dss3d,
       const DGtal::Color & color3d )
{
  viewer << color3d << dss3d;
}

template <typename Point1, typename Point2>
void assign( Point1 & p1, const Point2 & p2 )
{
  p1[ 0 ] = p2[ 0 ];
  p1[ 1 ] = p2[ 1 ];
  p1[ 2 ] = p2[ 2 ];
}

template <typename KSpace, typename StandardDSS6Computer, typename space, typename kspace >
void displayDSS3dTangent( PolyscopeViewer<space, kspace> & viewer,
        const KSpace & ks, const StandardDSS6Computer & dss3d,
        const DGtal::Color & color3d )
{
  typedef typename StandardDSS6Computer::Point3d Point;
  typedef typename StandardDSS6Computer::PointR3d PointR3d;
  typedef DGtal::PointVector<3,double> PointD3d;
  Point directionZ3;
  PointR3d interceptR, thicknessR;
  PointD3d P1, P2, direction;
  dss3d.getParameters( directionZ3, interceptR, thicknessR );
  
  PointD3d intercept;
  intercept[0] = (double) NumberTraits<int>::castToInt64_t ( interceptR[0].first ) / (double) NumberTraits<int>::castToInt64_t ( interceptR[0].second );
  intercept[1] = (double) NumberTraits<int>::castToInt64_t ( interceptR[1].first ) / (double) NumberTraits<int>::castToInt64_t ( interceptR[1].second );
  intercept[2] = (double) NumberTraits<int>::castToInt64_t ( interceptR[2].first ) / (double) NumberTraits<int>::castToInt64_t ( interceptR[2].second );
  
  PointD3d thickness;
  thickness[0] = (double) NumberTraits<int>::castToInt64_t ( thicknessR[0].first ) / (double) NumberTraits<int>::castToInt64_t ( thicknessR[0].second );
  thickness[1] = (double) NumberTraits<int>::castToInt64_t ( thicknessR[1].first ) / (double) NumberTraits<int>::castToInt64_t ( thicknessR[1].second );
  thickness[2] = (double) NumberTraits<int>::castToInt64_t ( thicknessR[2].first ) / (double) NumberTraits<int>::castToInt64_t ( thicknessR[2].second );
  
  assign( direction, directionZ3 );
  direction /= direction.norm();
  assign( P1, *dss3d.begin() );
  assign( P2, *(dss3d.end()-1) );
  double t1 = (P1 - intercept).dot( direction );
  double t2 = (P2 - intercept).dot( direction );

  PointD3d Q1 = intercept + t1 * direction;
  PointD3d Q2 = intercept + t2 * direction;
  viewer.drawColor(color3d);
  viewer.drawLine( DGtal::Z3i::RealPoint(Q1[ 0 ]-0.5, Q1[ 1 ]-0.5, Q1[ 2 ]-0.5),
      DGtal::Z3i::RealPoint(Q2[ 0 ]-0.5, Q2[ 1 ]-0.5, Q2[ 2 ]-0.5));
}

template <typename KSpace, typename StandardDSS6Computer, typename space, typename kspace >
void displayProj2d( PolyscopeViewer<space, kspace> & viewer,
        const KSpace & ks, const StandardDSS6Computer & dss3d,
        const DGtal::Color & color2d )
{
  typedef typename StandardDSS6Computer::ArithmeticalDSSComputer2d ArithmeticalDSSComputer2d;
  typedef typename ArithmeticalDSSComputer2d::ConstIterator ConstIterator2d;
  typedef typename ArithmeticalDSSComputer2d::Point Point2d;
  typedef typename KSpace::Cell Cell;
  typedef typename KSpace::Point Point3d;
  Point3d b = ks.lowerBound();
  for ( DGtal::Dimension i = 0; i < 3; ++i )
    {
      const ArithmeticalDSSComputer2d & dss2d = dss3d.arithmeticalDSS2d( i );
      for ( ConstIterator2d itP = dss2d.begin(), itPEnd = dss2d.end(); itP != itPEnd; ++itP )
  {
    Point2d p = *itP;
    Point3d q;
    switch (i) {
    case 0: q = Point3d( 2*b[ i ]  , 2*p[ 0 ]+1, 2*p[ 1 ]+1 ); break;
    case 1: q = Point3d( 2*p[ 0 ]+1, 2*b[ i ]  , 2*p[ 1 ]+1 ); break;
    case 2: q = Point3d( 2*p[ 0 ]+1, 2*p[ 1 ]+1, 2*b[ i ]   ); break;
    }
    Cell c = ks.uCell( q );
    viewer << color2d << c;
  }
    }
}

template <typename KSpace, typename StandardDSS6Computer, typename space, typename kspace >
void displayDSS2d( PolyscopeViewer<space, kspace> & viewer,
       const KSpace & ks, const StandardDSS6Computer & dss3d,
       const DGtal::Color & color2d )
{
  typedef typename StandardDSS6Computer::ConstIterator ConstIterator3d;
  typedef typename StandardDSS6Computer::ArithmeticalDSSComputer2d ArithmeticalDSSComputer2d;
  typedef typename ArithmeticalDSSComputer2d::ConstIterator ConstIterator2d;
  typedef typename ArithmeticalDSSComputer2d::Point Point2d;
  typedef typename KSpace::Cell Cell;
  typedef typename KSpace::Point Point3d;
  typedef DGtal::PointVector<2,double> PointD2d;

  Point3d b = ks.lowerBound();
  for ( DGtal::Dimension i = 0; i < 3; ++i )
    {
      const typename ArithmeticalDSSComputer2d::Primitive & dss2d
  = dss3d.arithmeticalDSS2d( i ).primitive();
      // draw 2D bounding boxes for each arithmetical dss 2D.
      std::vector<PointD2d> pts2d;
      pts2d.push_back( dss2d.project(dss2d.back(), dss2d.Uf()) );
      pts2d.push_back( dss2d.project(dss2d.back(), dss2d.Lf()) );
      pts2d.push_back( dss2d.project(dss2d.front(), dss2d.Lf()) );
      pts2d.push_back( dss2d.project(dss2d.front(), dss2d.Uf()) );
      std::vector<Point3d> bb;
      Point3d p3;
      for ( unsigned int j = 0; j < pts2d.size(); ++j )
  {
    switch (i) {
    case 0: p3[0] = (double) b[ i ]-0.5; p3[1] = pts2d[ j ][ 0 ];  p3[2] = pts2d[ j ][ 1 ]; break;
    case 1: p3[0] = pts2d[ j ][ 0 ];  p3[1] = (double) b[ i ]-0.5; p3[2] = pts2d[ j ][ 1 ];     break;
    case 2: p3[0] = pts2d[ j ][ 0 ];  p3[1] = pts2d[ j ][ 1 ];     p3[2] = (double) b[ i ]-0.5; break;
    }
    bb.push_back( p3 );
  }
      for ( unsigned int j = 0; j < pts2d.size(); ++j ){
  viewer.drawColor(color2d);
  viewer.drawLine( DGtal::Z3i::RealPoint(bb[ j ][0], bb[ j ][1], bb[ j ][2]),
                        DGtal::Z3i::RealPoint(bb[ (j+1)%4 ][0], bb[ (j+1)%4 ][1], bb[ (j+1)%4 ][2]));
      }
    } // for ( DGtal::Dimension i = 0; i < 3; ++i )
}

/**
 * segmentation test
 *
 */
template <typename KSpace, typename PointIterator, typename space, typename kspace >
bool displayCover( PolyscopeViewer<space, kspace> & viewer,
       const KSpace & ks, PointIterator b, PointIterator e,
       bool dss3d, bool proj2d, bool dss2d, bool tangent,
       int nbColors )
{
  typedef typename PointIterator::value_type Point;
  typedef StandardDSS6Computer<PointIterator,int,4> SegmentComputer;
  typedef SaturatedSegmentation<SegmentComputer> Decomposition;
  typedef typename Decomposition::SegmentComputerIterator SegmentComputerIterator;
  typedef typename SegmentComputer::ArithmeticalDSSComputer2d ArithmeticalDSSComputer2d;
  SegmentComputer algo;
  Decomposition theDecomposition(b, e, algo);

  HueShadeColorMap<int> cmap_hue( 0, nbColors, 1 );

  unsigned int c = 0;
  for ( SegmentComputerIterator i = theDecomposition.begin();
        i != theDecomposition.end(); ++i)
    {
      SegmentComputer ms3d(*i);
      const ArithmeticalDSSComputer2d & dssXY = ms3d.arithmeticalDSS2dXY();
      const ArithmeticalDSSComputer2d & dssXZ = ms3d.arithmeticalDSS2dXZ();
      const ArithmeticalDSSComputer2d & dssYZ = ms3d.arithmeticalDSS2dYZ();
      Point f = *ms3d.begin();
      Point l = *(ms3d.end() - 1);
      trace.info() << "- " << c
                   << " MS3D,"
                   << " [" << f[ 0 ] << "," << f[ 1 ] << ","<< f[ 2 ] << "]"
                   << "->[" << l[ 0 ] << "," << l[ 1 ] << ","<< l[ 2 ] << "]"
                   << ", XY("
                   << dssXY.a() << "," << dssXY.b() << "," << dssXY.mu()
                   << "), XZ("
                   << dssXZ.a() << "," << dssXZ.b() << "," << dssXZ.mu()
                   << "), YZ("
                   << dssYZ.a() << "," << dssYZ.b() << "," << dssYZ.mu()
                   << ")" << std::endl;
      //trace.info() << ms3d << std::endl;  // information

      Color color = cmap_hue( c );
      if ( tangent ) displayDSS3dTangent( viewer, ks, ms3d, color );
      if ( dss3d )   displayDSS3d( viewer, ks, ms3d, color );
      if ( dss2d )   displayDSS2d( viewer, ks, ms3d, color );
      if ( proj2d )  displayProj2d( viewer, ks, ms3d, CURVE2D_COLOR );
      c++;
    }
  return true;
}

/**
   Main function.

   @param argc the number of parameters given on the line command.

   @param argv an array of C-string, such that argv[0] is the name of
   the program, argv[1] the first parameter, etc.
*/
int main(int argc, char **argv)
{
  typedef SpaceND<3,int> Z3;
  typedef KhalimskySpaceND<3,int> K3;
  typedef Z3::Point Point;
  typedef Z3::RealPoint RealPoint;
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  int b {0};
  std::string viewBox {"WIRED"};
  bool curve3d {false};
  bool curve2d {false};
  bool cover3d {false};
  bool cover2d {false};
  bool tangent {false};
  int nbColors {3};
  
  app.description("Display a 3D curve given as the <input> filename (with possibly projections and/or tangent information) by using QGLviewer.\n Example:\n 3dCurveViewer -C -b 1 -3 -2 -c ${DGtal}/examples/samples/sinus.dat\n");
  app.add_option("-i,--input,1", inputFileName, "the name of the text file containing the list of 3D points (x y z per line)." )
  ->required()
  ->check(CLI::ExistingFile);
  app.add_option("--box,-b",b, "specifies the the tightness of the bounding box around the curve with a given integer displacement <arg> to enlarge it (0 is tight)", true);
  app.add_option("--viewBox,-v",viewBox, "displays the bounding box, <arg>=WIRED means that only edges are displayed, <arg>=COLORED adds colors for planes (XY is red, XZ green, YZ, blue).", true )
   -> check(CLI::IsMember({"WIRED", "COLORED"}));
  
  app.add_flag("--curve3d,-C", curve3d, "displays the 3D curve.");
  app.add_flag("--curve2d,-c", curve2d, "displays the 2D projections of the 3D curve on the bounding box.");
  app.add_flag("--cover3d,-3", curve2d, "displays the 3D tangential cover of the curve.");
  app.add_flag("--cover2d,-2", cover2d, "displays the 2D projections of the 3D tangential cover of the curve" );
  app.add_option("--nbColors,-n", nbColors, "sets the number of successive colors used for displaying 2d and 3d maximal segments (default is 3: red, green, blue)", true);

  app.add_flag("--tangent,-t", tangent, "displays the tangents to the curve" );


  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  
  // Create curve 3D.
  vector<Point> sequence;
  fstream inputStream;
  inputStream.open ( inputFileName.c_str(), ios::in);
  try {
    sequence = PointListReader<Point>::getPointsFromInputStream( inputStream );
    if ( sequence.size() == 0) throw IOException();
  }
  catch (DGtal::IOException & ioe) {
    trace.error() << "Size is null." << std::endl;
  }
  inputStream.close();

  // start viewer
  PolyscopeViewer<> viewer;
  trace.beginBlock ( "Tool 3dCurveViewer" );

  // ----------------------------------------------------------------------
  // Create domain and curve.
  Point lowerBound = sequence[ 0 ];
  Point upperBound = sequence[ 0 ];
  for ( unsigned int j = 1; j < sequence.size(); ++j )
    {
      lowerBound = lowerBound.inf( sequence[ j ] );
      upperBound = upperBound.sup( sequence[ j ] );
    }
  lowerBound -= Point::diagonal( b );
  upperBound += Point::diagonal( b+1 );
  K3 ks; ks.init( lowerBound, upperBound, true );
  GridCurve<K3> gc( ks );
  try {
    gc.initFromPointsVector( sequence );
  } catch (DGtal::ConnectivityException& /*ce*/) {
    throw ConnectivityException();
  }

  // ----------------------------------------------------------------------
  // Display axes.
  if ( viewBox != "" )
    displayAxes<Point,RealPoint, Z3i::Space, Z3i::KSpace>( viewer, lowerBound, upperBound, viewBox );
  // Display 3D tangential cover.
  bool res = displayCover( viewer, ks, sequence.begin(), sequence.end(),
                           cover3d,  curve2d, cover2d, tangent, nbColors );
  // Display 3D curve points.
  if ( curve3d )
    viewer << CURVE3D_COLOR
     << gc.getPointsRange()
     << sequence.back(); // curiously, last point is not displayed.

  // ----------------------------------------------------------------------
  // User "interaction".
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();

  // Displays everything.
  viewer.show();

  return res ? 0 : 1;
}
