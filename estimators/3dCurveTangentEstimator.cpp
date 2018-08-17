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
 * @file 3dCurveTangentEstimator.cpp
 * @ingroup visualisationTools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 *
 * @date 2015/06/16
 *
 * This file is part of the DGtal library
 */

/**
 * Description of 3dCurveTangentEstimator <p>
 *
 * Display a 3D curve given as the <input> filename (with possibly
 * projections and/or tangent information) by using QGLviewer.
 */

#include <iostream>
#include <iterator>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/base/Exceptions.h"
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/geometry/curves/GridCurve.h"
#include "DGtal/geometry/curves/Naive3DDSSComputer.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/estimation/VoronoiCovarianceMeasure.h"
#include "DGtal/geometry/curves/estimation/LambdaMST3D.h"
#include "DGtal/geometry/curves/estimation/FunctorsLambdaMST.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"


using namespace DGtal;
using namespace std;




/**
 @page Doc3dCurveTangentEstimator 3dCurveTangentEstimator
 
 @brief This program estimates the tangent vector to a set of 3D integer points, which are supposed to approximate a 3D curve.

 @b Usage: ./estimators/3dCurveTangentEstimator [options] --input <filename>


This program estimates the tangent vector to a set of 3D integer points, which are supposed to approximate a 3D curve. This set of points is given as a list of points in file <input>.
The tangent estimator uses either the digital Voronoi Covariance Measure (VCM) or the 3D lambda-Maximal Segment Tangent (L-MST).
This program can also displays the curve and tangent estimations, and it can also extract maximal digital straight segments (2D and 3D).

@note It is not compulsory for the points to be ordered in sequence, except if you wish to compute maximal digital straight segments. In this case, you can select the connectivity of your curve between 6 (standard) or 26 (naive).


@b Allowed @b options are :

 @code
  -h [ --help ]                         display this message
  -i [ --input ] arg                    the name of the text file containing 
                                        the list of 3D points: (x y z) per 
                                        line.
  -V [ --view ] arg (=OFF)              toggles display ON/OFF
  -b [ --box ] arg (=0)                 specifies the the tightness of the 
                                        bounding box around the curve with a 
                                        given integer displacement <arg> to 
                                        enlarge it (0 is tight)
  -v [ --viewBox ] arg (=WIRED)         displays the bounding box, <arg>=WIRED 
                                        means that only edges are displayed, 
                                        <arg>=COLORED adds colors for planes 
                                        (XY is red, XZ green, YZ, blue).
  -T [ --connectivity ] arg (=6)        specifies whether it is a 6-connected 
                                        curve or a 26-connected curve: arg=6 | 
                                        26.
  -C [ --curve3d ]                      displays the 3D curve
  -c [ --curve2d ]                      displays the 2D projections of the 3D 
                                        curve on the bounding box
  -3 [ --cover3d ]                      displays the 3D tangential cover of the
                                        curve
  -2 [ --cover2d ]                      displays the 2D projections of the 3D 
                                        tangential cover of the curve
  -n [ --nbColors ] arg (=3)            sets the number of successive colors 
                                        used for displaying 2d and 3d maximal 
                                        segments (default is 3: red, green, 
                                        blue)
  -t [ --tangent ]                      displays the tangents to the curve
  -R [ --big-radius ] arg (=10)         the radius parameter R in the VCM 
                                        estimator.
  -r [ --small-radius ] arg (=3)        the radius parameter r in the VCM 
                                        estimator.
  -m [ --method ] arg (=VCM)            the method of tangent computation: VCM 
                                        (default), L-MST.
  -o [ --output ] arg (=3d-curve-tangent-estimations)
                                        the basename of the output text file 
                                        which will contain points and tangent 
                                        vectors: (x y z tx ty tz) per line

 @endcode

 @b Example: 
 This command line show an example of tangent estimation with the VCM estimator.
 @code
  3dCurveTangentEstimator -i ${DGtal}/examples/samples/sinus.dat -V ON -c -R 20 -r 3 -T 6
 @endcode


 You should obtain such a result:
 @image html res3dCurveTangentEstimator.png "Resulting tangent vectors (red) obtained with the VCM estimator."
 
 @see
 @ref 3dCurveTangentEstimator.cpp

 */


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
void displayAxes( Viewer3D<space, kspace> & viewer,
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
    viewer.setLineColor(AXIS_COLOR);
    viewer.addLine( DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p0[ 2 ]),
		    DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p0[ 2 ]),  AXIS_LINESIZE );
    viewer.addLine( DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p0[ 2 ]),
		    DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p0[ 2 ]),  AXIS_LINESIZE );
    viewer.addLine( DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p0[ 2 ]),
		    DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p1[ 2 ]),  AXIS_LINESIZE );
    viewer.addLine( DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p0[ 2 ]),
		    DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p0[ 2 ]),  AXIS_LINESIZE );
    viewer.addLine( DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p0[ 2 ]),
		    DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p1[ 2 ]),  AXIS_LINESIZE );
    viewer.addLine( DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p0[ 2 ]),
		    DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p0[ 2 ]),  AXIS_LINESIZE );
    viewer.addLine( DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p0[ 2 ]),
		    DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p1[ 2 ]),  AXIS_LINESIZE );
    viewer.addLine( DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p1[ 2 ]),
		    DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p1[ 2 ]),  AXIS_LINESIZE );
    viewer.addLine( DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p1[ 2 ]),
		    DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p1[ 2 ]),  AXIS_LINESIZE );
    viewer.addLine( DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p0[ 2 ]),
		    DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p1[ 2 ]),  AXIS_LINESIZE );
    viewer.addLine( DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p1[ 2 ]),
		    DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p1[ 2 ]),  AXIS_LINESIZE );
    viewer.addLine( DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p1[ 2 ]),
		    DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p1[ 2 ]),  AXIS_LINESIZE );
  }
  if ( mode == "COLORED" )
  {
    viewer.setFillColor(XY_COLOR);
    viewer.addQuad(DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p1[ 2 ]),
		   DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p1[ 2 ]),
		   DGtal::Z3i::RealPoint(p0[ 0 ], p0[ 1 ], p1[ 2 ]),
		   DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p1[ 2 ]) );
    viewer.setFillColor(XZ_COLOR);
    viewer.addQuad(DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p1[ 2 ]),
		   DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p1[ 2 ]),
		   DGtal::Z3i::RealPoint(p0[ 0 ], p1[ 1 ], p0[ 2 ]),
		   DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p0[ 2 ]));
    viewer.setFillColor(YZ_COLOR);
    viewer.addQuad(DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p1[ 2 ]),
		   DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p1[ 2 ]),
		   DGtal::Z3i::RealPoint(p1[ 0 ], p0[ 1 ], p0[ 2 ]),
		   DGtal::Z3i::RealPoint(p1[ 0 ], p1[ 1 ], p0[ 2 ]));
  }
}

template <typename KSpace, typename Naive3DDSSComputer, typename space, typename kspace >
void displayDSS3d( Viewer3D<space, kspace> & viewer,
		   const KSpace & ks, const Naive3DDSSComputer & dss3d,
		   const DGtal::Color & color3d )
{
  viewer << CustomColors3D( color3d, color3d ) << dss3d;
}

template <typename Point1, typename Point2>
void assign( Point1 & p1, const Point2 & p2 )
{
  p1[ 0 ] = p2[ 0 ];
  p1[ 1 ] = p2[ 1 ];
  p1[ 2 ] = p2[ 2 ];
}

template <typename KSpace, typename Naive3DDSSComputer, typename space, typename kspace >
void displayDSS3dTangent( Viewer3D<space, kspace> & viewer,
			  const KSpace & ks, const Naive3DDSSComputer & dss3d,
			  const DGtal::Color & color3d )
{
  typedef typename Naive3DDSSComputer::Point3d Point;
  typedef typename Naive3DDSSComputer::Integer Integer;
  typedef typename Naive3DDSSComputer::PointR3d PointR3d;
  typedef typename Display3D<>::BallD3D PointD3D;
  Point directionZ3;
  PointR3d interceptR, thicknessR;
  Z3i::RealPoint P1, P2, direction, intercept, thickness;  
  dss3d.getParameters( directionZ3, interceptR, thicknessR );
  
  intercept[0] = (double) NumberTraits<Integer>::castToDouble ( interceptR[0].first ) / (double) NumberTraits<Integer>::castToDouble ( interceptR[0].second );
  intercept[1] = (double) NumberTraits<Integer>::castToDouble ( interceptR[1].first ) / (double) NumberTraits<Integer>::castToDouble ( interceptR[1].second );
  intercept[2] = (double) NumberTraits<Integer>::castToDouble ( interceptR[2].first ) / (double) NumberTraits<Integer>::castToDouble ( interceptR[2].second );
  thickness[0] = (double) NumberTraits<Integer>::castToDouble ( thicknessR[0].first ) / (double) NumberTraits<Integer>::castToDouble ( thicknessR[0].second );
  thickness[1] = (double) NumberTraits<Integer>::castToDouble ( thicknessR[1].first ) / (double) NumberTraits<Integer>::castToDouble ( thicknessR[1].second );
  thickness[2] = (double) NumberTraits<Integer>::castToDouble ( thicknessR[2].first ) / (double) NumberTraits<Integer>::castToDouble ( thicknessR[2].second );
  
  assign( direction, directionZ3 );
  direction /= direction.norm();
  assign( P1, *dss3d.begin() );
  assign( P2, *(dss3d.end()-1) );
  double t1 = (P1 - intercept).dot( direction );
  double t2 = (P2 - intercept).dot( direction );
  
  Z3i::RealPoint Q1 = intercept + t1 * direction;
  Z3i::RealPoint Q2 = intercept + t2 * direction;
  viewer.setLineColor(color3d);
  viewer.addLine( Z3i::RealPoint(Q1[ 0 ]-0.5, Q1[ 1 ]-0.5, Q1[ 2 ]-0.5),
		  Z3i::RealPoint(Q2[ 0 ]-0.5, Q2[ 1 ]-0.5, Q2[ 2 ]-0.5),
		  MS3D_LINESIZE );
}

template <typename KSpace, typename StandardDSS6Computer, typename space, typename kspace >
void displayProj2d6( Viewer3D<space, kspace> & viewer,
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
    switch (i) 
    {
      case 0: q = Point3d( 2*b[ i ]  , 2*p[ 0 ]+1, 2*p[ 1 ]+1 ); break;
      case 1: q = Point3d( 2*p[ 0 ]+1, 2*b[ i ]  , 2*p[ 1 ]+1 ); break;
      case 2: q = Point3d( 2*p[ 0 ]+1, 2*p[ 1 ]+1, 2*b[ i ]   ); break;
    }
    Cell c = ks.uCell( q );
    viewer << CustomColors3D( color2d, color2d ) << c;
  }
    }
}

template <typename KSpace, typename StandardDSS6Computer, typename space, typename kspace >
void displayDSS2d6( Viewer3D<space, kspace> & viewer,
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
  typedef typename Display3D<>::BallD3D PointD3D;
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
      std::vector<PointD3D> bb;
      PointD3D p3;
      for ( unsigned int j = 0; j < pts2d.size(); ++j )
  {
    switch (i) 
    {
      case 0: p3.center[0] = (double) b[ i ]-0.5; p3.center[1] = pts2d[ j ][ 0 ];  p3.center[2] = pts2d[ j ][ 1 ]; break;
      case 1: p3.center[0] = pts2d[ j ][ 0 ];  p3.center[1] = (double) b[ i ]-0.5; p3.center[2] = pts2d[ j ][ 1 ];     break;
      case 2: p3.center[0] = pts2d[ j ][ 0 ];  p3.center[1] = pts2d[ j ][ 1 ];     p3.center[2] = (double) b[ i ]-0.5; break;
    }
    bb.push_back( p3 );
  }
      for ( unsigned int j = 0; j < pts2d.size(); ++j ){
  viewer.setLineColor(color2d);
  viewer.addLine( DGtal::Z3i::RealPoint(bb[ j ].center[0], bb[ j ].center[1], bb[ j ].center[2]),
                        DGtal::Z3i::RealPoint(bb[ (j+1)%4 ].center[0], bb[ (j+1)%4 ].center[1], bb[ (j+1)%4 ].center[2]),
      MS3D_LINESIZE );
      }
    } // for ( DGtal::Dimension i = 0; i < 3; ++i )
}


//Why not to just project 3D curve?
template <typename KSpace, typename Naive3DDSSComputer, typename space, typename kspace >
void displayProj2d26( Viewer3D<space, kspace> & viewer,
		    const KSpace & ks, const Naive3DDSSComputer & dss3d,
		    const DGtal::Color & color2d )
{
  typedef typename Naive3DDSSComputer::ArithmeticalDSSComputer2d ArithmeticalDSSComputer2d;
  typedef typename ArithmeticalDSSComputer2d::ConstIterator ConstIterator2d;
  typedef typename ArithmeticalDSSComputer2d::Point Point2d;
  typedef typename KSpace::Cell Cell;
  typedef typename KSpace::Point Point3d;
  Point3d b = ks.lowerBound();
  
  const ArithmeticalDSSComputer2d & dssXY = dss3d.arithmeticalDSS2dXY();
  const ArithmeticalDSSComputer2d & dssXZ = dss3d.arithmeticalDSS2dXZ();
  const ArithmeticalDSSComputer2d & dssYZ = dss3d.arithmeticalDSS2dYZ();
  
  bool validXY = dss3d.validArithmeticalDSS2d ( 2 );
  bool validXZ = dss3d.validArithmeticalDSS2d ( 1 );
  bool validYZ = dss3d.validArithmeticalDSS2d ( 0 );
  
  if ( validXY && validXZ ) 
  { //XY-plane, XZ-plane
    for ( ConstIterator2d itXY = dssXY.begin(), itXZ = dssXZ.begin(), itPEnd = dssXY.end(); itXY != itPEnd; ++itXY, ++itXZ )
    {
      Point2d p1 = *itXY, p2 = *itXZ;
      Point3d q1, q2, q3;
      q1 = Point3d ( 2*b[ 0 ], 2*p1[ 1 ]+1, 2*p2[ 1 ]+1 );
      q2 = Point3d ( 2*p2[ 0 ]+1, 2*b[ 1 ], 2*p2[ 1 ]+1 );
      q3 = Point3d ( 2*p1[ 0 ]+1, 2*p1[ 1 ]+1, 2*b[ 2 ] );
      Cell c1 = ks.uCell( q1 ); Cell c2 = ks.uCell( q2 ); Cell c3 = ks.uCell( q3 );
      viewer << CustomColors3D( color2d, color2d ) << c1;
      viewer << CustomColors3D( color2d, color2d ) << c2;
      viewer << CustomColors3D( color2d, color2d ) << c3;
    }
  }  
  else 
  {
    if ( validYZ && validXY ) 
    { //XY-plane, YZ-plane
      for ( ConstIterator2d itYZ = dssYZ.begin(), itXY = dssXY.begin(), itPEnd = dssXY.end(); itXY != itPEnd; ++itXY, ++itYZ )
      {
	Point2d p1 = *itYZ, p2 = *itXY;
	Point3d q1, q2, q3;
	q1 = Point3d ( 2*b[ 0 ], 2*p1[ 0 ]+1, 2*p1[ 1 ]+1 );
	q2 = Point3d ( 2*p2[ 0 ]+1, 2*b[ 1 ], 2*p1[ 1 ]+1 );
	q3 = Point3d ( 2*p2[ 0 ]+1, 2*p2[ 1 ]+1, 2*b[ 2 ] );
	Cell c1 = ks.uCell( q1 ); Cell c2 = ks.uCell( q2 ); Cell c3 = ks.uCell( q3 );
	viewer << CustomColors3D( color2d, color2d ) << c1;
	viewer << CustomColors3D( color2d, color2d ) << c2;
	viewer << CustomColors3D( color2d, color2d ) << c3;
      }
    } 
    else
    {
      for ( ConstIterator2d itYZ = dssYZ.begin(), itXZ = dssXZ.begin(), itPEnd = dssXZ.end(); itXZ != itPEnd; ++itXZ, ++itYZ )
      {
	Point2d p1 = *itYZ, p2 = *itXZ;
	Point3d q1, q2, q3;
	q1 = Point3d ( 2*b[ 0 ], 2*p1[ 0 ]+1, 2*p1[ 1 ]+1 );
	q2 = Point3d ( 2*p2[ 0 ]+1, 2*b[ 1 ], 2*p2[ 1 ]+1 );
	q3 = Point3d ( 2*p2[ 0 ]+1, 2*p1[ 0 ]+1, 2*b[ 2 ] );
	Cell c1 = ks.uCell( q1 ); Cell c2 = ks.uCell( q2 );Cell c3 = ks.uCell( q3 );
	viewer << CustomColors3D( color2d, color2d ) << c1;
	viewer << CustomColors3D( color2d, color2d ) << c2;
	viewer << CustomColors3D( color2d, color2d ) << c3;
      }
    }
  }
}


template <typename KSpace, typename Naive3DDSSComputer, typename space, typename kspace >
void displayDSS2d26( Viewer3D<space, kspace> & viewer,
		   const KSpace & ks, const Naive3DDSSComputer & dss3d,
		   const DGtal::Color & color2d )
{
  typedef typename Naive3DDSSComputer::ConstIterator ConstIterator3d;
  typedef typename Naive3DDSSComputer::ArithmeticalDSSComputer2d ArithmeticalDSSComputer2d;
  typedef typename ArithmeticalDSSComputer2d::ConstIterator ConstIterator2d;
  typedef typename ArithmeticalDSSComputer2d::Point Point2d;
  typedef typename KSpace::Cell Cell;
  typedef typename KSpace::Point Point3d;
  typedef DGtal::PointVector<2,double> PointD2d;
  typedef typename Display3D<>::BallD3D PointD3D;
  Point3d b = ks.lowerBound();
  for ( DGtal::Dimension i = 0; i < 3; ++i )
  {
    if ( !dss3d.validArithmeticalDSS2d ( i ) )
      continue;
    
    const typename ArithmeticalDSSComputer2d::Primitive & dss2d
    = dss3d.arithmeticalDSS2d( i ).primitive();
    // draw 2D bounding boxes for each arithmetical dss 2D.
    std::vector<PointD2d> pts2d;
    pts2d.push_back( dss2d.project(dss2d.back(), dss2d.Uf()) );
    pts2d.push_back( dss2d.project(dss2d.back(), dss2d.Lf()) );
    pts2d.push_back( dss2d.project(dss2d.front(), dss2d.Lf()) );
    pts2d.push_back( dss2d.project(dss2d.front(), dss2d.Uf()) );
    std::vector<PointD3D> bb;
    PointD3D p3;
    for ( unsigned int j = 0; j < pts2d.size(); ++j )
    {
      switch (i) 
      {
	case 0: p3.center[0] = (double) b[ i ]-0.5; p3.center[1] = pts2d[ j ][ 0 ];  p3.center[2] = pts2d[ j ][ 1 ]; break;
	case 1: p3.center[0] = pts2d[ j ][ 0 ];  p3.center[1] = (double) b[ i ]-0.5; p3.center[2] = pts2d[ j ][ 1 ];     break;
	case 2: p3.center[0] = pts2d[ j ][ 0 ];  p3.center[1] = pts2d[ j ][ 1 ];     p3.center[2] = (double) b[ i ]-0.5; break;
      }
      bb.push_back( p3 );
    }
    for ( unsigned int j = 0; j < pts2d.size(); ++j ){
      viewer.setLineColor(color2d);
      viewer.addLine( DGtal::Z3i::RealPoint(bb[ j ].center[0], bb[ j ].center[1], bb[ j ].center[2]),
		      DGtal::Z3i::RealPoint(bb[ (j+1)%4 ].center[0], bb[ (j+1)%4 ].center[1], bb[ (j+1)%4 ].center[2]),
		      MS3D_LINESIZE );
    }
  } // for ( DGtal::Dimension i = 0; i < 3; ++i )
}

/**
 * Displays the tangential cover of 6-connected curves (i.e. standard curves).
 */
template <typename KSpace, typename PointIterator, typename space, typename kspace >
bool displayCover6( Viewer3D<space, kspace> & viewer,
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

  viewer << SetMode3D( algo.className(), "BoundingBox" );
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
      Color color = cmap_hue( c );
      if ( tangent ) displayDSS3dTangent( viewer, ks, ms3d, color );
      if ( dss3d )   displayDSS3d( viewer, ks, ms3d, color );
      if ( dss2d )   displayDSS2d6( viewer, ks, ms3d, color );
      if ( proj2d )  displayProj2d6( viewer, ks, ms3d, CURVE2D_COLOR );
      c++;
    }
  return true;
}

/**
 * Displays the tangential cover of 26-connected curves (i.e. naive
 * curves). Note that is still experimental.
 */
template <typename KSpace, typename PointIterator, typename space, typename kspace >
bool displayCover26( Viewer3D<space, kspace> & viewer,
		   const KSpace & ks, PointIterator b, PointIterator e,
		   bool dss3d, bool proj2d, bool dss2d, bool tangent,
		   int nbColors )
{
  typedef typename PointIterator::value_type Point;
  typedef Naive3DDSSComputer<PointIterator,int, 8 > SegmentComputer;
  typedef SaturatedSegmentation<SegmentComputer> Decomposition;
  typedef typename Decomposition::SegmentComputerIterator SegmentComputerIterator;
  typedef typename SegmentComputer::ArithmeticalDSSComputer2d ArithmeticalDSSComputer2d;
  SegmentComputer algo;
  Decomposition theDecomposition(b, e, algo);
  
  viewer << SetMode3D( algo.className(), "BoundingBox" );
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
	 if ( dss2d )   displayDSS2d26( viewer, ks, ms3d, color );
	 if ( proj2d )  displayProj2d26( viewer, ks, ms3d, CURVE2D_COLOR );
	 c++;
       }
       return true;
}

template < typename PointIterator, typename Space, typename TangentSequence >
void ComputeVCM ( const double & R, const double & r,
		  const PointIterator & begin, const PointIterator & end, TangentSequence & tangents, const std::string & output )
{
  using namespace DGtal;
  typedef ExactPredicateLpSeparableMetric < Space, 2 > Metric; // L2-metric
  typedef VoronoiCovarianceMeasure < Space, Metric > VCM;
  typedef HyperRectDomain < Space > Domain;
  typedef functors::HatPointFunction < typename Space::Point, double > KernelFunction;
  typedef EigenDecomposition < 3, double > LinearAlgebraTool;
  typedef LinearAlgebraTool::Matrix Matrix;
  
  fstream outputStream;
  outputStream.open ( ( output + ".vcm" ).c_str(), std::ios::out );
  outputStream << "# VCM estimation R=" << R << " r=" << r << " chi=hat" << endl; 
   
  Metric l2;
  VCM vcm ( R, ceil( r ), l2, true );
  vcm.init( begin, end );
  KernelFunction chi( 1.0, r );
  Matrix vcm_r, evec;
  typename Space::RealVector eval;
  for ( PointIterator it = begin; it != end; ++it )
  {
    // Compute VCM and diagonalize it.
    vcm_r = vcm.measure( chi, *it );
    LinearAlgebraTool::getEigenDecomposition ( vcm_r, evec, eval );
    typename Space::RealVector tangent = evec.column( 0 );
    tangents.push_back ( tangent );
    outputStream << (*it)[0]   << " " << (*it)[1]   << " " << (*it)[2] << " "
    << tangent[0] << " " << tangent[1] << " " << tangent[2] << endl;
  }
  outputStream.close();
}

template < typename PointIterator, typename Space, typename TangentSequence >
void ComputeLMST6 ( const PointIterator & begin, const PointIterator & end, TangentSequence & tangents, const std::string & output  )
{
  typedef typename PointIterator::value_type Point;
  typedef StandardDSS6Computer<PointIterator,int, 4 > SegmentComputer;
  typedef SaturatedSegmentation<SegmentComputer> Decomposition;
  typedef typename Decomposition::SegmentComputerIterator SegmentComputerIterator;
  typedef typename SegmentComputer::ArithmeticalDSSComputer2d ArithmeticalDSSComputer2d;
  SegmentComputer algo;
  Decomposition theDecomposition ( begin, end, algo );
  
  fstream outputStream;
  outputStream.open ( ( output + ".lmst" ).c_str(), std::ios::out );
  outputStream << "X Y Z X Y Z" << endl;
  
  LambdaMST3D < Decomposition > lmst64;
  lmst64.attach ( theDecomposition );
  lmst64.init ( begin, end );
  lmst64.eval ( begin, end, std::back_inserter ( tangents ) );
  typename TangentSequence::iterator itt = tangents.begin();
  for ( PointIterator it = begin; it != end; ++it, ++itt )
  {
    typename Space::RealVector & tangent = (*itt);
    if ( tangent.norm() != 0. )
      tangent = tangent.getNormalized();
    outputStream << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << " "
    << tangent[0] << " " << tangent[1] << " " << tangent[2] << endl;
  }
}

template < typename PointIterator, typename Space, typename TangentSequence >
void ComputeLMST26 ( const PointIterator & begin, const PointIterator & end, TangentSequence & tangents, const std::string & output  )
{
  typedef typename PointIterator::value_type Point;
  typedef Naive3DDSSComputer<PointIterator,int, 8 > SegmentComputer;
  typedef SaturatedSegmentation<SegmentComputer> Decomposition;
  typedef typename Decomposition::SegmentComputerIterator SegmentComputerIterator;
  typedef typename SegmentComputer::ArithmeticalDSSComputer2d ArithmeticalDSSComputer2d;
  SegmentComputer algo;
  Decomposition theDecomposition ( begin, end, algo );
  
  fstream outputStream;
  outputStream.open ( ( output + ".lmst" ).c_str(), std::ios::out );
  outputStream << "X Y Z X Y Z" << endl;
  
  LambdaMST3D < Decomposition > lmst64;
  lmst64.attach ( theDecomposition );
  lmst64.init ( begin, end );
  lmst64.eval ( begin, end, std::back_inserter ( tangents ) );
  typename TangentSequence::iterator itt = tangents.begin();
  for ( PointIterator it = begin; it != end; ++it, ++itt )
  {
    typename Space::RealVector & tangent = (*itt);
    if ( tangent.norm() != 0. )
      tangent = tangent.getNormalized();
    outputStream << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << " "
    << tangent[0] << " " << tangent[1] << " " << tangent[2] << endl;
  }
}

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

/**
 * Main function.
 * 
 * @param argc the number of parameters given on the line command.
 * 
 * @param argv an array of C-string, such that argv[0] is the name of
 * the program, argv[1] the first parameter, etc.
 */
int main(int argc, char **argv)
{
  using namespace std;
  using namespace DGtal;
  typedef SpaceND<3,int>                         Z3;
  typedef KhalimskySpaceND<3,int>                K3;
  typedef Z3::Point                              Point;
  typedef Z3::RealPoint                          RealPoint;
  typedef Z3::RealVector                         RealVector;
  typedef HyperRectDomain<Z3>                    Domain;
  typedef typename vector<Point>::const_iterator PointIterator;
  
  // specify command line ----------------------------------------------
  QApplication application(argc,argv); // remove Qt arguments.
  po::options_description general_opt("Specific allowed options are ");
  general_opt.add_options()
  ("help,h", "display this message")
  ("input,i", po::value<string>(), "the name of the text file containing the list of 3D points: (x y z) per line." )
  ("view,V", po::value<string>()->default_value( "OFF" ), "toggles display ON/OFF" )
  ("box,b",  po::value<int>()->default_value( 0 ), "specifies the the tightness of the bounding box around the curve with a given integer displacement <arg> to enlarge it (0 is tight)" )
  ("viewBox,v",  po::value<string>()->default_value( "WIRED" ), "displays the bounding box, <arg>=WIRED means that only edges are displayed, <arg>=COLORED adds colors for planes (XY is red, XZ green, YZ, blue)." )
  ("connectivity,T", po::value<string>()->default_value( "6" ), "specifies whether it is a 6-connected curve or a 26-connected curve: arg=6 | 26.")
  ("curve3d,C", "displays the 3D curve")
  ("curve2d,c", "displays the 2D projections of the 3D curve on the bounding box")
  ("cover3d,3", "displays the 3D tangential cover of the curve" )
  ("cover2d,2", "displays the 2D projections of the 3D tangential cover of the curve" )
  ("nbColors,n",  po::value<int>()->default_value( 3 ), "sets the number of successive colors used for displaying 2d and 3d maximal segments (default is 3: red, green, blue)" )
  ("tangent,t", "displays the tangents to the curve" )
  ("big-radius,R", po::value<double>()->default_value( 10.0 ), "the radius parameter R in the VCM estimator." )
  ("small-radius,r", po::value<double>()->default_value( 3.0 ), "the radius parameter r in the VCM estimator." )
  ("method,m", po::value<string>()->default_value( "VCM" ), "the method of tangent computation: VCM (default), L-MST." )
  ("output,o", po::value<string>()->default_value( "3d-curve-tangent-estimations"), "the basename of the output text file which will contain points and tangent vectors: (x y z tx ty tz) per line" )
  ;
  po::positional_options_description pos_opt;
  pos_opt.add("input", 1);
  
  // parse command line ----------------------------------------------
  bool parseOK=true;
  po::variables_map vm;
  try 
  {
    po::command_line_parser clp( argc, argv );
    clp.options( general_opt ).positional( pos_opt );
    po::store( clp.run(), vm );
  } 
  catch( const exception& ex ) 
  {
    parseOK = false;
    trace.info() << "Error checking program options: "<< ex.what() << endl;
  }
  po::notify( vm );
  if( !parseOK || vm.count("help")||argc<=1)
  {
    cout << "Usage: " << argv[0] << " [options] --input  <filename>\n"
    << "This program estimates the tangent vector to a set of 3D integer points, which are supposed to approximate a 3D curve. "
    << "This set of points is given as a list of points in file <input>."
    << endl 
    << "The tangent estimator uses either the digital Voronoi Covariance Measure (VCM) or the 3D lambda-Maximal Segment Tangent (L-MST)." 
    << endl
    << "This program can also displays the curve and tangent estimations, "
    << "and it can also extract maximal digital straight segments (2D and 3D)."
    << endl
    << "Note: It is not compulsory for the points to be ordered in sequence, except if you wish to compute maximal digital straight segments." 
    << endl
    << "In this case, you can select the connectivity of your curve between 6 (standard) or 26 (naive).\n"
    << general_opt << "\n\n";
    cout << "Example:\n"
    << "3dCurveTangentEstimator -i ${DGtal}/examples/samples/sinus.dat -V ON -c -R 20 -r 3 -T 6\n";
    return 0;
  }
  
  // process command line ----------------------------------------------
  string input = vm["input"].as<string>();
  int b = vm["box"].as<int>();
  // Create curve 3D.
  vector<Point> sequence;
  fstream inputStream;
  inputStream.open ( input.c_str(), ios::in);
  try 
  {
    sequence = PointListReader<Point>::getPointsFromInputStream( inputStream );
    if ( sequence.size() == 0) throw IOException();
  }
  catch (DGtal::IOException & ioe) 
  {
    trace.error() << "Size is null." << endl;
  }
  inputStream.close();
  
  // start viewer
  bool view = vm[ "view" ].as<string>() == "ON";
  
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
  upperBound += Point::diagonal( b + 1 );
  K3 ks; ks.init( lowerBound, upperBound, true );
  Viewer3D<Z3,K3> viewer( ks );
  trace.beginBlock ( "Tool 3dCurveTangentEstimator" );
  
  std::vector< RealVector > tangents;
  string output = vm["output"].as<string>();
  string method = vm["method"].as<string>();
  
  if ( method == "VCM" )
  {
    // input points of the curve are in sequence vector.
    const double R = vm["big-radius"].as<double>();
    trace.info() << "Big radius   R = " << R << endl;
    const double r = vm["small-radius"].as<double>();
    trace.info() << "Small radius r = " << r << endl;
    
    ComputeVCM < PointIterator, Z3, std::vector< RealVector > > ( R, r, sequence.begin(), sequence.end(), tangents, output );
  }
  else if ( method == "L-MST" )
  {
    if (vm[ "connectivity" ].as<string>() == "6")
      ComputeLMST6  < PointIterator, Z3, std::vector< RealVector > > ( sequence.begin(), sequence.end(), tangents, output );
    else
      ComputeLMST26  < PointIterator, Z3, std::vector< RealVector > > ( sequence.begin(), sequence.end(), tangents, output );
  }
  else
  {
    trace.info() << "Wrong method! Try: VCM or L-MST" << endl;
    abort();
  }
  
  if ( view ) 
  {
    for ( unsigned int i = 0; i < tangents.size(); i++ )
    {
      // Display normal
      RealPoint p = sequence[i]; 
      RealVector tangent = tangents[i]; 
      viewer.setFillColor( Color ( 255, 0, 0, 255 ) );
      viewer.setLineColor( Color ( 255, 0, 0, 255 ) );
      viewer.addLine( p + 2.0 * tangent, p - 2.0 * tangent,  5.0 );
      viewer.setFillColor( Color ( 100, 100, 140, 255 ) );
      viewer.setLineColor( Color ( 100, 100, 140, 255 ) );
      viewer.addBall( p, 0.125, 8 );
    }
  }
  
  GridCurve<K3> gc( ks );
  try 
  {
    gc.initFromPointsVector( sequence );
  } 
  catch (DGtal::ConnectivityException& /*ce*/) 
  {
    trace.warning() << "[ConnectivityException] GridCurve only accepts a sequence of face adjacent points. Try connectivity=6 instead." << endl;
  }
  
  // ----------------------------------------------------------------------
  // Displays everything.
  bool res = true;
  if ( view )
  {
    viewer.show();
    // Display axes.
    if ( vm.count( "viewBox" ) )
      displayAxes<Point,RealPoint, Z3i::Space, Z3i::KSpace>( viewer, lowerBound, upperBound, vm[ "viewBox" ].as<std::string>() );
    // Display 3D tangential cover.
    res = vm[ "connectivity" ].as<string>() == "6"
    ? displayCover6( viewer, ks, sequence.begin(), sequence.end(),
		     vm.count( "cover3d" ),
		     vm.count( "curve2d" ),
		     vm.count( "cover2d" ),
		     vm.count( "tangent" ),
		     vm["nbColors"].as<int>() )
    : displayCover26( viewer, ks, sequence.begin(), sequence.end(),
		      vm.count( "cover3d" ),
		      vm.count( "curve2d" ),
		      vm.count( "cover2d" ),
		      vm.count( "tangent" ),
		      vm["nbColors"].as<int>() );
    // Display 3D curve points.
    if ( vm.count( "curve3d" ) )
    {
      viewer << CustomColors3D( CURVE3D_COLOR, CURVE3D_COLOR );
      for ( vector<Point>::const_iterator it = sequence.begin(), itE = sequence.end();
	   it != itE; ++it )
	   viewer << *it;  
    }
    // ----------------------------------------------------------------------
    // User "interaction".
    viewer << Viewer3D<Z3,K3>::updateDisplay;
    application.exec();
  }
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  
  return res ? 0 : 1;
}
