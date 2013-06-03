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
 * @file patternTriangulation.cpp
 * @ingroup visualisationTools
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2013/05/29
 *
 * This file is part of the DGtal library
 */

/**
 * Description of patternTriangulation <p>
 *
 * Display the Delaunay triangulation of a digital pattern
 * described by its rational slope.
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
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/arithmetic/LighterSternBrocot.h"
#include "DGtal/arithmetic/IntegerComputer.h"
#include "DGtal/arithmetic/Pattern.h"

#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/Color.h"

using namespace DGtal;
using namespace std;

namespace po = boost::program_options;
///////////////////////////////////////////////////////////////////////////////

/** 
 * A model of point functor that 
 * maps an odd convergent to a convex hull vertex
 *
 * @tparam Point a model of point
 */
template <typename Point>
struct OddConvexHullMap {

public: 
  typedef Point Vector; 

public: 
/** 
 * Constructor. 
 *
 * @param aShift shift vector: 
 * (1,-1) for standard pattern by default
 */
  OddConvexHullMap(const Vector& aShift = Vector(1,-1))
    : myS(aShift) {}
  
/** 
 * Main method 
 *
 * @param aP point to transform (convergent)
 * @return the transformed point (convex hull vertex)
 */
  Point operator()(const Point& aP) const
  {
    return aP + myS; 
  }

private: 
  Vector myS; 

}; 

/** 
 * A model of point functor that 
 * maps an odd convergent to a convex hull vertex
 *
 * @tparam Point a model of point
 */
template <typename Point>
struct EvenConvexHullMap {

public: 
  typedef Point Vector; 

public: 
/** 
 * Constructor. 
 *
 * @param aPoint last convergent
 * @param aShift shift vector:
 *  (1,-1) for standard pattern by default
 */
  EvenConvexHullMap(const Point& aPoint, const Vector& aShift = Vector(1,-1))
    : myZn(aPoint), myS(aShift) {}
  
/** 
 * Main method 
 *
 * @param aP point to transform (convergent)
 * @return the transformed point (convex hull vertex)
 */
  Point operator()(const Point& aP) const
  {
    return myZn - aP + myS; 
  }

private: 
  Point myZn; 
  Vector myS; 

}; 

///////////////////////////////////////////////////////////////////////////////
/** 
 * Procedure that displays an edge described by 
 * its two extremities. 
 *
 * @param aBoard board to display
 * @param aP first end of the segment
 * @param aQ second end of the segment
 * 
 * @tparam Board a model of 2d board
 * @tparam Point a model of point
 */
template <typename Board, typename Point>
void drawSegment(Board& aBoard, const Point& aP, const Point& aQ)
{
  aBoard.drawLine(aP[0], aP[1], aQ[0], aQ[1]);   
  aBoard << aP << aQ; 
}

/** 
 * Procedure that displays the edges of a polyline
 * described by a range of point.
 * 
 * @param aBoard board to display
 * @param itb begin iterator of the points range
 * @param ite end iterator of the points range 
 * @param aF a functor from a point to another point
 * 
 * @tparam Board a model of 2d board
 * @tparam Iterator at least a model of forward iterator
 * @tparam Functor a model of unary functor mapping 
 * points to points
 */
template <typename Board, typename Iterator, typename Functor>
void drawSegments(Board& aBoard, 
		  const Iterator& itb, const Iterator& ite, 
		  const Functor& aF)
{
  Iterator it = itb; 
  if (it != ite) 
    {
      Iterator prevIt = it; 
      ++it; 
      for ( ; it != ite; ++prevIt, ++it) 
	{
	  drawSegment( aBoard, aF(*prevIt), aF(*it) ); 
	}
    }
}


///////////////////////////////////////////////////////////////////////////////
/** 
 * Procedure that displays the convex hull
 * and possibly the farthest-point Delaunay 
 * triangulation of a given pattern.
 * 
 * @param aP irreductible pattern
 * @param aD number of repetitions of the pattern
 * @param withFDT boolean equal to true to 
 * display the internal edges of the 
 * farthest-point Delaunay triangulation
 *
 * @tparam Pattern a model of pattern
 * @tparam Integer a model of integer
 */
template <typename Pattern, typename Integer>
void displayConvexHull(const Pattern& aP, const Integer& aD, bool withFDT = false)
{
  typedef typename Pattern::Fraction Fraction; 
  typedef typename Pattern::Quotient Quotient; 
  typedef typename Pattern::Vector2I Vector;  
  typedef std::vector<Vector> Convergents; 
  typedef typename Pattern::Point2I Point;  
  
  Board2D aBoard;

  //list of convergents
  Convergents oddConvergents, evenConvergents;

  //fill the lists
  Fraction zn = aP.slope();
  // aD > 1
  if (aD > 1)
    {
      Fraction znm1 = zn.father(); 
      znm1.selfDisplay(std::cout);   std::cout << std::endl;  
      if (zn.odd())
	oddConvergents.push_back(Vector(znm1.q(),znm1.p())); 
      else 
	evenConvergents.push_back(Vector(znm1.q(),znm1.p())); 
    }
  // aD >= 1
  for (Quotient i = 1; i <= zn.k(); ++i)
    {
      Fraction zk = zn.reduced(i); 
      zk.selfDisplay(std::cout); 
      std::cout << " " << zn.k() - i << std::endl;
      if (((zn.k() - i)%2) == 1 )
	{ //odd
	  oddConvergents.push_back(Vector(zk.q(),zk.p())); 
	}
      else 
	{
	  evenConvergents.push_back(Vector(zk.q(),zk.p())); 
	}
    }

  //display the domain
  Point Zn = Point(aD*zn.q(),aD*zn.p()); 
  aBoard << Z2i::Domain( Point(0,0), Zn ); 
  aBoard << SetMode( Zn.className(), "Grid" );
  aBoard.setLineStyle(Board2D::Shape::SolidStyle); 
  aBoard.setPenColor(DGtal::Color::Black);

  //display the two main lists
  OddConvexHullMap<Point> oddH;
  drawSegments( aBoard, oddConvergents.begin(), oddConvergents.end(), oddH ); 

  EvenConvexHullMap<Point> evenH(Zn);
  drawSegments( aBoard, evenConvergents.begin(), evenConvergents.end(), evenH ); 

  //display the four last segments
  drawSegment( aBoard, Point(0,0), Zn ); 
  if (evenConvergents.size() != 0)
    {
      drawSegment( aBoard, evenH(*evenConvergents.rbegin()), Zn ); 
      if (oddConvergents.size() != 0)
	{
	  drawSegment( aBoard, Point(0,0), oddH(*oddConvergents.rbegin()) );
	  drawSegment( aBoard, 
		       oddH(*oddConvergents.begin()), 
		       evenH(*evenConvergents.begin()) ); 
	}
      else 
	{
	  drawSegment( aBoard, Point(0,0), evenH(*evenConvergents.rbegin()) );
	}
    }

  if (withFDT)
    {//display internal edges

      //first internal edge
      if (evenConvergents.size() != 0)
	drawSegment( aBoard, 
		     Point(0,0), 
		     evenH(*evenConvergents.rbegin()) ); 
      //other internal edges
      typedef typename Convergents::const_reverse_iterator RI; 
      RI oddRit = oddConvergents.rbegin(); 
      RI evenRit = evenConvergents.rbegin(); 
      bool hasToContinue = true; 
      while (hasToContinue) 
	{
	  if (oddRit != oddConvergents.rend())
	    {	      
	      drawSegment( aBoard, 
			   oddH(*oddRit), 
			   evenH(*evenRit) ); 
	    }
	  else 
	    hasToContinue = false; 
	  ++evenRit; 

	  if ( (hasToContinue) && (evenRit != evenConvergents.rend()) )
	    {	      
	      drawSegment( aBoard, 
			   oddH(*oddRit), 
			   evenH(*evenRit) ); 
	    }
	  else 
	    hasToContinue = false; 
	  ++oddRit; 
	}

      aBoard.saveEPS("FDT.eps");
    }
  else
    {
      aBoard.saveEPS("CH.eps");
    }  
}

///////////////////////////////////////////////////////////////////////////////
// /** 
//  * //TO REMOVE 
//  * Procedure that displays the main facet
//  * of a pattern.  
//  *
//  * @param aBoard board to display
//  * @param aP irreductible pattern
//  * 
//  * @tparam Board a model of 2d board
//  * @tparam Pattern a model of patter
//  */
// template <typename Board, typename Pattern>
// void drawMainFacet(Board& aBoard, const Pattern& aP)
// {
//   typedef typename Pattern::Fraction Fraction; 
//   typedef typename Pattern::Quotient Quotient; 
//   typedef typename Pattern::Point2I Point;  
//   typedef typename Pattern::Vector2I Vector;  

//   Point O = aP.U(NumberTraits<Quotient>::ZERO); 
//   Point Zn = aP.U(NumberTraits<Quotient>::ONE); 
//   Vector B = aP.bezout(); 

//   aBoard.drawTriangle( O[0], O[1], Zn[0], Zn[1], B[0], B[1] ); 
// }

/** 
 * Procedure that displays a triangle
 *
 * @param aBoard board to display
 * @param aP a first point
 * @param aQ a second point
 * @param aR a third point
 * 
 * @tparam Board a model of 2d board
 * @tparam Point a model of point
 */
template <typename Board, typename Point>
void drawTriangle(Board& aBoard, 
		  const Point& aP, const Point& aQ, const Point& aR)
{
  aBoard.drawTriangle( aP[0], aP[1], aQ[0], aQ[1], aR[0], aR[1] ); 
}

//pb mySlope is private et non protected...
// template <typename TFraction>
// struct MyPattern : public DGtal::Pattern<TFraction> 
// {
// public: 
//   typedef TFraction Fraction;
//   typedef MyPattern<TFraction> Self;
//   typedef typename Fraction::Integer Integer;
//   typedef typename Fraction::Quotient Quotient;

//   typedef IntegerComputer<Integer> IC;
//   typedef typename IC::Point2I Point2I;
//   typedef typename IC::Vector2I Vector2I;

// public: 
//   MyPattern( Integer p, Integer q ) : mySlope(p, q); 
//   MyPattern( Fraction f = Fraction( 0, 0 ) ) : mySlope(f); 

// public: 
//   Vector2I positiveBezout() const
//   {
//     bezout(); 
//   }
// }; 

// //EN FAIRE DES FONCTEURS
// template<typename Pattern>
// typename Pattern::Vector2I
// getPositiveBezout(const Pattern& aP)
// {
//   return aP.bezout(); 
// }

// template<typename Pattern>
// typename Pattern::Vector2I
// getNegativeBezout(const Pattern& aP)
// {
//   return aP.slope().even() 
//     ? aP.U( 1 ) - aP.previousPattern().U( 1 ) 
//     : aP.previousPattern().U( 1 );
// }

/** 
 * A helper class that provides two methods 
 * that returns positive or negative Bezout
 * point of a pattern. 
 * The positive Bezout point lies on the left side
 * whereas the negative Bezout point lies on the
 * right side of the straight line that has the slope
 * of the pattern. 
 */
struct InDirectBezoutComputer
{
  /** 
   * @return the positive Bezout point 
   * 
   * @param aP any pattern
   *
   * @tparam Pattern a model of pattern
   */
  template<typename Pattern>
  typename Pattern::Vector2I
  getPositiveBezout(const Pattern& aP) const
  {
    return aP.bezout(); 
  }

  /** 
   * @return the negative Bezout point 
   * 
   * @param aP any pattern
   *
   * @tparam Pattern a model of pattern
   */
  template<typename Pattern>
  typename Pattern::Vector2I
  getNegativeBezout(const Pattern& aP) const
  {
    return aP.slope().even() 
      ? aP.U( 1 ) - aP.previousPattern().U( 1 ) 
      : aP.previousPattern().U( 1 );
  }

}; 
 
/** 
 * A helper class that provides two methods 
 * that returns positive or negative Bezout
 * point of a pattern. 
 * The positive Bezout point lies on the right side
 * whereas the negative Bezout point lies on the
 * left side of the straight line that has the slope
 * of the pattern. 
 */
struct DirectBezoutComputer
{
  /** 
   * @return the positive Bezout point 
   * 
   * @param aP any pattern
   *
   * @tparam Pattern a model of pattern
   */
  template<typename Pattern>
  typename Pattern::Vector2I
  getPositiveBezout(const Pattern& aP) const
  {
    return aP.slope().even() 
      ? aP.U( 1 ) - aP.previousPattern().U( 1 ) 
      : aP.previousPattern().U( 1 );
  }

  /** 
   * @return the negative Bezout point 
   * 
   * @param aP any pattern
   *
   * @tparam Pattern a model of pattern
   */
  template<typename Pattern>
  typename Pattern::Vector2I
  getNegativeBezout(const Pattern& aP) const
  {
    return aP.bezout(); 
  }

}; 

/** 
 * Procedure that displays the upper part
 * of the closest-point Delaunay triangulation
 * of a given pattern.
 * 
 * @param aBoard board on which the triangulation
 * is displayed
 * @param aP irreductible pattern
 * @param aStartingPoint first point of the pattern
 * 
 * @tparam Board a model of board
 * @tparam Pattern a model of pattern
 * @tparam BezoutComputer a type having
 * getPositiveBezout and getNegativeBezout methods
 * in order to adapt patterns. 
 */
template <typename Board, typename Pattern, typename BezoutComputer>
void displayPartialCDT(Board& aBoard, 
		       const Pattern& aP, 
		       const typename Pattern::Point2I& aStartingPoint,
		       const BezoutComputer& aBC )
{
  typedef typename Pattern::Fraction Fraction; 
  typedef typename Pattern::Vector2I Vector;  
  typedef typename Pattern::Point2I Point;  
  
  Fraction f = aP.slope();
  // std::cout << f.p() << "/" << f.q() << std::endl;  
  if ( (f.p() >= 1)&&(f.q() >= 1) )
    {
      Point O = aStartingPoint; 
      Point Zn = aStartingPoint + aP.v();
      Vector v1 = aBC.getNegativeBezout(aP); 
      Point B = aStartingPoint + v1; 
      drawTriangle(aBoard, O, Zn, B); 

      if ( (f.p() > 1)||(f.q() > 1) )
	{
	  //recursive calls 
	  // - first pattern
	  displayPartialCDT(aBoard, Pattern(Fraction(v1[1],v1[0])), O, aBC); 
	  // - second pattern 
	  Vector v2 = aBC.getPositiveBezout(aP); 
	  displayPartialCDT(aBoard, Pattern(Fraction(v2[1],v2[0])), B, aBC); 
	}
    }
  // else 
  //   std::cout << " stop " << std::endl;   
}

/** 
 * Procedure that displays the upper part
 * of the closest-point Delaunay triangulation
 * of a given pattern.
 * 
 * @param aP irreductible pattern
 * @param aD number of repetitions of the pattern
 * 
 * @tparam Pattern a model of pattern
 * @tparam Integer a model of integer
 */
template <typename Pattern, typename Integer>
void displayPartialCDT(const Pattern& aP, const Integer& aD)
{
  std::cout << "partial CDT" << std::endl; 
  
  Board2D aBoard; 

  //display the domain
  typedef typename Pattern::Point2I Point;  
  typedef typename Pattern::Quotient Quotient;  
  typedef typename Pattern::Fraction Fraction;  
  Fraction zn = aP.slope(); 
  Point Zn = Point(aD*zn.q(),aD*zn.p()); 
  aBoard << Z2i::Domain( Point(0,0), Zn ); 
  aBoard << SetMode( Zn.className(), "Grid" );
  aBoard.setLineStyle(Board2D::Shape::SolidStyle); 
  aBoard.setPenColor(DGtal::Color::Black);

  //display the upper part of the triangulation
  for (Integer i = 0; i < aD; ++i)
    {
      displayPartialCDT( aBoard, aP, Point(i*zn.q(), i*zn.p()), InDirectBezoutComputer() );
    }

  //display the lower parts of the triangulation
  // - getting the reversed patterns
  typedef typename Pattern::Vector2I Vector;  
  typedef std::vector<Vector> Convergents; 
  Convergents oddConvergents, evenConvergents;
  // TO FACTORIZE
  // // aD > 1
  // if (aD > 1)
  //   {
  //     Fraction znm1 = zn.father(); 
  //     znm1.selfDisplay(std::cout);   std::cout << std::endl;  
  //     if (zn.odd())
  // 	oddConvergents.push_back(Vector(znm1.q(),znm1.p())); 
  //     else 
  // 	evenConvergents.push_back(Vector(znm1.q(),znm1.p())); 
  //   }
  // aD >= 1
  for (Quotient i = 1; i <= zn.k(); ++i)
    {
      Fraction zk = zn.reduced(i); 
      zk.selfDisplay(std::cout); 
      std::cout << " " << zn.k() - i << std::endl;
      if (((zn.k() - i)%2) == 1 )
	{ //odd
	  oddConvergents.push_back(Vector(zk.q(),zk.p())); 
	}
      else //even
	{
	  evenConvergents.push_back(Vector(zk.q(),zk.p())); 
	}
    }
  // - displaying the reversed patterns...
  //  - of odd slope: 
  Point rPStartingPoint; 
  OddConvexHullMap<Point> oddH;
  {
    typedef typename Convergents::const_reverse_iterator Iterator; 
    Iterator itb = evenConvergents.rbegin(); 
    Iterator ite = evenConvergents.rend(); 
    Point aStartingPoint = oddH( *oddConvergents.rbegin() ); 
    //
    Point p = aStartingPoint; 
    std::cout << p << std::endl; 
    Iterator it = itb; 
    if (it != ite)
      {
	for (++it; it != ite; ++it) 
	  {
	    std::cout << p << "(" << (*it)[1] << ", " <<  (*it)[0] << ")" << std::endl; 
	    displayPartialCDT( aBoard, Pattern( (*it)[1], (*it)[0] ), p, DirectBezoutComputer() ); 
	    p += *it; 
	  }
      }
    rPStartingPoint = p; 
  }
  //  - of even slope: 
  EvenConvexHullMap<Point> evenH(Zn);
  {
    typedef typename Convergents::const_reverse_iterator Iterator; 
    Iterator itb = oddConvergents.rbegin(); 
    Iterator ite = oddConvergents.rend(); 
    Point aStartingPoint = evenH( *evenConvergents.rbegin() ); 
    //
    Point p = aStartingPoint;
    std::cout << p << std::endl; 
    Iterator it = itb; 
    for ( ; it != ite; ++it) 
      {
	std::cout << p << std::endl; 
	p -= *it; 
	std::cout << p << "(" << (*it)[1] << ", " <<  (*it)[0] << ")" << std::endl; 
	displayPartialCDT( aBoard, Pattern( (*it)[1], (*it)[0] ), p, DirectBezoutComputer() ); 
      }
  }
  //  - of same slope: 
  {
    for (Integer i = 0; i < (aD-1); ++i)
      {
	Point p = rPStartingPoint + i*aP.v(); 
	displayPartialCDT( aBoard, aP, p, DirectBezoutComputer() );
      }
  }

 
  aBoard.saveEPS("CDT.eps");
}

// template <typename Board, typename Point, typename Iterator>
// void displayReversedPatterns(Board& aBoard, 
// 			     const Point& aStartingPoint, 
// 			     const Iterator& itb, const Iterator& ite) 
// {
//   Point p = aStartingPoint; 
//   Iterator it = itb; 
//   for ( ; it != ite; ++it) 
//     {
//       displayPartialCDT( aBoard, Pattern( (*it)[1], (*it)[2] ), p, DirectBezoutComputer() ); 
//       p += *it; 
//     }
// }

///////////////////////////////////////////////////////////////////////////////
/** 
 * Missing parameter error message.
 * 
 * @param param 
 */
void missingParam(std::string param)
{
  trace.error() << " Parameter: " << param << " is required...";
  trace.info() << std::endl;
  exit(1);
}


///////////////////////////////////////////////////////////////////////////////
/**
   Main function.

   @param argc the number of parameters given on the line command.

   @param argv an array of C-string, such that argv[0] is the name of
   the program, argv[1] the first parameter, etc.
*/
int main(int argc, char **argv)
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("aparam,a",   po::value<int>(), "pattern a parameter" )
    ("bparam,b",   po::value<int>(), "pattern b parameter" )
    ("delta,d",   po::value<int>()->default_value(1), "number of repetitions" )
    ("triangulation,t",   po::value<string>()->default_value("CH"), "output:\n\tClosest-point Delaunay triangulation {CDT}\n\tFarthest-point Delaunay triangulation {FDT}\n\tConvex hull {CH}" );

  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
 
  po::notify(vm);    
  if(!parseOK || vm.count("help")||argc<=1)
    {
      trace.info()<< "Draw the Delaunay triangulation of a pattern using DGtal library" <<std::endl 
		  << "Basic usage: "<<std::endl
		  << "\t" << argv[0] << " -a 5 -b 8 "<<std::endl
		  << general_opt << "\n";
      return 0;
    }


  if (!(vm.count("aparam"))) missingParam("-a");
  if (!(vm.count("bparam"))) missingParam("-b");

  typedef DGtal::int32_t Integer;
  Integer a = vm["aparam"].as<int>();
  Integer b = vm["bparam"].as<int>();
  if ( (a < 0) || (b <= 0) )
    {
      trace.error() << " parameters a and b should be strictly positive.";
      trace.info() << std::endl;
      return 0;
    }
  if (a >= b)
    {
      trace.error() << " parameter a should be strictly less than b.";
      trace.info() << std::endl;
      return 0;
    }
  Integer d = vm["delta"].as<int>(); 
  if (d <= 0)
    {
      trace.error() << " parameter d should be strictly positive";
      trace.info() << std::endl;
      return 0;
    }
  trace.info() << "a=" << a << ", b=" << b << ", d=" << d << std::endl; 

  typedef DGtal::int32_t Quotient;
  typedef LighterSternBrocot<Integer, Quotient, StdMapRebinder> SB;
  typedef SB::Fraction Fraction; // the type for fractions
  typedef Pattern<Fraction> Pattern; // the type for patterns
 
  Pattern pattern( a, b );
  //MyPattern essai(a, b); 

  string type = vm["triangulation"].as<string>();
  if (type == "CDT")
    {
      displayPartialCDT(pattern, d); 
    }
  else if (type == "FDT")
    {
      displayConvexHull(pattern, d, true); 
    }
  else if (type == "CH")
    {
      displayConvexHull(pattern, d); 
    }
  else 
    {
      trace.error() << " unknown output type. Try -h option. ";
      trace.info() << std::endl;
      return 0;
    }

  return 1;
}
