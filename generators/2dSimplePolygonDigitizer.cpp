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
 * @file 2dSimplePolygonDigitizer.cpp
 * @ingroup Generators
 * @author Phuc Ngo (\c hoai-diem-phuc.ngo@loria.fr)
 * LORIA - Lorraine Univeristy , France
 *
 * @date 2021/04/06
 *
 * Gauss Digitization of a simple polyline.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <string>
#include <iterator>
#include <fstream>

#include "CLI11.hpp"

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

//image
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/boards/Board2D.h"

//contour
#include "DGtal/io/readers/PointListReader.h"

using namespace DGtal;


/**
 @page 2dSimplePolygonDigitizer 2dSimplePolygonDigitizer
 @brief Compute the Gauss Digitization of a simple closed polyline (no hole or self-intersection).
 @ingroup generatorstools

 
 The digitizer compute the set of integer points inside the input polyline.
 
 @b Usage: 2dSimplePolygonDigitizer [input] [output]
 
 @b Allowed @b options @b are:
 
 @code
 
 Positionals:
 1 TEXT:FILE REQUIRED                  Input polyline filename (sdp).
 
 Options:
 
 Positionals:
 1 TEXT:FILE REQUIRED                  Input polyline file name.
 
 Options:
 -h,--help                             Print this help message and exit
 -i,--input TEXT:FILE REQUIRED         Input sdp filename.
 -o,--output TEXT=result.pgm           the output image filename (pgm or svg)
 @endcode
 
 @b Example:
 @code
 $2dSimplePolygonDigitizer -i ${DGtal}/tests/samples/flower-30-8-3.sdp -o sample.pgm
 @endcode
 You will obtain such image:
 @image html  resPolygonDigitizer.png "Resulting image" 
 
 @b Example @b with @b more @b complex  @b contours:
  
  The file located in $DGtal/examples/samples/contourS.sdp
 
  @code
  $ 2dSimplePolygonDigitizer -i  $DGtal/examples/samples/contourS.sdp -o sample2.pgm
  @endcode
  
  You will obtain such image:
  @image html  resPolygonDigitizer2.png "Resulting image"
 
 @see @ref 2dSimplePolygonDigitizer
 @ref 2dSimplePolygonDigitizer.cpp
 
 */

std::pair<Z2i::Point, Z2i::Point > getBoundingBox(const std::vector<Z2i::Point>& vecPts);
std::vector<Z2i::Point> GaussDigization(const std::vector<Z2i::Point>& polygon);

int main( int argc, char** argv )
{
    
    // parse command line using CLI ----------------------------------------------
    CLI::App app;
    std::string inputFileName;
    std::string outputFileName="result.pgm";
    
    app.description("compute the set of integer points inside the input polyline.\n Basic example:\n \t 2dSimplePolygonDigitizer  -i  inputPolyline.sdp -o contourDisplay.pgm");
    app.add_option("-i,--input,1", inputFileName, "Input filename (freeman chain of sdp)." )
    ->required()
    ->check(CLI::ExistingFile);
    app.add_option("-o,--output,2", outputFileName, "Output filename");
    
    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------
    
    std::string inputExt = inputFileName.substr(inputFileName.find_last_of(".")+1);
    if(inputExt != "sdp"){
        trace.error() << "input file should be sdp file" << std::endl;
        return EXIT_FAILURE;
    }
    
    std::vector<Z2i::Point> poly = PointListReader<Z2i::Point>::getPointsFromFile(inputFileName);
    std::vector<Z2i::Point> gd = GaussDigization(poly);
    
    std::string outputExt = outputFileName.substr(outputFileName.find_last_of(".")+1);
    std::string outputBaseName = outputFileName.substr(0, outputFileName.find_last_of("."));
    
    if(outputExt == "pgm") {
        //Write results to pgm image
        std::pair<Z2i::Point, Z2i::Point> bb = getBoundingBox(gd);
        typedef ImageSelector <Z2i::Domain, unsigned char>::Type Image;
        Z2i::Point p1 =  bb.first-Z2i::Point(1,1);
        Z2i::Point p2 =  bb.second+Z2i::Point(1,1);
        Image image(Z2i::Domain(p1,p2));
        for(std::vector<Z2i::Point>::const_iterator it=gd.begin(); it!=gd.end(); it++) {
            Z2i::Point p ((*it)[0],(*it)[1]);
            image.setValue(p,255);
        }
        PGMWriter<Image>::exportPGM(outputFileName,image);
    }
    else {
        //Write results to svg image
        outputFileName = outputBaseName + ".svg";
        DGtal::Board2D board;
        board << DGtal::SetMode("PointVector", "Both");
        for(std::vector<Z2i::Point>::const_iterator it=gd.begin(); it!=gd.end(); it++) {
            board << *it;
        }
        board.setPenColor(DGtal::Color::Red);
        board.setLineWidth(2);
        for(size_t it=0; it<poly.size()-1; it++) {
            DGtal::Z2i::RealPoint p1 (poly.at(it)[0],poly.at(it)[1]);
            DGtal::Z2i::RealPoint p2 (poly.at(it+1)[0],poly.at(it+1)[1]);
            board.drawLine(p1[0],p1[1],p2[0],p2[1]);
        }
        DGtal::Z2i::RealPoint p1 (poly.front()[0],poly.front()[1]);
        DGtal::Z2i::RealPoint p2 (poly.back()[0],poly.back()[1]);
        board.drawLine(p1[0],p1[1],p2[0],p2[1]);

        board.saveSVG(outputFileName.c_str());
    }
    return EXIT_SUCCESS;
}

// Auxiliaires functions
// Compute the bounding box of the polyline
std::pair<Z2i::Point, Z2i::Point > getBoundingBox(const std::vector<Z2i::Point>& vecPts) {
    int minX = vecPts.at(0)[0];
    int maxX = vecPts.at(0)[0];
    int minY = vecPts.at(0)[1];
    int maxY = vecPts.at(0)[1];
    
    for(size_t it=1; it<vecPts.size(); it++) {
        int x = vecPts.at(it)[0];
        int y = vecPts.at(it)[1];
        if(minX > x) minX = x;
        if(maxX < x) maxX = x;
        if(minY > y) minY = y;
        if(maxY < y) maxY = y;
    }
    Z2i::Point tl (minX,minY);
    Z2i::Point rb (maxX,maxY);
    return std::make_pair(tl,rb);
}

/* Source adapted from : https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/ */
// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
template<typename T>
T max(T v1, T v2) {
    return v1 > v2 ? v1 : v2;
}
template<typename T>
T min(T v1, T v2) {
    return v1 > v2 ? v2 : v1;
}
template<typename TPoint>
bool onSegment(TPoint p, TPoint q, TPoint r)
{
    if (q[0] <= max(p[0], r[0]) && q[0] >= min(p[0], r[0]) && q[1] <= max(p[1], r[1]) && q[1] >= min(p[1], r[1]))
        return true;
    return false;
}
 
// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
template<typename TPoint>
int orientation(TPoint p, TPoint q, TPoint r)
{
    int val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1]);
    if (val == 0) return 0; // colinear
    return (val > 0)? 1: 2; // clock or counterclock wise
}
 
// The function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
template<typename TPoint>
bool doIntersect(TPoint p1, TPoint q1, TPoint p2, TPoint q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);
 
    // General case
    if (o1 != o2 && o3 != o4)
        return true;
 
    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;
 
    // p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;
 
    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;
 
    // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;
 
    return false; // Doesn't fall in any of the above cases
}

//Find direction that does not pass through any vertex of the polygon
template<typename TPoint>
TPoint findDirection(const std::vector<TPoint>& polygon, TPoint p, TPoint d1, TPoint d2 )
{
    TPoint d = d1;
    bool belong = true;
    //Find the ray that doest pass by any vertices of polygon
    do {
        belong = false;
        for(size_t it=0; it<polygon.size(); it++) {
            // Check p is a vertex of polygon, then it is inside
            if(polygon.at(it)==p)
                return p;
            else {
                // Check if 'polygon[i]' belongs to the line segment from 'p' to 'extreme',
                // if it is then launch another ray using dichotomy search
                if(orientation(p, polygon.at(it), d)==0) {
                    d[0] = (d[0] + d2[0])/2;
                    d[1] = (d[1] + d2[1])/2;
                    belong = true;
                }
            }
        }
    } while (belong);
    return d;
}

// Returns true if the point p lies inside the polygon[] with n vertices
template<typename TPoint>
bool isInsidePolygon(const std::vector<TPoint>& polygon, TPoint p)
{
    int n = polygon.size();
    // There must be at least 3 vertices in polygon[]
    if (n < 3) return false;
 
    // Create a point for line segment from p to infinite
    // Find the direction without containing any polygon vertices
    TPoint extreme1 (p[0], INT_MAX/1e6);
    TPoint extreme2 (INT_MAX/1e6, p[1]);
    TPoint extreme = findDirection(polygon, p, extreme1, extreme2);

    // Check p is a vertex of polygon, then it is inside
    if(extreme==p)
        return true;
    // Count intersections of the above line with sides of polygon
    int count = 0, i = 0;
    do
    {
        int next = (i+1)%n;
 
        // Check if the line segment from 'p' to 'extreme' intersects
        // with the line segment from 'polygon[i]' to 'polygon[next]'
        if (doIntersect(polygon[i], polygon[next], p, extreme))
        {
            // If the point 'p' is colinear with line segment 'i-next',
            // then check if it lies on segment. If it lies, return true,
            // otherwise false
            if (orientation(polygon[i], p, polygon[next]) == 0)
            return onSegment(polygon[i], p, polygon[next]);
 
            count++;
        }
        i = next;
    } while (i != 0);
 
    // Return true if count is odd, false otherwise
    return count&1; // Same as (count%2 == 1)
}

// Gauss digitizer
std::vector<Z2i::Point> GaussDigization(const std::vector<Z2i::Point>& polygon)
{
    typedef Z2i::Point TPoint;
    std::vector<TPoint> gd;
    std::pair<TPoint, TPoint> bb = getBoundingBox(polygon);
    int minx = int(bb.first[0])-1;
    int miny = int(bb.first[1])-1;
    int maxx = int(bb.second[0])+1;
    int maxy = int(bb.second[1])+1;
    for(int x = minx; x<maxx; x++) {
        for(int y = miny; y<maxy; y++) {
            TPoint p(x,y);
            bool ok = isInsidePolygon<TPoint>(polygon,p);
            if(ok) {
                TPoint pi(p[0],p[1]);
                gd.push_back(pi);
            }
        }
    }
    return gd;
}
