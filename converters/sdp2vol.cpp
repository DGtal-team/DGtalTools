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
 * @file sdp2vol.cpp
 * @ingroup converters
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2013/07/21
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
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/writers/GenericWriter.h"

#include "CLI11.hpp"

using namespace std;
using namespace DGtal;

/**
 @page sdp2vol sdp2vol
 @brief  Converts digital set of points into a volumic file.

@b Usage: sdp2vol [input] [output]

@b Allowed @b options @b are:

@code

Positionals:
  1 TEXT:FILE REQUIRED                  Sequence of 3d Discrete points (.sdp).
  2 TEXT                                Vol file  (.vol, .longvol, .pgm3d) 

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         Sequence of 3d Discrete points (.sdp).
  -o,--output TEXT                      Vol file  (.vol, .longvol, .pgm3d) 
  -f,--foregroundVal INT                value which will represent the foreground object in the resulting image (default 128)
  -b,--backgroundVal INT                value which will represent the background outside the  object in the resulting image (default 0)
  -d,--domain INT x 6                   customizes the domain of the resulting image xmin ymin zmin xmax ymax zmax (computed automatically by default)
  --invertY                             Invert the Y axis (image flip in the y direction)
@endcode

@b Example:
@code
  $ sdp2vol -i volumePoints.sdp -o volume.vol -d 0 0 0 10 10 10
@endcode

@see sdp2vol.cpp

*/





int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;

// parse command line using CLI ----------------------------------------------
   CLI::App app;
   std::string inputSDP;
   std::string outputFileName {"result.vol"};
   DGtal::int64_t rescaleInputMin {0};
   DGtal::int64_t rescaleInputMax {255};

   int foregroundVal {128};
   int backgroundVal {0};
   bool invertY {false};
   std::vector<int> domainCoords;

   app.description("Convert digital set of points into a volumic file.\n Example:\n sdp2vol -i volumePoints.sdp -o volume.vol -d 0 0 0 10 10 10 \n");
   app.add_option("-i,--input,1", inputSDP, "Sequence of 3d Discrete points (.sdp)." )
    ->required()
    ->check(CLI::ExistingFile);
   app.add_option("-o,--output,2", outputFileName, "Vol file  (.vol, .longvol, .pgm3d) ");
   app.add_option("-f,--foregroundVal", foregroundVal, "value which will represent the foreground object in the resulting image (default 128)");
   app.add_option("-b,--backgroundVal", backgroundVal, "value which will represent the background outside the  object in the resulting image (default 0)");
   app.add_option("-d,--domain", domainCoords, "customizes the domain of the resulting image xmin ymin zmin xmax ymax zmax (computed automatically by default)")
     ->expected(6);
   
   app.add_flag("--invertY", invertY,  "Invert the Y axis (image flip in the y direction)");
   

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  vector<unsigned int> vPos;
  vPos.push_back(0);
  vPos.push_back(1);
  vPos.push_back(2);
  trace.info() << "Reading input SDP file: " << inputSDP  ; 
  std::vector<Z3i::Point> vectPoints=  PointListReader<Z3i::Point>::getPointsFromFile(inputSDP, vPos); 
  trace.info() << " [done] " << std::endl ; 

  Z3i::Point ptLower;
  Z3i::Point ptUpper;
 
  struct BBCompPoints
  {
    explicit BBCompPoints(unsigned int d): myDim(d){};
    bool operator() (const Z3i::Point &p1, const Z3i::Point &p2){return p1[myDim]<p2[myDim];};
    unsigned int myDim;
  };
  if(domainCoords.size() != 6 )
  {
    unsigned int marge = 1;
    for(unsigned int i=0; i< 4; i++)
    {
      BBCompPoints cmp_points(i);
      ptUpper[i] = (*(std::max_element(vectPoints.begin(), vectPoints.end(), cmp_points)))[i]+marge;
      ptLower[i] = (*(std::min_element(vectPoints.begin(),  vectPoints.end(), cmp_points)))[i]-marge;
    }
  }
  else
  {
    ptLower = Z3i::Point(domainCoords[0],domainCoords[1], domainCoords[2]);
    ptUpper = Z3i::Point(domainCoords[3],domainCoords[4], domainCoords[5]);
  }
  
  Image3D::Domain imageDomain(ptLower, ptUpper);
  trace.info() << "domain: "<<imageDomain<<std::endl;
  Image3D imageResult(imageDomain); 
  for(Image3D::Domain::ConstIterator iter = imageResult.domain().begin();
      iter!= imageResult.domain().end();
      iter++)
  {
    imageResult.setValue(*iter, backgroundVal);
  }    
  
  for(unsigned int i=0; i<vectPoints.size(); i++)
  {
    if(invertY)
    {
      vectPoints[i][1]=ptUpper[1]-vectPoints[i][1];
    }
    if(imageResult.domain().isInside(vectPoints[i]))
    {
      imageResult.setValue(vectPoints[i], foregroundVal);
    }
    else
    {
      trace.warning() << "point " << vectPoints[i] << " outside the domain (ignored in the resulting volumic image)" << std::endl;  
    }
  }
  trace.info()<< "Exporting resulting volumic image ... ";
  GenericWriter<Image3D>::exportFile(outputFileName, imageResult);
  trace.info() << " [done]"<<std::endl;
  return EXIT_SUCCESS;  
}
