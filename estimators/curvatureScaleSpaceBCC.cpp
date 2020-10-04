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
 * @file curvatureScaleSpace.cpp
 * @ingroup estimators
 *
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 *
 * @date 2014/04/02
 *
 * Output the image of curvature scale space
 * 
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "DGtal/base/Common.h"

#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"

//Grid curve
#include "DGtal/geometry/curves/FreemanChain.h"

//Estimators
#include "DGtal/geometry/curves/BinomialConvolver.h"


// Generation of resulting image
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/PPMWriter.h"

#include "CLI11.hpp"

#include <vector>
#include <string>
#include <iomanip>


using namespace DGtal;



/**
 @page curvatureScaleSpaceBCC curvatureScaleSpaceBCC
 
 @brief Generate the Curvature Scale Space image using a binomial convolver based estimator.

The x axis is associated to the contour point and the y axis to the scale. The colors represent the curvature values included between the cutoff values (set to 10 by default).

 @b Usage:  curvatureScaleSpaceBCC --input <filename>  --output <filename> 


 @b Allowed @b options @b are : 
 @code
  Positionals:
  1 TEXT:FILE REQUIRED                  Input FreemanChain file name
  2 TEXT REQUIRED                       Set the output name

  Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         Input FreemanChain file name
  -o,--output TEXT REQUIRED             Set the output name
  --gridStepInit FLOAT=1                Grid step initial
  --gridStepIncrement FLOAT=1           Grid step increment 
  --gridStepFinal FLOAT=1               Grid step final
  -c,--curvatureCutOff FLOAT=10         set the curvature limits to better display
 @endcode

 @b Example: 
 @code
$ curvatureScaleSpaceBCC -i ${DGtal}/examples/samples/contourS.fc --gridStepInit 0.001 --gridStepIncrement  0.0005 --gridStepFinal 0.1 -o cssResu.ppm
 @endcode


 You should obtain such a result:
 @image html resCurvatureScaleSpaceBCC.png "Resulting visualization."
 
 @see
 @ref curvatureScaleSpaceBCC.cpp

 */




void
computeCurvatureBCC(double h, const FreemanChain<Z2i::Integer> &fc, std::vector<double> &resCurvature, 
                    bool isClosed){
  std::vector<Z2i::Point> vectPoints;
  FreemanChain<Z2i::Integer>::getContourPoints( fc, vectPoints ); 
  
  typedef BinomialConvolver<std::vector<Z2i::Point>::const_iterator, double> MyBinomialConvolver;
  typedef CurvatureFromBinomialConvolverFunctor< MyBinomialConvolver, double >   CurvatureBCFct;
  BinomialConvolverEstimator< MyBinomialConvolver, CurvatureBCFct> BCCurvatureEstimator;
  BCCurvatureEstimator.init( h, vectPoints.begin(), vectPoints.end(), isClosed );
  resCurvature.clear();
  resCurvature.resize(vectPoints.size());
  BCCurvatureEstimator.eval( vectPoints.begin(), vectPoints.end(), resCurvature.begin() ); 
}



///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
 
  // parse command line CLI ----------------------------------------------
  CLI::App app;
  std::string fileName;
  std::string outputFileName;
  double h_initial {1.0};
  double h_increment {1.0};
  double h_final {1.0};
  double curvatureCutOff {10};

  app.description("Generate the Curvature Scale Space image using a binomial convolver based estimator.\n Typical use example:\n \t curvatureScaleSpaceBCC --input <filename>  --output <filename>\n");
  app.add_option("-i,--input,1",fileName,"Input FreemanChain file name")->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output,2",outputFileName,"Set the output name")->required();  
  app.add_option("--gridStepInit", h_initial, "Grid step initial", true);
  app.add_option("--gridStepIncrement", h_increment, "Grid step increment ",true);
  app.add_option("--gridStepFinal", h_final, "Grid step final", true);
  app.add_option("--curvatureCutOff,-c", curvatureCutOff, "set the curvature limits to better display", true);

  app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------------  

  typedef ImageContainerBySTLVector<Z2i::Domain, double > Image2D;

    
  std::vector< DGtal::FreemanChain<Z2i::Integer>  > vectFcs =  PointListReader< Z2i::Point >:: getFreemanChainsFromFile<Z2i::Integer> (fileName);     
  bool isClosed = vectFcs.at(0).isClosed(); 
  
  
  // Preparing resulting image:
  unsigned int height =  (int)((h_final-h_initial)/h_increment);
  // We add one point if the freemnchain is open.
  Z2i::Domain domain (Z2i::Point(0,0), Z2i::Point(vectFcs.at(0).size()+(isClosed? 0: 1), height-1));
  Image2D cssImage(domain);    
  HueShadeColorMap<double>  gradCurvature (-curvatureCutOff, curvatureCutOff);
  
  trace.progressBar(0, height);
  double h= h_initial;
  for(double l= 0; l < height; l++ ){
    // Binomial estimator
    trace.progressBar(l, height);
    std::vector<double> curvaturesBCC;
    computeCurvatureBCC(h, vectFcs.at(0), curvaturesBCC, isClosed);
    // Output
    unsigned int j = 0;
    for ( std::vector<double>::const_iterator it = curvaturesBCC.begin(), it_end = curvaturesBCC.end();
          it != it_end; ++it, ++j ) {
      double c = *it;
      c = c<-curvatureCutOff? -curvatureCutOff: c;
      c = c>curvatureCutOff? curvatureCutOff: c;
      cssImage.setValue(Z2i::Point(j, l), c); 
    }      
    h=h+h_increment;
  }
  trace.progressBar(height, height);
  trace.info() <<std::endl;
  DGtal::GenericWriter<Image2D, 2, DGtal::Color, HueShadeColorMap<double> >::exportFile(outputFileName, cssImage, gradCurvature );
  return 0;
}

