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
#include "DGtal/geometry/curves/GridCurve.h"

//Estimators
#include "DGtal/geometry/curves/BinomialConvolver.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantMeanCurvatureEstimator.h"

// Generation of resulting image
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/PPMWriter.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <vector>
#include <string>
#include <iomanip>


using namespace DGtal;




template< typename RealPoint >
struct OptionsIntegralInvariant
{
  double alpha; // <! Alpha parameter for the convolution kernel. 1/3 by default
  double radius; // <! Radius of the convolution kernel.
  RealPoint center; // <! Center of the shape.
  bool lambda_optimized;
};

/**
 *
 * @return Euclidean radius for the convolver of Integral Invariant estimators
 */
template< typename ConstIteratorOnPoints, typename Point >
unsigned int suggestedSizeIntegralInvariant( const double h,
                                             const Point& center,
                                             const ConstIteratorOnPoints& itb,
                                             const ConstIteratorOnPoints& ite )
{
  typedef typename Point::Component TValue;

  ConstIteratorOnPoints it = itb;
  Point p( *it );
  Point distance = p - center;
  TValue minRadius = distance.norm();
  ++it;

  for ( ; it != ite; ++it )
  {
    p = *it;
    distance = p - center;
    if ( distance.norm() < minRadius )
    {
      minRadius = distance.norm();
    }
  }

  return minRadius * h;
}



///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("FreemanChain,f", po::value<std::string>(), "FreemanChain file name")
    ("gridStepInit", po::value<double>()->default_value(1.0), "Grid step initial")
    ("gridStepIncrement", po::value<double>()->default_value(1.0), "Grid step increment ")
    ("output,o ", po::value<std::string>(), "set the output name ")
    ("gridStepFinal", po::value<double>()->default_value(1.0), "Grid step final")
    ("curvatureCutOff,c", po::value<double>()->default_value(10.0), "set the curvature limits to better display");
  

  typedef ImageContainerBySTLVector<Z2i::Domain, double > Image2D;

 
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }
  po::notify(vm);    
  if(!parseOK || vm.count("help")||argc<=1 || (!(vm.count("FreemanChain"))) )
    {
      trace.info()<< "Generate the Curvature Scale Sapce image using a binomial convolver based estimator." <<std::endl
                  << "The x axis is associated to the contour point and the y axis to the scale. The color represent the curvature values included between the cutoff values (set to 10 by default)."
                  <<std::endl << "Basic usage: "<<std::endl
      << "\t curvatureScaleSpace -f ${DGtal}/examples/samples/contourS.fc --gridStepInit 0.001 --gridStepIncrement  0.0005 --gridStepFinal 0.05 -o cssResu.ppm "<<std::endl
      << general_opt << "\n";
      return 0;
    }
  double h_initial = vm["gridStepInit"].as<double>();
  double h_increment = vm["gridStepIncrement"].as<double>();
  double h_final = vm["gridStepFinal"].as<double>();
  double curvatureCutOff = vm["curvatureCutOff"].as<double>();
  
  
  if(vm.count("FreemanChain")){
    std::string fileName = vm["FreemanChain"].as<std::string>();

    typedef Z2i::Space Space; 
    typedef Space::Point Point; 
    typedef Space::Integer Integer;  
    typedef FreemanChain<Integer> FreemanChain; 
    typedef std::vector< Point > Storage;
    typedef Storage::const_iterator ConstIteratorOnPoints; 
    typedef HueShadeColorMap<unsigned char> Hue;

    std::vector< FreemanChain > vectFcs =  PointListReader< Point >:: getFreemanChainsFromFile<Integer> (fileName); 
    
    bool isClosed = vectFcs.at(0).isClosed(); 
    std::cout << "# grid curve " << 1 << "/" << vectFcs.size() << " "
              << ( (isClosed)?"closed":"open" ) << std::endl;
    
    Storage vectPts; 
    FreemanChain::getContourPoints( vectFcs.at(0), vectPts ); 
    
    // Preparing resulting image:
    unsigned int height =  (int)((h_final-h_initial)/h_increment);
    Z2i::Domain domain (Z2i::Point(0,0), Z2i::Point(vectPts.size()-1, height-1));
    Image2D cssImage(domain);    
    HueShadeColorMap<double>  gradCurvature (-curvatureCutOff, curvatureCutOff);
    
    trace.progressBar(0, height);
    unsigned int nbLine=0;
    for(double h= h_initial; h <= h_final; h=h+h_increment){
      // Binomial estimator
      trace.progressBar(nbLine, height);
      
      typedef BinomialConvolver<ConstIteratorOnPoints, double> MyBinomialConvolver;
      typedef CurvatureFromBinomialConvolverFunctor< MyBinomialConvolver, double >   CurvatureBCFct;
      BinomialConvolverEstimator< MyBinomialConvolver, CurvatureBCFct> BCCurvatureEstimator;
      
      BCCurvatureEstimator.init( h, vectPts.begin(), vectPts.end(), isClosed );
      std::vector <double> curvatures( vectPts.size() ); 
      BCCurvatureEstimator.eval( vectPts.begin(), vectPts.end(), curvatures.begin() ); 
      
      // Output
      unsigned int j = 0;
      for ( ConstIteratorOnPoints it = vectPts.begin(), it_end = vectPts.end();
            it != it_end; ++it, ++j ) {
        double c = curvatures[j];
        c = c<-curvatureCutOff? -curvatureCutOff: c;
        c = c>curvatureCutOff? curvatureCutOff: c;
        cssImage.setValue(Z2i::Point(j,nbLine), c); 
      }
      nbLine++;
    }
    
    trace.progressBar(height, height);
    trace.info() <<std::endl;
    
    DGtal::GenericWriter<Image2D, 2, double, HueShadeColorMap<double> >::exportFile(vm["output"].as<std::string>(), cssImage, gradCurvature );       
  }
 
  return 0;
}

