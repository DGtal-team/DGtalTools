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
 * @file curvatureBC.cpp
 * @ingroup estimators
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr ) 
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS,
 * France
 *
 * @date 2011/07/13
 *
 * Output the curvature of the Freeman code of a grid curve
 * using the Binomial convolver
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


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <vector>
#include <string>
#include <iomanip>

using namespace DGtal;



/**
 @page curvatureBC curvatureBC
 
 @brief Estimatates curvature using a binomial convolver.

 @b Usage: curvatureBC [options] --input  \<fileName\> 


 @b Allowed @b options @b are : 
 @code
  -h [ --help ]              display this message
  -i [ --input ] arg         input FreemanChain file name
  -s [ --GridStep ] arg (=1) Grid step
 @endcode

 @b Example: 
 
We consider as input shape the freeman chain of the DGtal/examples/sample directory. The contour can be displayed with @ref displayContours : 
 @code
$  displayContours -i $DGtal/examples/samples/contourS.fc -o contourS.png  --drawPointOfIndex 0
 @endcode

The curvature can be computed as follows:
@code 
$ curvatureBC -i $DGtal/examples/samples/contourS.fc  > curvatureBC.dat
$ gnuplot 
gnuplot> plot 'curvatureBC.dat' w lines title "curvature with Binomial Convolution estimator"
@endcode



 You should obtain such a result:

 | contour  | curvature  | 
 | :------: | :--------: |   
 | ![ ](resCurvatureBCcontour.png)  | ![ ](resCurvatureBCcurvature.png)  |
 | CCW oriented (index 0=blue pt)| resulting curvature |

 @see
 @ref curvatureBC.cpp

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
    ("GridStep,step", po::value<double>()->default_value(1.0), "Grid step");
  
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }
  po::notify(vm);    
  if(!parseOK || vm.count("help")||argc<=1 || (!(vm.count("input"))) )
    {
      trace.info()<< "Curvature using a binomial convolver " <<std::endl << "Basic usage: "<<std::endl
      << "\t curvatureBC [options] --input  <fileName> "<<std::endl
      << general_opt << "\n";
      return 0;
    }
  
  
  double h = vm["GridStep"].as<double>();  


 
  if(vm.count("input")){
    std::string fileName = vm["input"].as<std::string>();

    typedef Z2i::Space Space; 
    typedef Space::Point Point; 
    typedef Space::Integer Integer;  
    typedef FreemanChain<Integer> FreemanChain; 
    typedef std::vector< Point > Storage;
    typedef Storage::const_iterator ConstIteratorOnPoints; 

    std::vector< FreemanChain > vectFcs =  PointListReader< Point >:: getFreemanChainsFromFile<Integer> (fileName); 
   

    for(unsigned int i=0; i<vectFcs.size(); i++){

      bool isClosed = vectFcs.at(i).isClosed(); 
      std::cout << "# grid curve " << i+1 << "/" << vectFcs.size() << " "
      << ( (isClosed)?"closed":"open" ) << std::endl;

      Storage vectPts; 
      FreemanChain::getContourPoints( vectFcs.at(i), vectPts ); 

      // Binomial
      std::cout << "# Curvature estimation from binomial convolution" << std::endl;
      typedef BinomialConvolver<ConstIteratorOnPoints, double> MyBinomialConvolver;
      std::cout << "# mask size = " << 
      MyBinomialConvolver::suggestedSize( h, vectPts.begin(), vectPts.end() ) << std::endl;
      typedef CurvatureFromBinomialConvolverFunctor< MyBinomialConvolver, double >
      CurvatureBCFct;
      BinomialConvolverEstimator< MyBinomialConvolver, CurvatureBCFct> BCCurvatureEstimator;

      BCCurvatureEstimator.init( h, vectPts.begin(), vectPts.end(), isClosed );

      std::vector <double> curvatures( vectPts.size() ); 
      BCCurvatureEstimator.eval( vectPts.begin(), vectPts.end(), curvatures.begin() ); 

      // Output
      std::cout << "# id curvature" << std::endl;  
      unsigned int j = 0;
      for ( ConstIteratorOnPoints it = vectPts.begin(), it_end = vectPts.end();
      it != it_end; ++it, ++j ) {
  std::cout << j << std::setprecision( 15 )
       << " " << curvatures[ j ] << std::endl;
      }

   }

 }
 
  return 0;
}

