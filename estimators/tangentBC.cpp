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
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2011/07/14
 *
 * Output the tangent of the Freeman code of a grid curve
 * using the Binomial convolver
 * 
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/PointListReader.h"

//Grid curve
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/geometry/curves/GridCurve.h"

//Estimators
#include "DGtal/geometry/curves/BinomialConvolver.h"

using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;



/**
 @page tangentBC tangentBC
 
 @brief  Estimates tangent using a binomial convolver.

 @b Usage:  tangentBC [options] --input  <fileName> 


 @b Allowed @b options @b are : 
 @code
  -h [ --help ]              display this message
  -i [ --input ] arg         input file name: FreemanChain (.fc) or a sequence of discrete points (.sdp).
  -s [ --GridStep ] arg (=1) Grid step
 @endcode

@note The file may contain several freeman chains.


 @b Example: 

We consider as input shape the freeman chain of the DGtal/examples/sample directory. The contour can be displayed with @ref displayContours : 
 @code
$  displayContours -i $DGtal/examples/samples/contourS.fc -o contourS.png  --drawPointOfIndex 0
 @endcode

The tangents can be computed as follows:
@code 
$ tangentBC -i $DGtal/examples/samples/contourS.fc  > tangentsBC.dat
$ gnuplot 
gnuplot> plot [] [-1.2:1.2]'tangentsBC.dat' using 1:3  w lines title "tangents with Binomial Convolution estimator"
@endcode

 You should obtain such a result:

 | contour  | curvature  | 
 | :------: | :--------: |   
 | ![ ](resCurvatureBCcontour.png)  | ![ ](resTangentBC.png)  |
 | CCW oriented (index 0=blue pt)| resulting tangent (angle) |
 
 @see
 @ref tangentBC.cpp

 */


int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "input file name: FreemanChain (.fc) or a sequence of discrete points (.sdp).")
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
      trace.info()<< "Tangent using a binomial convolver " <<std::endl << "Basic usage: "<<std::endl
      << "\t tangentBC [options] --input  <fileName> "<<std::endl
      << general_opt << "\n"
      << "NB: the file may contain several freeman chains." << "\n";
      return 0;
    }
  
  
  double h = vm["GridStep"].as<double>();  


 
  if(vm.count("input")){
    std::string fileName = vm["input"].as<std::string>();
    std::string extension =  fileName.substr( fileName.find_last_of(".") + 1 );
    bool isSDP = extension == "sdp";
    typedef Z2i::Space Space; 
    typedef Space::Point Point; 
    typedef PointVector<2, double> RealPoint; 
    typedef Space::Integer Integer;  
    typedef FreemanChain<Integer> FreemanChain; 
    typedef std::vector< Point > Storage;
    typedef Storage::const_iterator ConstIteratorOnPoints; 

    std::vector< FreemanChain > vectFcs;
    if(!isSDP)
      {
        vectFcs =   
          PointListReader< Point >:: getFreemanChainsFromFile<Integer> (fileName); 
      }
    for(unsigned int i=0; i<vectFcs.size() || (i==0 && isSDP); i++){
      Storage vectPts; 
      bool isClosed;
      if(!isSDP)
        {
          isClosed = vectFcs.at(i).isClosed(); 
          std::cout << "# grid curve " << i << "/" << vectFcs.size() << " "
                    << ( (isClosed)?"closed":"open" ) << std::endl;
          FreemanChain::getContourPoints( vectFcs.at(i), vectPts ); 
        }
      else
        {
          vectPts = PointListReader<Z2i::Point>::getPointsFromFile(fileName);
          Z2i::Point pf =vectPts[0];
          Z2i::Point pl =vectPts[vectPts.size()-1];
          isClosed = (pf[0]-pl[0])+(pf[1]-pl[1]) <= 1;
        }
      
      // Binomial
      std::cout << "# Curvature estimation from binomial convolution" << std::endl;
      typedef BinomialConvolver<ConstIteratorOnPoints, double> MyBinomialConvolver;
      std::cout << "# mask size = " << 
      MyBinomialConvolver::suggestedSize( h, vectPts.begin(), vectPts.end() ) << std::endl;
      typedef 
        TangentFromBinomialConvolverFunctor< MyBinomialConvolver, RealPoint >
        TangentBCFct;
      BinomialConvolverEstimator< MyBinomialConvolver, TangentBCFct> 
        BCTangentEstimator;
      
      BCTangentEstimator.init( h, vectPts.begin(), vectPts.end(), isClosed );

      std::vector<RealPoint> tangents( vectPts.size() ); 
      BCTangentEstimator.eval( vectPts.begin(), vectPts.end(), 
             tangents.begin() ); 

      // Output
      std::cout << "# id tangent.x tangent.y angle(atan2(y,x)) x y" << std::endl;  
      unsigned int j = 0;
      for ( ConstIteratorOnPoints 
        it = vectPts.begin(), it_end = vectPts.end();
      it != it_end; ++it, ++j ) 
  {
    double x = tangents[ j ][ 0 ];
    double y = tangents[ j ][ 1 ];
    std::cout << j << std::setprecision( 15 )
              << " " << x << " " << y 
              << " " << atan2( y, x ) << " " << (*it)[0] << " " << (*it)[1]
              << std::endl;
  }

    }

  }
  return 0;
}

