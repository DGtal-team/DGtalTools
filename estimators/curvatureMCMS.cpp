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
 * @file curvatureMCMS.cpp
 * @ingroup estimators
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr ) 
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS,
 * France
 *
 * @date 2011/07/19
 *
 * Output the curvature of the Freeman code of a grid curve
 * using the maximal segments cover
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
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h"
#include "DGtal/geometry/curves/estimation/SegmentComputerEstimators.h"


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <vector>
#include <string>
#include <iomanip>


using namespace DGtal;
using namespace std;


/**
 @page curvatureMCMS curvatureMCMS
 
 @brief Estimates curvature using length of most centered segment computers.

 @b Usage:  curvatureMCMS [options] --input  <fileName> 


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
$ curvatureMCMS -i $DGtal/examples/samples/contourS.fc  > curvatureMCMS.dat
$ gnuplot 
gnuplot> plot 'curvatureMCMS.dat' w lines title "curvature using length of most centered segment computers estimator"
@endcode



 You should obtain such a result:

 | contour  | curvature  | 
 | :------: | :--------: |   
 | ![ ](resCurvatureBCcontour.png)  | ![ ](resCurvatureMCMS.png)  |
 | CCW oriented (index 0=blue pt)| resulting curvature |
 
 @see
 @ref curvatureMCMS.cpp

 */



///////////////////////////////////////////////////////////////////////////////
/*
  Curvature estimation from the length of the most centered maximal segments
  @param h grid step
  @param itb begin iterator on points
  @param ite end iterator on points
  @param ito output iterator on curvature
 */
template<typename I, typename O>
void estimationFromLength( double h, const I& itb, const I& ite, const O& ito )
{
  typedef ArithmeticalDSSComputer<I,int,4> SegmentComputer;
  SegmentComputer sc;
  typedef CurvatureFromDSSLengthEstimator<SegmentComputer> SCEstimator;
  SCEstimator f; 
  typedef MostCenteredMaximalSegmentEstimator<SegmentComputer,SCEstimator> Estimator;
  Estimator CurvatureEstimator(sc, f);

  CurvatureEstimator.init( h, itb, ite );
  CurvatureEstimator.eval( itb, ite, ito ); 
}

/*
  Curvature estimation from the length and width
  of the most centered maximal segments
  @param h grid step
  @param itb begin iterator on points
  @param ite end iterator on points
  @param ito output iterator on curvature
 */
template<typename I, typename O>
void estimationFromLengthAndWidth( double h, const I& itb, const I& ite, const O& ito )
{
  typedef ArithmeticalDSSComputer<I,int,4> SegmentComputer;
  SegmentComputer sc;
  typedef CurvatureFromDSSEstimator<SegmentComputer> SCEstimator;
  SCEstimator f; 
  typedef MostCenteredMaximalSegmentEstimator<SegmentComputer,SCEstimator> Estimator;
  Estimator CurvatureEstimator(sc, f);

  CurvatureEstimator.init( h, itb, ite );
  CurvatureEstimator.eval( itb, ite, ito ); 
}


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
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
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }

  po::notify(vm);    
  if(!parseOK || vm.count("help")||argc<=1 || (!(vm.count("input"))) )
    {
      trace.info() << "Estimates curvature using length of most centered segment computers. " << std::endl; 
      trace.info() << "Basic usage: " << std::endl
		   << "\t curvatureMCMS [options] --input  <fileName> "<< std::endl
		   << general_opt << "\n";
      return 0;
    }
  
  
  double h = vm["GridStep"].as<double>();  
 
  if(vm.count("input")){
    string fileName = vm["input"].as<string>();

    typedef Z2i::Space Space; 
    typedef Space::Point Point; 
    typedef Space::Integer Integer;  
    typedef Z2i::KSpace KSpace; 
    typedef FreemanChain<Integer> FreemanChain; 

    vector< FreemanChain > vectFcs =  PointListReader< Point >::getFreemanChainsFromFile<Integer> (fileName);  

    for(unsigned int i=0; i<vectFcs.size(); i++){

      // Freeman chain
      FreemanChain fc = vectFcs.at(i); 
      // Create GridCurve
      GridCurve<> gridcurve;
      gridcurve.initFromPointsRange( fc.begin(), fc.end() );

      cout << "# grid curve " << i+1 << "/" 
	   << gridcurve.size() << " "
	   << ( (gridcurve.isClosed())?"closed":"open" ) << endl;

      // Create range of incident points
      typedef GridCurve<KSpace>::PointsRange Range;
      Range r = gridcurve.getPointsRange();//building range
 
      // Estimation
      cout << "# Curvature estimation from maximal segments" << endl; 
      std::vector<double> estimations1; 
      std::vector<double> estimations2; 
      if (gridcurve.isOpen())
        { 
	  cout << "# open grid curve" << endl;
	  estimationFromLength( h, r.begin(), r.end(), back_inserter(estimations1) ); 
	  estimationFromLengthAndWidth( h, r.begin(), r.end(), back_inserter(estimations2) ); 
        }
      else
        { 
	  cout << "# closed grid curve" << endl;
	  estimationFromLength( h, r.c(), r.c(), back_inserter(estimations1) ); 
	  estimationFromLengthAndWidth( h, r.c(), r.c(), back_inserter(estimations2) ); 
        }

      // Output
      cout << "# id curvatureFromLength curvatureFromLengthAndWidth" << endl;  
      unsigned int j = 0;
      for ( Range::ConstIterator it = r.begin(), itEnd = r.end();
      it != itEnd; ++it, ++j ) {
	cout << j << setprecision( 15 )
	     << " " << estimations1[ j ]
	     << " " << estimations2[ j ] << endl;
      }

   }

 }
 
  return 0;
}

