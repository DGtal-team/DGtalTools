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
 * @file volImageMetrics.cpp
 *
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2013/07/20
 *
 *
 * This file is part of the DGtal library.
 */


/**
 @page volImageMetrics volImageMetrics
 
 @brief  Applies basic image measures (RMSE, PSNR) between two volumetric images A and B.

 @b Usage:  volImageMetrics --volA \<volAFilename\> --volB \<volBFilename\> 


 @b Allowed @b options @b are : 
 @code
  -h [ --help ]         display this message.
  -a [ --volA ] arg     Input filename of volume A (vol format, and other pgm3d
                        can also be used).
  -b [ --volB ] arg     Input filename of volume B (vol format, and other pgm3d
                        can also be used).
 @endcode

 @b Example: 

 @code
 # generating another input vol file using tutorial example (eroded.vol):
 $DGtal/build/examples/tutorial-examples/FMMErosion
 # compare the two images:
 $ volImageMetrics -a eroded.vol -b $DGtal/examples/samples/cat10.vol 
 @endcode


 You should obtain such an output:
@verbatim
# Image based measures (generated with volImageMetrics) given with the image A: eroded.voland the image B: /Users/kerautre/EnCours/DGtal/examples/samples/cat10.vol
#  RMSE PSNR 
 33.9411 171.331
@endverbatim
 
 @see  \ref volImageMetrics.cpp

 */


#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/math/Statistic.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <limits>


using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;

typedef ImageContainerBySTLVector < Z3i::Domain,  int > Image3D;
typedef ImageContainerBySTLVector < Z2i::Domain,  int > Image2D;


double
getRMSE(const Image3D & imageA, const Image3D &imageB){
  double sumDiff=0;
  for(Image3D::Domain::ConstIterator it = imageA.domain().begin(); it!=imageA.domain().end(); it++){
    sumDiff+=(imageA(*it)-imageB(*it))*(imageA(*it)-imageB(*it));
  }
  return sqrt(sumDiff/imageA.domain().size());
}


double
getPSNR(const Image3D & imageA, const Image3D &imageB, double rmsd){
  unsigned long long int d =  std::numeric_limits<Image3D::Value>::max();
  return 10.0*log10(d*d/rmsd);
}




int main(int argc, char**argv)
{
  
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "volA,a", po::value<std::string>(), "Input filename of volume A (vol format, and other pgm3d can also be used)." )
    ( "volB,b", po::value<std::string>(), "Input filename of volume B (vol format, and other pgm3d can also be used)." );
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);    
  
  if ( !parseOK || vm.count ( "help" ) || ! vm.count("volA")||! vm.count("volB")  )
    {
      trace.info() << "Apply basic image measures (RMSE, PSNR) between two volumetric images A and B."<<std::endl
		   << std::endl << "Basic usage: "<<std::endl
		   << "\t volImageMetrics --volA <volAFilename> --volB <volBFilename> "<<std::endl
		   << general_opt << "\n"
		   << "Typical use :\n  volImageMetrics -a imageA.vol  -b imageB.vol \n" ;

      return 0;
    }

  if(! vm.count("volA")||! vm.count("volB"))
    {
      trace.error() << " The two volume filename are needed to be defined" << endl;      
      return 0;
    }
 
  std::string volAFilename = vm["volA"].as<std::string>();
  std::string volBFilename = vm["volB"].as<std::string>();
  
  Image3D imageA = GenericReader<Image3D>::import(volAFilename);
  Image3D imageB = GenericReader<Image3D>::import(volBFilename);
 
  
  std::cout << "# Image based measures (generated with volImageMetrics) given with the image A: "<< volAFilename<< " and the image B: "<< volBFilename << endl;
  std::cout << "#  RMSE PSNR "<< endl;    
  
  double rmse= getRMSE(imageA, imageB);
  double psnr= getPSNR(imageA, imageB, rmse);
    
  std::cout << " " << rmse << " " << psnr << endl;
  
  return 1;
}

