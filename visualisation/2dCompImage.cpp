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
 * @file
 * @ingroup visualisation
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Lorraine, France
 *
 * @date 2018/01/11
 *
 * Source file of the tool 2dCompImage
 *
 * This file is part of the DGtal library/DGtalTools Project.
 */

///////////////////////////////////////////////////////////////////////////////
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/GenericWriter.h>


#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/io/writers/PPMWriter.h>
#include <DGtal/math/Statistic.h>
#include "DGtal/io/colormaps/GradientColorMap.h"
#include <DGtal/math/Statistic.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace DGtal;
///////////////////////////////////////////////////////////////////////////////

typedef ImageContainerBySTLVector < Z2i::Domain,  unsigned char > Image2D;
typedef ImageContainerBySTLVector < Z2i::Domain,  unsigned int > Image2DErr;
typedef GradientColorMap<unsigned int, CMAP_JET, 1 > JetMap;

namespace po = boost::program_options;


/**
 @page Doc2dCompImage 2dCompImage
 
 @brief Compare images and displays differences (squared and absolute differences). 


 @b Usage:  2dCompImage --imageA <imageA>.pgm --imageB <imageB>.pgm --imageError <name> 

 @b Allowed @b options @b are :
 
 @code
  -h [ --help ]                    display this message
  -a [ --imageA ] arg              Input filename of image A.
  -b [ --imageB ] arg              Input filename of image B.
  -e [ --imageError ] arg          Output error image basename (will generate 
                                   two images <basename>MSE.ppm and 
                                   <basename>MAE.ppm).
  -S [ --setMaxColorValueMSE ] arg Set the maximal color value for the scale 
                                   display of MSE (else the scale is set the 
                                   maximal MSE value).
  -A [ --setMaxColorValueMAE ] arg Set the maximal color value for the scale 
                                   display of MAE (else the scale is set from 
                                   the maximal MAE value).
 @endcode

 @b Example: 
 Typical use example:
 @code
       2dCompImage -a imageA.pgm -b imageB.pgm -e errorImage -S 100 
 @endcode

You should obtain such a visualisation:
 @image html res2dCompImage.png "resulting visualisation of absolute error between two images."

 @see
 @ref 2dCompImage.cpp


 */




/**
 * Compute statistics on the mean absolute difference between pixel values of image A and B. 
 **/
Statistic< int>
getMAEstats(const Image2D & imageA, const Image2D &imageB, Image2DErr &imageMAE)
{
  Statistic< int> stat(false);
  for(auto const &point: imageA.domain())
  {
    unsigned int error = abs((imageA(point)-imageB(point)));
    stat.addValue(error);
    imageMAE.setValue(point, error);     
  }
  stat.terminate();
  return stat;
}


/**
 * Compute statistics on the mean squared difference between pixel values of image A and B. 
 **/
Statistic< int>
getMSEstats(const Image2D & imageA, const Image2D &imageB, Image2DErr &imageMSE)
{
  Statistic< int> stat(false);
  for(Image2D::Domain::ConstIterator it = imageA.domain().begin();
      it!=imageA.domain().end(); it++)
  {
    int error = (imageA(*it)-imageB(*it))*(imageA(*it)-imageB(*it));
    stat.addValue(error);
    imageMSE.setValue(*it, error);     
  }
  stat.terminate();
  return stat;
}



/**
 * Displays std statistics values
 **/
template<typename StatisticT>
void
displayStats(const StatisticT &aStat, const string &type)
{
  std::cout << "# Stats on " << type << ": mean max min unbiased_variance nb_samples  " << std::endl; 
  std::cout << aStat.mean() << " " << aStat.max() << " " << aStat.min() << " "
            << " " << aStat.unbiasedVariance() << " " << aStat.samples() << std::endl;  
}


int main( int argc, char** argv )
{
  // parse command line -------------------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("imageA,a", po::value<std::string >(), "Input filename of image A." )
    ("imageB,b", po::value<std::string >(), "Input filename of image B." )
    ("imageError,e", po::value<std::string >(), "Output error image basename (will generate two images <basename>MSE.ppm and <basename>MAE.ppm)." )
    ("setMaxColorValueMSE,S", po::value<int>(), "Set the maximal color value for the scale display of MSE (else the scale is set the maximal MSE value)." )
    ("setMaxColorValueMAE,A", po::value<int>(), "Set the maximal color value for the scale display of MAE (else the scale is set from the maximal MAE value).");


  bool parseOK=true;
  po::variables_map vm;
  try
    {
      po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }
  catch(const std::exception& ex)
    {
      parseOK=false;
      trace.info()<< "Error checking program options: "<< ex.what()<< endl;
    }
  

  // check if min arguments are given and tools description ------------------
  po::notify(vm);
  if( !parseOK || vm.count("help")|| argc<=1 || !vm.count("imageA")
      ||  !vm.count("imageB") )
    {
      std::cout << "Usage: " << argv[0] << " --imageA <imageA>.pgm --imageB <imageB>.pgm -imageError <name> \n"
                << "Compare images and displays differences (squared and absolute differences).  \n"
                << general_opt << "\n"
                << "Typical use example:\n \t 2dCompImage -a imageA.pgm -b imageB.pgm -e errorImage -S 100 \n";
      return 0;
    }  
 if(!vm.count("imageA") ||  !vm.count("imageB"))
    {
      trace.error() << " The two images filename are needed to be defined" << endl;      
      return 0;
    }


  //  recover the  args -----------------------------------------------
  string inputFileNameA = vm["imageA"].as<string>();
  string inputFileNameB = vm["imageB"].as<string>();
  string basenameOutput  =vm["imageError"].as<string>();


  
  //  Input images ----------------------------------------------------
  Image2D imageA = GenericReader<Image2D>::import(inputFileNameA);
  Image2D imageB = GenericReader<Image2D>::import(inputFileNameB);

  
  Image2DErr imageErr (imageA.domain());  


  //  Absolute Error between imageA and imageB -------------------------
  Statistic<int> statMA = getMAEstats(imageA, imageB, imageErr); 
  int maxVal = statMA.max();
  if(vm.count("setMaxColorValueMAE"))
  {
    maxVal = vm["setMaxColorValueMAE"].as<int>();
  }
  JetMap jmapMA(0, maxVal);
  displayStats(statMA, "Absolute errror");
  stringstream maeName; maeName << basenameOutput;
  maeName << "MAE.ppm";
  PPMWriter<Image2DErr, JetMap>::exportPPM(maeName.str(), imageErr,  jmapMA);


  //  Squared Error between imageA and imageB -------------------------
  Statistic<int> statSE = getMSEstats(imageA, imageB, imageErr); 
  maxVal = statMA.max();
  if(vm.count("fixMaxColorValueMSE"))
  {
    maxVal = vm["fixMaxColorValueMSE"].as<int>();
  }    
  JetMap jmapSE(0, maxVal);
  displayStats(statSE, "Squared errror");
  stringstream mseName; mseName << basenameOutput;
  mseName << "MSE.ppm";
  PPMWriter<Image2DErr, JetMap>::exportPPM(mseName.str(), imageErr,  jmapSE);

 

  return 0;
}

