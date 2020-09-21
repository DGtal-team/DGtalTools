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

#include "CLI11.hpp"


///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace DGtal;
///////////////////////////////////////////////////////////////////////////////

typedef ImageContainerBySTLVector < Z2i::Domain,  unsigned char > Image2D;
typedef ImageContainerBySTLVector < Z2i::Domain,  unsigned int > Image2DErr;
typedef GradientColorMap<unsigned int, CMAP_JET, 1 > JetMap;



/**
 @page Doc2dCompImage 2dCompImage
 
 @brief Compare images and displays differences (squared and absolute differences). 


 @b Usage:  2dCompImage --imageA <imageA>.pgm --imageB <imageB>.pgm --imageError <name> 

 @b Allowed @b options @b are :
 
 @code

 Positionals:
   1 TEXT:FILE REQUIRED                  Input filename of image A.
   2 TEXT:FILE REQUIRED                  Input filename of image B.

 Options:
   -h,--help                             Print this help message and exit
   -a,--imageA TEXT:FILE REQUIRED        Input filename of image A.
   -b,--imageB TEXT:FILE REQUIRED        Input filename of image B.
   -e,--imageError TEXT                  Output error image basename (will generate two images <basename>MSE.ppm and <basename>MAE.ppm).
   -S,--setMaxColorValueMSE INT          Set the maximal color value for the scale display of MSE (else the scale is set the maximal MSE value).
   -A,--setMaxColorValueMAE INT          Set the maximal color value for the scale display of MAE (else the scale is set from the maximal MAE value).
   
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
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileNameA;
  std::string inputFileNameB;
  std::string basenameOutput;
  int maxValueMAE;
  int maxValueMSE;
  
  app.description("Compare images and displays differences (squared and absolute differences).\n Typical use example:\n \t 2dCompImage -a imageA.pgm -b imageB.pgm -e errorImage -S 100 \n");
  app.add_option("-a,--imageA,1", inputFileNameA, "Input filename of image A." )
  ->required()
  ->check(CLI::ExistingFile);
  app.add_option("-b,--imageB,2", inputFileNameA, "Input filename of image B." )
   ->required()
   ->check(CLI::ExistingFile);

  
  app.add_option("--imageError,-e",basenameOutput, "Output error image basename (will generate two images <basename>MSE.ppm and <basename>MAE.ppm).");
  
  auto setMaxMSEOpt = app.add_option("--setMaxColorValueMSE,-S",maxValueMSE, "Set the maximal color value for the scale display of MSE (else the scale is set the maximal MSE value).");
  auto setMaxMAEOpt = app.add_option("--setMaxColorValueMAE,-A",maxValueMAE, "Set the maximal color value for the scale display of MAE (else the scale is set from the maximal MAE value).");

  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  //  Input images ----------------------------------------------------
  Image2D imageA = GenericReader<Image2D>::import(inputFileNameA);
  Image2D imageB = GenericReader<Image2D>::import(inputFileNameB);
  
  Image2DErr imageErr (imageA.domain());  


  //  Absolute Error between imageA and imageB -------------------------
  Statistic<int> statMA = getMAEstats(imageA, imageB, imageErr); 
  int maxVal = statMA.max();
  if(setMaxMAEOpt->count() > 0)
  {
    maxVal = maxValueMAE;
  }
  JetMap jmapMA(0, maxVal);
  displayStats(statMA, "Absolute errror");
  stringstream maeName; maeName << basenameOutput;
  maeName << "MAE.ppm";
  PPMWriter<Image2DErr, JetMap>::exportPPM(maeName.str(), imageErr,  jmapMA);

  //  Squared Error between imageA and imageB -------------------------
  Statistic<int> statSE = getMSEstats(imageA, imageB, imageErr); 
  maxVal = statMA.max();
  if(setMaxMSEOpt->count()>0)
  {
    maxVal = maxValueMSE;
  }    
  JetMap jmapSE(0, maxVal);
  displayStats(statSE, "Squared error");
  stringstream mseName; mseName << basenameOutput;
  mseName << "MSE.ppm";
  PPMWriter<Image2DErr, JetMap>::exportPPM(mseName.str(), imageErr,  jmapSE);

  return 0;
}

