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
 * @file volCrop.cpp
 * @ingroup volumetric/voltools
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2013/08/06
 *
 * 
 *
 * This file is part of the DGtalTools.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ConstImageAdapter.h>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


using namespace std;
using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;



/**
 @page volFlip volFlip
 
 @brief  Flips 2D slice image of an 3D vol image (mirror transformation).

 @b Usage:  volFlip --input \<volFileName\> --imagePlane 0 1 --flipDimension 0 --o \<volOutputFileName\> (vol, longvol, p3d format)



 @b Allowed @b options @b are : 
 @code
  -h [ --help ]                     display this message.
  -i [ --input ] arg                Input vol file.
  --imagePlane arg                  arg=  {0,1,2} x {0,1,2} defines the axis of
                                    the slice image which will be transformed 
                                    (by default arg= 0 1  i.e. the slice image 
                                    defined in the X,Y plane (Z=cst)
  --flipDimension arg               specify which axis will be used to apply 
                                    the flip.
  -o [ --output ] arg (=output.vol) Output filename.
 @endcode

 @b Example: 

 @code
 $ volFlip --imagePlane 0 1 --flipDimension 0 -i ${DGtal}/examples/samples/lobster.vol -o flippedXxyLobster.vol 
 @endcode


 You should obtain such a result:
 @image html resVolFlip.png "(a) source image (b) flipped version."
 
 @see
 @ref volFlip.cpp

 */


/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam ( std::string param )
{
  trace.error() <<" Parameter: "<<param<<" is required..";
  trace.info() <<std::endl;
  exit ( 1 );
}


int main(int argc, char**argv)
{

  // parse command line ----------------------------------------------
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<std::string>(), "Input vol file." )
    ( "imagePlane",po::value<std::vector <unsigned int> >()->multitoken(),
      "arg=  {0,1,2} x {0,1,2} defines the axis of the slice image which will be transformed (by default arg= 0 1  i.e. the slice image defined in the X,Y plane (Z=cst)" )
    ( "flipDimension", po::value<unsigned int>(), "specify which axis will be used to apply the flip." )
    ( "output,o", po::value<string>()->default_value("output.vol"),"Output filename." );
  bool parseOK=true;
  
  po::variables_map vm;
  
  
  try{
    po::store ( po::parse_command_line ( argc, argv, general_opt ), vm );
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify ( vm );
  if (!parseOK || ! ( vm.count ( "input" ) ) || ! ( vm.count ( "output" ) ) || ! ( vm.count ( "imagePlane" ) ) 
      || ! ( vm.count ( "flipDimension" ) ) || vm.count ( "help" ))
    {
      trace.info() << "Flip 2D slice image of an 3D vol image (mirror transformation)"<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\t volFlip --input <volFileName> --imagePlane 0 1 --flipDimension 0 --o <volOutputFileName> (vol, longvol, p3d format)"<<std::endl
                   << general_opt << "\n";
      std::cout << "Example:\n"
		<< "volFlip --imagePlane 0 1 --flipDimension 0 -i ${DGtal}/examples/samples/lobster.vol -o flippedXxyLobster.p3d \n The resulting Z slice images (Z= cst) of flippedXxyLobster.p3d will appears flipped according the x axis.  ";
      return 0;
    }


  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
  if ( ! ( vm.count ( "imagePlane" ) ) || vm["imagePlane"].as<std::vector<unsigned int > >().size()!=2 ) missingParam ( "--imagePlane" );
  if ( ! ( vm.count ( "flipDimension" ) ) ) missingParam ( "--flipDimension" );


  std::string inputFilename = vm["input"].as<std::string>();
  std::string outputFileName = vm["output"].as<std::string>();
  
  unsigned int dimFirstImg = vm["imagePlane"].as<std::vector<unsigned int > >().at(0);
  unsigned int dimSecondImg = vm["imagePlane"].as<std::vector<unsigned int > >().at(1);
  unsigned int dimFlip = vm["flipDimension"].as<unsigned int>();
  
  unsigned int normalImgDim = (dimFirstImg!=2 && dimSecondImg!=2)? 2 :( (dimFirstImg!=1 && dimSecondImg!=1)? 1:  0    );  
    
    
  trace.beginBlock("Loading file");
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  Image3D;
  Image3D  imageSRC =  GenericReader<Image3D>::import ( inputFilename );
  trace.endBlock();
  Image3D  imageRes(imageSRC.domain());
  for(int i=0; i <= imageSRC.domain().upperBound()[normalImgDim]; i++){
    Point startPoint(0,0, 0);
    startPoint[normalImgDim]=i;
    for( Domain::ConstSubRange::ConstIterator 
	   it = imageSRC.domain().subRange(dimFirstImg, dimSecondImg, 
					   startPoint).begin(),
	   itend =  imageSRC.domain().subRange(dimFirstImg, dimSecondImg, startPoint).end();
	 it != itend; ++it){
      Point pt = *it;
      pt[dimFlip]= imageSRC.domain().upperBound()[dimFlip] - pt[dimFlip] ;
      imageRes.setValue(*it, imageSRC(pt)); 
    }
  }  
  

  trace.beginBlock("Exporting...");
  bool res =  VolWriter< Image3D>::exportVol(outputFileName, imageRes);
  trace.endBlock();
  if (res) return 0; else return 1;
}




