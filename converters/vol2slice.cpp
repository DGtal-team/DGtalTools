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
 * @file vol2slice.cpp
 * @ingroup Converters
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2013/05/07
 *
 *
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/format.hpp>



using namespace std;
using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

/**
 @page vol2slice vol2slice
 @brief  Convert a volumetric file (.vol, .longvol, .pgm3d) into a set of 2D slice  images.

@b Usage: vol2slice [input] [output]

@b Allowed @b options @b are:

@code
  -h [ --help ]                      display this message
  -i [ --input ] arg                 vol file (.vol, .longvol .p3d, .pgm3d and 
                                     if WITH_ITK is selected: dicom, dcm, mha, 
                                     mhd). For longvol, dicom, dcm, mha or mhd 
                                     formats, the input values are linearly 
                                     scaled between 0 and 255.
  -o [ --output ] arg                base_name.extension:  extracted 2D slice 
                                     volumetric files (will result n files 
                                     base_name_xxx.extension) 
  -f [ --setFirstSlice ] arg (=0)    Set the first slice index to be extracted.
  -l [ --setLastSlice ] arg          Set the last slice index to be extracted 
                                     (by default set to maximal value according
                                     to the given volume).
  -s [ --sliceOrientation ] arg (=2) specify the slice orientation for which 
                                     the slice are defined (by default =2 (Z 
                                     direction))
  --rescaleInputMin arg (=0)         min value used to rescale the input 
                                     intensity (to avoid basic cast into 8  
                                     bits image).
  --rescaleInputMax arg (=255)       max value used to rescale the input 
                                     intensity (to avoid basic cast into 8 bits
                                     image).

@endcode

@b Example:
@code 
   # Export Z slice images (-s 2): 
   $ vol2slice -i ${DGtal}/examples/samples/lobster.vol -o slice.pgm  -f 10 -l 15 -s 2
@endcode

You should obtain such a visualization:
@image html resVol2slice.png "resulting visualisation."

@see
@ref vol2slice.cpp

*/


int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char > Image2D;
  typedef DGtal::ConstImageAdapter<Image3D, Image2D::Domain, DGtal::functors::Projector< DGtal::Z3i::Space>,
				   Image3D::Value,  DGtal::functors::Identity >  SliceImageAdapter;


  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
    ("output,o", po::value<std::string>(), "base_name.extension:  extracted 2D slice volumetric files (will result n files base_name_xxx.extension) " )
    ("setFirstSlice,f", po::value<unsigned int>()->default_value(0), "Set the first slice index to be extracted.") 
    ("setLastSlice,l", po::value<unsigned int>(), "Set the last slice index to be extracted (by default set to maximal value according to the given volume).") 
    ("sliceOrientation,s", po::value<unsigned int>()->default_value(2), "specify the slice orientation for which the slice are defined (by default =2 (Z direction))" )
    ("rescaleInputMin", po::value<DGtal::int64_t>()->default_value(0), "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).")
    ("rescaleInputMax", po::value<DGtal::int64_t>()->default_value(255), "max value used to rescale the input intensity (to avoid basic cast into 8 bits image).");


  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);


  if( !parseOK || !  vm.count("input")||! vm.count("output") || vm.count("help"))
    {
      std::cout << "Usage: " << argv[0] << " [inputs] [output]\n"
		<< "Convert a volumetric file (.vol, .longvol, .pgm3d) into a set of 2D slice  images."
		<< general_opt << "\n";
      std::cout << "Example: to extract all slices defined in Y plane (y=cst): \n"
		<< "vol2slice -i image3d.vol -s 1  -o slice.pgm \n"
                << "see slice2vol"<< endl;
      return 0;
    }

  if(! vm.count("input")||! vm.count("output"))
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;
      return 0;
    }




  std::string inputFileName = vm["input"].as<std::string>();
  std::string outputFileName = vm["output"].as<std::string>();
  std::string outputExt = outputFileName.substr(outputFileName.find_last_of(".")+1);
  std::string outputBasename = outputFileName.substr(0, outputFileName.find_last_of("."));
  unsigned int sliceOrientation = vm["sliceOrientation"].as<unsigned int>();
  DGtal::int64_t rescaleInputMin = vm["rescaleInputMin"].as<DGtal::int64_t>();
  DGtal::int64_t rescaleInputMax = vm["rescaleInputMax"].as<DGtal::int64_t>();

  
  trace.info()<< "Importing volume file base name:  " << outputBasename << " extension: " << outputExt << " ..." ;
  typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
  Image3D input3dImage =  GenericReader< Image3D >::importWithValueFunctor( inputFileName,RescalFCT(rescaleInputMin,
                                                                                                    rescaleInputMax,
                                                                                                    0, 255) );
  
  trace.info()<< "[done]" << endl;
  
  unsigned int startSlice=0;
  unsigned int endSlice=input3dImage.domain().upperBound()[sliceOrientation];

  if(vm.count("setFirstSlice")){
    startSlice = vm["setFirstSlice"].as<unsigned int>();
  }
  if(vm.count("setLastSlice")){
    endSlice = vm["setLastSlice"].as<unsigned int>();    
  }

  //Processing each slice
#pragma omp parallel for schedule(dynamic)
   for( unsigned int i=startSlice; i <= endSlice; i++){
    trace.info() << "Exporting slice image "<< i ;
    DGtal::functors::Projector<DGtal::Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(sliceOrientation);
    DGtal::Z2i::Domain domain2D(invFunctor(input3dImage.domain().lowerBound()),
  				invFunctor(input3dImage.domain().upperBound()));
    DGtal::functors::Projector<DGtal::Z3i::Space> aSliceFunctor(i); aSliceFunctor.initAddOneDim(sliceOrientation);
    const DGtal::functors::Identity identityFunctor{};
    SliceImageAdapter sliceImage( input3dImage, domain2D, aSliceFunctor, identityFunctor );
    stringstream outName; outName << outputBasename << "_" <<  boost::format("%|05|")% i <<"."<< outputExt ;
    trace.info() << ": "<< outName.str() ;
    GenericWriter<SliceImageAdapter>::exportFile(outName.str(), sliceImage);
    trace.info() << " [done]"<< endl;
  }


  return 0;
}
