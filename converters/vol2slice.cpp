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
 * @file slice2vol.cpp
 * @ingroup surfaceTools
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



using namespace std;
using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z3i::Domain, unsigned char > Image3D;
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char > Image2D;
  typedef DGtal::ConstImageAdapter<Image3D, Image2D::Domain, DGtal::Projector< DGtal::Z3i::Space>,
				   Image3D::Value,  DGtal::DefaultFunctor >  SliceImageAdapter;


  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input-file,i", po::value<std::string >(), "input volumetric file (.vol, .longvol, .pgm3d) " )
    ("output-files,o", po::value<std::string>(), "base_name.extension:  extracted 2D slice volumetric files (will result n files base_name_xxx.extension) " )
    ("sliceOrientation,s", po::value<unsigned int>()->default_value(2), "specify the slice orientation for which the slice are defined (by default =2 (Z direction))" );


  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);


  if( !parseOK || vm.count("help"))
    {
      std::cout << "Usage: " << argv[0] << " [input-files] [output-file]\n"
		<< "Convert a volumetric file (.vol, .longvol, .pgm3d) into a set of 2D slice  images."
		<< general_opt << "\n";
      std::cout << "Example: to extract all slices defined in Y plane (y=cst): \n"
		<< "vol2slice -i image3d.vol -s 1  -o slice.pgm \n"
                << "see slice2vol"<< endl;
      return 0;
    }

  if(! vm.count("input-file")||! vm.count("output-files"))
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;
      return 0;
    }




  std::string inputFileName = vm["input-file"].as<std::string>();
  std::string outputFileName = vm["output-files"].as<std::string>();
  std::string outputExt = outputFileName.substr(outputFileName.find_last_of(".")+1);
  std::string outputBasename = outputFileName.substr(0, outputFileName.find_last_of("."));
  unsigned int sliceOrientation = vm["sliceOrientation"].as<unsigned int>();

  trace.info()<< "Importing volume file base name:  " << outputBasename << " extension: " << outputExt << " ..." ;
  Image3D input3dImage = GenericReader<Image3D>::import(inputFileName);
  trace.info()<< "[done]" << endl;


  //Processing each slice
#pragma omp parallel for schedule(dynamic)
  for( unsigned int i=0; i <= input3dImage.domain().upperBound()[sliceOrientation]; i++){
    trace.info() << "Exporting slice image "<< i ;
    DGtal::Projector<DGtal::Z2i::Space>  invFunctor; invFunctor.initRemoveOneDim(sliceOrientation);
    DGtal::Z2i::Domain domain2D(invFunctor(input3dImage.domain().lowerBound()),
  				invFunctor(input3dImage.domain().upperBound()));
    DGtal::Projector<DGtal::Z3i::Space> aSliceFunctor(i); aSliceFunctor.initAddOneDim(sliceOrientation);
    SliceImageAdapter sliceImage(input3dImage, domain2D, aSliceFunctor, DGtal::DefaultFunctor());
    stringstream outName; outName << outputBasename << "_" <<  boost::format("%|05|")% i <<"."<< outputExt ;
    trace.info() << ": "<< outName.str() ;
    GenericWriter<SliceImageAdapter>::exportFile(outName.str(), sliceImage);
    trace.info() << " [done]"<< endl;
  }


  return 0;
}
