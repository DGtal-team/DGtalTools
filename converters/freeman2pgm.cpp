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
 * @file freeman2pgm.cpp
 * @ingroup Tools
 * @author Bertrand Kerautret (\c kerautre@loria.fr) and Jacques-Olivier Lachaud 
 * LORIA (CNRS, UMR 7503), University of Nancy, France 
 * (backport from ImaGene)
 * @date 2012/19/05
 *
 * convert freeman chain to a Sequence of Discrete Points.  
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

//image
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/writers/GenericWriter.h"

//contour
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/topology/helpers/Surfaces.h"

//boost
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

//STL
#include <vector>
#include <string>

using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")  
    ("FreemanChain,f", po::value<std::string>(), "FreemanChain file name")
    ("output,o", po::value<std::string>()->default_value("result.pgm"), " the output filename")
    ("space,s", po::value<std::vector <int> >()->multitoken(), "Define the space from its bounding box (lower and upper coordinates) \
else the space is automatically defined from the freemanchain bounding boxes.");
  
  typedef KhalimskySpaceND<2, int>::Cell Cell;
  typedef KhalimskySpaceND<2, int>::SCell SCell;
  typedef FreemanChain<Z2i::Integer> FreemanChain;
  typedef DGtal::KhalimskySpaceND< 2, int > KSpace;
  typedef DGtal::ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image2D ; 
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }
  
 
  po::notify(vm);    
  if(!parseOK||vm.count("help")||argc<=1 || (!(vm.count("FreemanChain")) ) )
    {
      trace.info()<< "Transform a freeman chain into a pgm file." <<std::endl << "Basic usage: "<<std::endl
		  << "\t freeman2pgm [options] --FreemanChain  <fileName>  "<<std::endl
		  << general_opt << "\n";
      return 0;
    }  
  
  if( vm.count("FreemanChain") ){
    std::string fileName = vm["FreemanChain"].as<std::string>();
    std::vector< FreemanChain > vectFcs =  PointListReader< Z2i::Point >::getFreemanChainsFromFile<Z2i::Integer> (fileName);    
    int minx, miny, maxx, maxy;
    if(!vm.count("space")){
      for(std::vector< FreemanChain >::const_iterator it = vectFcs.begin(); it!= vectFcs.end(); it++){
        FreemanChain fc = *it;
        int t_minx, t_miny, t_maxx, t_maxy;
        fc.computeBoundingBox(t_minx, t_miny, t_maxx, t_maxy);
        minx = t_minx > minx? minx: t_minx;
        miny = t_miny > miny? miny: t_miny;
        maxx = t_maxx < maxx? maxx: t_maxx;
        maxy = t_maxy < maxy? maxy: t_maxy;
      }
    }else{
      std::vector<int> vectSpace = vm["space"].as<std::vector<int> > ();
      if(vectSpace.size()!=4){
        trace.error() << " Option : --space: you need to enter the two lower and upper point of the space."<<std::endl;
        return 0;
      }
      minx = vectSpace[0];
      miny = vectSpace[1];
      maxx = vectSpace[2];
      maxy = vectSpace[3];
      
    }
    KSpace aKSpace;
    aKSpace.init(Z2i::Point(minx, miny), Z2i::Point(maxx, maxy), false);
    std::set<SCell> boundarySCell;
    std::set<Cell> interiorCell;
    for(std::vector< FreemanChain >::const_iterator it = vectFcs.begin(); it!= vectFcs.end(); it++){
      FreemanChain fc = *it;
      FreemanChain::getContourSCell(aKSpace, fc, boundarySCell, true); 
    }
    
    Surfaces<KSpace>::uComputeInterior(aKSpace, boundarySCell, interiorCell, false);  
    Image2D imageResult (Z2i::Domain(Z2i::Point(minx, miny), Z2i::Point(maxx, maxy))); 
    for(std::set<Cell>::const_iterator it = interiorCell.begin(); it!=interiorCell.end(); it++){
      imageResult.setValue(aKSpace.uCoords(*it), 255);
    }
    imageResult >> vm["output"].as<std::string>();
  }

}

