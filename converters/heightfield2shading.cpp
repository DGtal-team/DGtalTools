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
 * @file heightfield2vol.cpp
 * @ingroup converters
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2015/03/18
 *
 * 
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

template<typename TImage, typename TImageVector>
void 
computerBasicNormalsFromHeightField(const TImage &anHeightMap, TImageVector &vectorField)
{
  for(typename TImage::Domain::ConstIterator it = anHeightMap.domain().begin(); 
      it != anHeightMap.domain().end(); it++){
    if(anHeightMap.domain().isInside(*it+Z2i::Point::diagonal(1))&&
       anHeightMap.domain().isInside(*it-Z2i::Point::diagonal(1))){
      double dx = (anHeightMap(*it-Z2i::Point(1,0))-anHeightMap(*it+Z2i::Point(1,0)))/2.0;
      double dy = (anHeightMap(*it-Z2i::Point(0,1))-anHeightMap(*it+Z2i::Point(0,1)))/2.0;
      Z3i::RealPoint n (dx, dy, 1);
      n /= n.norm();
      vectorField.setValue(*it,n);
    }
  }
}



template<typename TImageVector>
void 
importNormals(std::string file, TImageVector &vectorField)
{
  std::vector<Z3i::RealPoint> vp = PointListReader<Z3i::RealPoint>::getPointsFromFile(file);
   for(unsigned int i = 0; i< vp.size(); i=i+2){
        Z3i::RealPoint p = vp.at(i);
        Z3i::RealPoint q = vp.at(i+1);
        Z3i::RealPoint n = (q-p)/(p-q).norm();
        vectorField.setValue(Z2i::Point(p[0], p[1]),n);
   }
}


// Defining an helper to get the 3D point functor from an 2DImage 
template<typename TImage2D, typename TPoint3D >
struct LambertianShadindFunctor{
  typedef  TPoint3D Point3D;
  typedef typename TImage2D::Value Value;
  
  /**
   *  Construct the predicat given a normal vector
   **/
  LambertianShadindFunctor(const Point3D &lightSourceDirection):
    myLightSourceDirection(lightSourceDirection/lightSourceDirection.norm()){}
  
  inline
  unsigned int operator()(const Point3D &aNormal)  const {
    //
    int intensity = aNormal.dot(myLightSourceDirection)*std::numeric_limits<typename TImage2D::Value>::max();
    return intensity>0? intensity:0;
  }
  
  Point3D  myLightSourceDirection;
};



int main( int argc, char** argv )
{
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char> Image2D;
  typedef ImageContainerBySTLVector < Z2i::Domain, Z3i::RealPoint> Image2DNormals;
  
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "heightfield file." )
    ("output,o", po::value<std::string>(), "output image.") 
    ("importNormal", po::value<std::string>(), "import normals from file.") 
    ("lx", po::value<double>(), "x light source direction.") 
    ("ly", po::value<double>(), "y light source direction." )
    ("lz", po::value<double>(), "z light source direction.");
    
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);    
  if( !parseOK || vm.count("help")||argc<=1)
    {
      std::cout << "Usage: " << argv[0] << " [input] [output]\n"
		<< "Render a 2D heightfield image into a shading image."
		<< general_opt << "\n";
      std::cout << "Example:\n"
		<< "heightfield2shading -i ${DGtal}/examples/samples/church.pgm -o volResu.vol -s 0.3 -z 50  \n";
      return 0;
    }
  
  if(! vm.count("input") ||! vm.count("output"))
    {
      trace.error() << " Input and output filename are needed to be defined" << endl;      
      return 0;
    }
  
  string inputFilename = vm["input"].as<std::string>();
  string outputFilename = vm["output"].as<std::string>();
  double lx = vm["lx"].as<double>();
  double ly = vm["ly"].as<double>();
  double lz = vm["lz"].as<double>();
  LambertianShadindFunctor<Image2D, Z3i::RealPoint> lShade (Z3i::RealPoint(lx,ly,lz));
  
  trace.info() << "Reading input file " << inputFilename ; 
  Image2D inputImage = DGtal::GenericReader<Image2D>::import(inputFilename);  
  Image2DNormals vectNormals (inputImage.domain());
  Image2D result (inputImage.domain());
  if(vm.count("importNormal")){
    std::string normalFileName = vm["importNormal"].as<string>();
    importNormals(normalFileName, vectNormals);
  }else{
    computerBasicNormalsFromHeightField(inputImage, vectNormals);
  }
  for(typename Image2D::Domain::ConstIterator it = inputImage.domain().begin(); 
      it != inputImage.domain().end(); it++){
    result.setValue(*it, lShade(vectNormals(*it))); 
                    
  }

  result >> outputFilename;
  return 0;  
}



