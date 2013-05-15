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
 * @file visuDistanceTransform.cpp
 * @ingroup Examples
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2011/01/04
 *
 * An example file named qglViewer.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <QtGui/qapplication.h>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/readers/PointListReader.h"


#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/PNMReader.h"
#include "DGtal/images/ImageSelector.h"



#include "specificClasses/Viewer3DImage.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


using namespace std;
using namespace DGtal;
using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{

  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain,  unsigned char > Image3D;
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain,  unsigned char > Image2D;

  
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input-file,i", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d) file or sdp (sequence of discrete points)" )
    ("grid", "draw slice images using grid mode. " ) 
    ("intergrid", "draw slice images using inter grid mode. " ) 
    ("boundingbox", "draw slice images using boudingbox " ) 
    ("thresholdImage", "threshold the image to define binary shape" ) 
    ("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min to define binary shape" ) 
    ("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max to define binary shape" )
    ("transparency,t",  po::value<uint>()->default_value(255), "transparency") ; 
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
      std::cout << "Usage: " << argv[0] << " [input-file]\n"
    << "Display volume file as a voxel set by using QGLviewer"
    << general_opt << "\n";
      return 0;
    }
  
  if(! vm.count("input-file"))
    {
      trace.error() << " The file name was defined" << endl;      
      return 0;
    }
  string inputFilename = vm["input-file"].as<std::string>();
  int thresholdMin = vm["thresholdMin"].as<int>();
  int thresholdMax = vm["thresholdMax"].as<int>();
  unsigned char transp = vm["transparency"].as<uint>();
 
  QApplication application(argc,argv);
 


  string extension = inputFilename.substr(inputFilename.find_last_of(".") + 1);
  if(extension!="vol" && extension != "p3d" && extension != "pgm3D" && extension != "pgm3d" && extension != "sdp"){
    trace.info() << "File extension not recognized: "<< extension << std::endl;
    return 0;
  }
  Viewer3DImage::ModeVisu mode;
  if(vm.count("boudingbox"))
    mode=Viewer3DImage::BoundingBox;
  else if(vm.count("grid"))
    mode=Viewer3DImage::Grid;
  else if(vm.count("intergrid"))
    mode=Viewer3DImage::InterGrid;
  else
    mode=Viewer3DImage::Empty;
  
  Viewer3DImage viewer(mode);
  viewer.setWindowTitle("simple Volume Viewer");
  viewer.show();
  
  
  if(extension=="vol" || extension=="pgm3d" || extension=="pgm3D"){
   Image3D image = (extension=="vol")? VolReader<Image3D>::importVol( inputFilename ): PNMReader<Image3D>::importPGM3D( inputFilename );
    trace.info() << "Image loaded: "<<image<< std::endl;
    viewer.setVolImage(&image);
    viewer << Z3i::Point(512, 512, 0);

    viewer << Viewer3D::updateDisplay;
    if(vm.count("thresholdImage")){
      Domain domain = image.domain();
      GradientColorMap<long> gradient( thresholdMin, thresholdMax);
      gradient.addColor(Color::Blue);
      gradient.addColor(Color::Green);
      gradient.addColor(Color::Yellow);
      gradient.addColor(Color::Red);
      for(Domain::ConstIterator it = domain.begin(), itend=domain.end(); it!=itend; ++it){
	unsigned char  val= image( (*it) );     
      
	Color c= gradient(val);
	if(val<=thresholdMax && val >=thresholdMin){
	  viewer <<  CustomColors3D(Color((float)(c.red()), (float)(c.green()),(float)(c.blue()), transp),
				    Color((float)(c.red()), (float)(c.green()),(float)(c.blue()), transp));     
	  viewer << *it;     
	}     
      }
      viewer<< Viewer3D::updateDisplay;
    }
  }else if(extension=="sdp"){
    vector<Z3i::Point> vectVoxels = PointListReader<Z3i::Point>::getPointsFromFile(inputFilename);
    for(int i=0;i< vectVoxels.size(); i++){
      viewer << vectVoxels.at(i);
    }


    
  }
  

  return application.exec();
}
