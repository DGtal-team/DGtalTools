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
 * @file 3dImageViewer.cpp
 * @ingroup visualisation
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

#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/SurfelAdjacency.h"

#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#ifdef WITH_ITK
#include "DGtal/io/readers/DicomReader.h"
#endif

#include "DGtal/images/ImageSelector.h"



#include "specificClasses/Viewer3DImage.cpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
   @page Doc3dImageViewer 3dImageViewer
 
   @brief Displays volume file as a voxel set by using QGLviewer.

   @b Usage:  3dImageViewer [options] input

   @b Allowed @b options @b are :
 
   @code
   -h [ --help ]                    display this message
   -i [ --input ] arg               vol file (.vol, .longvol .p3d, .pgm3d and if
                                    WITH_ITK is selected: dicom, dcm, mha, mhd) 
                                    or sdp (sequence of discrete points). For 
                                    longvol, dicom, dcm, mha or mhd formats, the 
                                    input values are linearly scaled between 0 and 255.
   --grid                           draw slice images using grid mode. 
   --intergrid                      draw slice images using inter grid mode. 
   --emptyMode                      remove the default boundingbox display 
   --thresholdImage                 threshold the image to define binary shape
   -m [ --thresholdMin ] arg (=0)   threshold min to define binary shape
   -M [ --thresholdMax ] arg (=255) threshold max to define binary shape
   -s [ --displaySDP ] arg          display a set of discrete points (.sdp)
   --SDPindex arg                   specify the sdp index (by default 0,1,2).
   --SDPball arg (=0.5)             use balls to display a set of discrete 
   points (used with displaySDP option)
   --displayMesh arg                display a Mesh given in OFF or OFS format. 
   --displayDigitalSurface          display the digital surface instead of 
   display all the set of voxels (used with 
   thresholdImage or displaySDP options)
   --colorizeCC                     colorize each Connected Components of the 
   surface displayed by displayDigitalSurface 
   option.
   -c [ --colorSDP ] arg            set the color  discrete points: r g b a 
   --colorMesh arg                  set the color of Mesh (given from 
   displayMesh option) : r g b a 
   -x [ --scaleX ] arg (=1)         set the scale value in the X direction 
   (default 1.0)
   -y [ --scaleY ] arg (=1)         set the scale value in the Y direction 
   (default 1.0)
   -z [ --scaleZ ] arg (=1)         set the scale value in the Z direction 
   (default 1.0)
  --rescaleInputMin arg (=0)       min value used to rescale the input 
                                   intensity (to avoid basic cast into 8  bits 
                                   image).
  --rescaleInputMax arg (=255)     max value used to rescale the input 
                                   intensity (to avoid basic cast into 8 bits 
                                   image).
   -t [ --transparency ] arg (=255) transparency
   @endcode


   @b Example: 
   With the image display you can also threshold the image and display a set of voxel:  
   @code
   3dImageViewer -i $DGtal/examples/samples/lobster.vol --thresholdImage -m 180
   @endcode

   You should obtain such a result:

   @image html res3dImageViewer.png "resulting visualisation of 3d image with thresholded set of voxels."
 

   @see
   @ref 3dImageViewer.cpp

*/

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{

  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain,  unsigned char > Image3D;
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain,  unsigned char > Image2D;


  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd) or sdp (sequence of discrete points). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
    ("grid", "draw slice images using grid mode. " )
    ("intergrid", "draw slice images using inter grid mode. " )
    ("emptyMode", "remove the default boundingbox display " )
    ("thresholdImage", "threshold the image to define binary shape" )
    ("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min to define binary shape" )
    ("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max to define binary shape" )
    ("displaySDP,s", po::value<std::string>(), "display a set of discrete points (.sdp)" )
    ("SDPindex", po::value<std::vector <unsigned int> >()->multitoken(), "specify the sdp index (by default 0,1,2).")
    ("SDPball", po::value<double>()->default_value(0.5), "use balls to display a set of discrete points (used with displaySDP option)")
    ("displayMesh", po::value<std::string>(), "display a Mesh given in OFF or OFS format. " )
    ("displayDigitalSurface", "display the digital surface instead of display all the set of voxels (used with thresholdImage or displaySDP options)" )
    ("colorizeCC", "colorize each Connected Components of the surface displayed by displayDigitalSurface option." )
    ("colorSDP,c", po::value<std::vector <int> >()->multitoken(), "set the color  discrete points: r g b a " )
    ("colorMesh", po::value<std::vector <int> >()->multitoken(), "set the color of Mesh (given from displayMesh option) : r g b a " )
    ("scaleX,x",  po::value<float>()->default_value(1.0), "set the scale value in the X direction (default 1.0)" )
    ("scaleY,y",  po::value<float>()->default_value(1.0), "set the scale value in the Y direction (default 1.0)" )
    ("scaleZ,z",  po::value<float>()->default_value(1.0), "set the scale value in the Z direction (default 1.0)")
    ("rescaleInputMin", po::value<DGtal::int64_t>()->default_value(0), "min value used to rescale the input intensity (to avoid basic cast into 8  bits image).")
    ("rescaleInputMax", po::value<DGtal::int64_t>()->default_value(255), "max value used to rescale the input intensity (to avoid basic cast into 8 bits image).")
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
      std::cout << "Usage: " << argv[0] << " [input]\n"
                << "Displays volume file as a voxel set by using QGLviewer"
                << general_opt << "\n";
      return 0;
    }

  if(! vm.count("input"))
    {
      trace.error() << " The file name was defined" << endl;
      return 0;
    }
  string inputFilename = vm["input"].as<std::string>();
  int thresholdMin = vm["thresholdMin"].as<int>();
  int thresholdMax = vm["thresholdMax"].as<int>();
  unsigned char transp = vm["transparency"].as<uint>();

  QApplication application(argc,argv);
  
  float sx = vm["scaleX"].as<float>();
  float sy = vm["scaleY"].as<float>();
  float sz = vm["scaleZ"].as<float>();

  double ballRadius = vm["SDPball"].as<double>();
  string extension = inputFilename.substr(inputFilename.find_last_of(".") + 1);
 
  Viewer3DImage<>::ModeVisu mode;
  if(vm.count("emptyMode"))
    mode=Viewer3DImage<>::Empty;
  else if(vm.count("grid"))
    mode=Viewer3DImage<>::Grid;
  else if(vm.count("intergrid"))
    mode=Viewer3DImage<>::InterGrid;
  else
    mode=Viewer3DImage<>::BoundingBox;

  Viewer3DImage<> viewer(mode);
  viewer.setWindowTitle("simple Volume Viewer");
  viewer.show();
  viewer.setGLScale(sx, sy, sz);

  DGtal::int64_t rescaleInputMin = vm["rescaleInputMin"].as<DGtal::int64_t>();
  DGtal::int64_t rescaleInputMax = vm["rescaleInputMax"].as<DGtal::int64_t>();

  typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
  Image3D image =  GenericReader< Image3D >::importWithValueFunctor( inputFilename,RescalFCT(rescaleInputMin,
                                                                                             rescaleInputMax,
                                                                                             0, 255) );
  Domain domain = image.domain();

  trace.info() << "Image loaded: "<<image<< std::endl;
  viewer.setVolImage(&image);
  
  // Used to display 3D surface
  Z3i::DigitalSet set3d(domain);



  viewer << Viewer3D<>::updateDisplay;
  if(vm.count("thresholdImage")){
    GradientColorMap<long> gradient( thresholdMin, thresholdMax);
    gradient.addColor(Color::Blue);
    gradient.addColor(Color::Green);
    gradient.addColor(Color::Yellow);
    gradient.addColor(Color::Red);
    for(Domain::ConstIterator it = domain.begin(), itend=domain.end(); it!=itend; ++it){
      unsigned char  val= image( (*it) );
      Color c= gradient(val);
      if(val<=thresholdMax && val >=thresholdMin){
        if(!vm.count("displayDigitalSurface")){
          viewer <<  CustomColors3D(Color((float)(c.red()), (float)(c.green()),(float)(c.blue()), transp),
                                    Color((float)(c.red()), (float)(c.green()),(float)(c.blue()), transp));
          viewer << *it;
        }
      }else{
        set3d.insert(*it);
      }
    }
  }

  if(vm.count("displaySDP")){
    if(vm.count("colorSDP")){
      std::vector<int> vcol= vm["colorSDP"].as<std::vector<int > >();
      if(vcol.size()<4){
        trace.error() << "Not enough parameter: color specification should contains four elements: red, green, blue and alpha values."  << std::endl;
        return 0;
      }
      Color c(vcol[0], vcol[1], vcol[2], vcol[3]);
      viewer << CustomColors3D(c, c);
    }

    vector<Z3i::Point> vectVoxels;
    if(vm.count("SDPindex")) {
      std::vector<unsigned int > vectIndex = vm["SDPindex"].as<std::vector<unsigned int > >();
      if(vectIndex.size()!=3){
        trace.error() << "you need to specify the three indexes of vertex." << std::endl;
        return 0;
      }
      vectVoxels = PointListReader<Z3i::Point>::getPointsFromFile(vm["displaySDP"].as<std::string>(), vectIndex);
    }else{
      vectVoxels = PointListReader<Z3i::Point>::getPointsFromFile(vm["displaySDP"].as<std::string>());
    }
    for(unsigned int i=0;i< vectVoxels.size(); i++){
      if(!vm.count("displayDigitalSurface")){
        if(vm.count("SDPball")){
          viewer.addBall (vectVoxels.at(i), ballRadius);
        }else{
          viewer << vectVoxels.at(i);
        }
      }else{
        set3d.insert(vectVoxels.at(i));
      }
    }
  }

  if(vm.count("displayMesh")){
    if(vm.count("colorMesh")){
      std::vector<int> vcol= vm["colorMesh"].as<std::vector<int > >();
      if(vcol.size()<4){
        trace.error() << "Not enough parameter: color specification should contains four elements: red, green, blue and alpha values." << std::endl;
        return 0;
      }
      Color c(vcol[0], vcol[1], vcol[2], vcol[3]);
      viewer.setFillColor(c);
    }

    DGtal::Mesh<Z3i::RealPoint> aMesh(!vm.count("colorMesh"));
    MeshReader<Z3i::RealPoint>::importOFFFile(vm["displayMesh"].as<std::string>(), aMesh);
    viewer << aMesh;
  }

  if(vm.count("displayDigitalSurface")){
    KSpace K;
    Point low = domain.lowerBound(); low[0]=low[0]-1; low[1]=low[1]-1; low[2]=low[2]-1;
    Point upp = domain.upperBound(); upp[0]=upp[0]+1; upp[1]=upp[1]+1; upp[2]=upp[2]+1;
    K.init(low, upp , true);
    SurfelAdjacency<3> SAdj( true );
    vector<vector<SCell> > vectConnectedSCell;
    trace.info() << "Extracting surface  set ... " ;
    Surfaces<KSpace>::extractAllConnectedSCell(vectConnectedSCell,K, SAdj, set3d, true);
    trace.info()<< " [done] " <<std::endl;
    GradientColorMap<long> gradient( 0, vectConnectedSCell.size());
    gradient.addColor(DGtal::Color::Red);
    gradient.addColor(DGtal::Color::Yellow);
    gradient.addColor(DGtal::Color::Green);
    gradient.addColor(DGtal::Color::Cyan);
    gradient.addColor(DGtal::Color::Blue);
    gradient.addColor(DGtal::Color::Magenta);
    gradient.addColor(DGtal::Color::Red);

    viewer << DGtal::SetMode3D(vectConnectedSCell.at(0).at(0).className(), "Basic");
    for(unsigned int i= 0; i <vectConnectedSCell.size(); i++){
      for(unsigned int j= 0; j <vectConnectedSCell.at(i).size(); j++){
        if(vm.count("colorizeCC")){
          DGtal::Color c= gradient(i);
          viewer << CustomColors3D(Color(250, 0,0, transp), Color(c.red(),
                                                                  c.green(),
                                                                  c.blue(), transp));
        }else  if(vm.count("colorSDP")){
          std::vector<int> vcol= vm["colorSDP"].as<std::vector<int > >();
          Color c(vcol[0], vcol[1], vcol[2], vcol[3]);
          viewer << CustomColors3D(c, c);
        }

        viewer << vectConnectedSCell.at(i).at(j);
      }
    }
  }
  
  viewer << Viewer3D<>::updateDisplay;
  DGtal::Z3i::Point size = image.domain().upperBound() - image.domain().lowerBound();
  DGtal::Z3i::Point center = image.domain().lowerBound()+size/2;
  unsigned int maxDist = std::max(std::max(size[2], size[1]), size[0]);
  viewer.camera()->setPosition(qglviewer::Vec(center[0],center[1], 
                                              center[2] + 2.0*maxDist));
  viewer.camera()->setSceneCenter(qglviewer::Vec(center[0],center[1],center[2]));
  


  return application.exec();
}
