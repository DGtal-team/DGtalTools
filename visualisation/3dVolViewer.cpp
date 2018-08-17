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
 * @file 3dVolViewer.cpp
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Lorraine, France
 *
 * @date 2011/01/04
 *
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/readers/MeshReader.h"

#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/images/ImageSelector.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
   @page Doc3dVolViewer 3dVolViewer
 
   @brief Displays volume file as a voxel set by using QGLviewer.

   The mode  specifies if you wish to see surface elements (BDRY), the inner
   voxels (INNER) or the outer voxels (OUTER) that touch the boundary.

   @b Usage:   3dVolViewer [input]

   @b Allowed @b options @b are :
 
   @code
   -h [ --help ]                      display this message
   -i [ --input ] arg                sdp (sequence of discrete points)  or vol
                                     file (.vol, .longvol .p3d, .pgm3d and 
                                     if WITH_ITK is selected: dicom, dcm, mha, 
                                     mhd) or sdp (sequence of discrete points).
                                     For longvol, dicom, dcm, mha or mhd 
                                     formats, the input values are linearly 
                                     scaled between 0 and 255.
   -m [ --thresholdMin ] arg (=0)     threshold min to define binary shape
   -M [ --thresholdMax ] arg (=255)   threshold max to define binary shape
   -n [ --numMaxVoxel ] arg (=500000) set the maximal voxel number to be 
   displayed.
  --rescaleInputMin arg (=0)         min value used to rescale the input 
                                     intensity (to avoid basic cast into 8  
                                     bits image).
  --rescaleInputMax arg (=255)       max value used to rescale the input 
                                     intensity (to avoid basic cast into 8 bits
                                     image).
   --displayMesh arg                display a Mesh given in OFF or OFS format. 
   --colorMesh arg                  set the color of Mesh (given from 
   displayMesh option) : r g b a 
   -d [ --doSnapShotAndExit]  filename,  save display snapshot into file. Notes that the camera setting is set by default according the last saved configuration (use SHIFT+Key_M to save current camera setting in the Viewer3D). If the camera setting was not saved it will use the default camera setting.
   -t [ --transparency ] arg (=255)   transparency
   @endcode


   @b Example: 


   @code
   $ 3dVolViewer -i $DGtal/examples/samples/lobster.vol -m 60 -t 10
   @endcode

   You should obtain such a result:

   @image html res3dVolViewer.png "Resulting visualization."
 

   @see
   @ref 3dVolViewer.cpp

*/




template < typename Space = DGtal::Z3i::Space, typename KSpace = DGtal::Z3i::KSpace>
struct ViewerSnap: DGtal::Viewer3D <Space, KSpace>
{

  ViewerSnap(bool saveSnap): Viewer3D<Space, KSpace>(), mySaveSnap(saveSnap){
  };

  virtual  void
  init(){
    DGtal::Viewer3D<>::init();
    if(mySaveSnap){
      QObject::connect(this, SIGNAL(drawFinished(bool)), this, SLOT(saveSnapshot(bool)));
    }
  };
  bool mySaveSnap;
};


///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "vol file (.vol, .longvol .p3d, .pgm3d and if WITH_ITK is selected: dicom, dcm, mha, mhd) or sdp (sequence of discrete points). For longvol, dicom, dcm, mha or mhd formats, the input values are linearly scaled between 0 and 255." )
    ("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min to define binary shape" )
    ("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max to define binary shape" )
    ("numMaxVoxel,n",  po::value<int>()->default_value(500000), "set the maximal voxel number to be displayed." )
    ("displayMesh", po::value<std::string>(), "display a Mesh given in OFF or OFS format. " )
    ("colorMesh", po::value<std::vector <int> >()->multitoken(), "set the color of Mesh (given from displayMesh option) : r g b a " )
    ("doSnapShotAndExit,d", po::value<std::string>(), "save display snapshot into file. Notes that the camera setting is set by default according the last saved configuration (use SHIFT+Key_M to save current camera setting in the Viewer3D). If the camera setting was not saved it will use the default camera setting." )
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
                << "Display volume file as a voxel set by using QGLviewer"<< endl
                << general_opt << "\n"
                << "Example: "<< std::endl
                << "    \t 3dVolViewer -i $DGtal/examples/samples/lobster.vol -m 60 -t 10" << endl;
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

  bool limitDisplay=false;
  if(vm.count("numMaxVoxel")){
    limitDisplay=true;
  }
  unsigned int numDisplayedMax = vm["numMaxVoxel"].as<int>();


  QApplication application(argc,argv);
  typedef ViewerSnap<> Viewer;
  
  Viewer viewer(vm.count("doSnapShotAndExit"));
  if(vm.count("doSnapShotAndExit")){
    viewer.setSnapshotFileName(QString(vm["doSnapShotAndExit"].as<std::string>().c_str()));
  }

  viewer.setWindowTitle("simple Volume Viewer");
  viewer.show();

  typedef ImageSelector<Domain, unsigned char>::Type Image;
  string extension = inputFilename.substr(inputFilename.find_last_of(".") + 1);
  if(extension != "sdp")
    {
      unsigned int numDisplayed=0;
      DGtal::int64_t rescaleInputMin = vm["rescaleInputMin"].as<DGtal::int64_t>();
      DGtal::int64_t rescaleInputMax = vm["rescaleInputMax"].as<DGtal::int64_t>();
  
      typedef DGtal::functors::Rescaling<DGtal::int64_t ,unsigned char > RescalFCT;
      Image image =  GenericReader< Image >::importWithValueFunctor( inputFilename,RescalFCT(rescaleInputMin,
                                                                                             rescaleInputMax,
                                                                                             0, 255) );
  
  
      trace.info() << "Image loaded: "<<image<< std::endl;
      Domain domain = image.domain();
      GradientColorMap<long> gradient( thresholdMin, thresholdMax);
      gradient.addColor(Color::Blue);
      gradient.addColor(Color::Green);
      gradient.addColor(Color::Yellow);
      gradient.addColor(Color::Red);
      for(Domain::ConstIterator it = domain.begin(), itend=domain.end(); it!=itend; ++it){
        unsigned char  val= image( (*it) );
        if(limitDisplay && numDisplayed > numDisplayedMax)
          break;
        Color c= gradient(val);
        if(val<=thresholdMax && val >=thresholdMin){
          viewer <<  CustomColors3D(Color((float)(c.red()), (float)(c.green()),(float)(c.blue()), transp),
                                    Color((float)(c.red()), (float)(c.green()),(float)(c.blue()), transp));
          viewer << *it;
          numDisplayed++;
        }
      }
    }else if(extension=="sdp"){
    vector<Z3i::RealPoint> vectVoxels = PointListReader<Z3i::RealPoint>::getPointsFromFile(inputFilename);
    for(unsigned int i=0;i< vectVoxels.size(); i++){
      viewer << vectVoxels.at(i);
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

  viewer << Viewer3D<>::updateDisplay;
  if(vm.count("doSnapShotAndExit")){
    // Appy cleaning just save the last snap
    if(!viewer.restoreStateFromFile())
      {
        viewer.update();
      }    
    std::string name = vm["doSnapShotAndExit"].as<std::string>();
    std::string extension = name.substr(name.find_last_of(".") + 1);
    std::string basename = name.substr(0, name.find_last_of("."));
    for(int i=0; i< viewer.snapshotCounter()-1; i++){
      std::stringstream s;
      s << basename << "-"<< setfill('0') << setw(4)<<  i << "." << extension;
      trace.info() << "erase temp file: " << s.str() << std::endl;
      remove(s.str().c_str());
    }
    std::stringstream s;
    s << basename << "-"<< setfill('0') << setw(4)<<  viewer.snapshotCounter()-1 << "." << extension;
    rename(s.str().c_str(), name.c_str());
    return 0;
  }

  return application.exec();
}
