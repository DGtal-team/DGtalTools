/**
 * This Source Code Form is subject to the terms of the Mozilla Public License,
 * v. 2.0. If a copy of the MPL was not distributed with this file, You can
 * obtain one at http://mozilla.org/MPL/2.0/.
**/

/**
 * @file criticalKernelsThinning3D.cpp
 * @author Pablo Hernandez-Cerdan (\c pablo.hernandez.cerdan@outlook.com )
 * Institute of Fundamental Sciences, Massey University, New Zealand
 *
 * @date 2017/11/04
 *
 * Apply a thinning using the Critical Kernel Framework.
 *
 * This file is part of the DGtal library.
 */

/**
 @page criticalKernelsThinning3D criticalKernelsThinning3D

 @brief Applies an criticalKernels thinning algorithm of a 3d image file (vol,longvol,pgm3d...) with 3D viewer.

 @b Usage: criticalKernelsThinning3D [options] --input <3dImageFileName>  {vol,longvol,pgm3d...}


 @b Allowed @b options @b are :
 @code
  -h [ --help ]                    display this message.
  -i [ --input ] arg               Input vol file.
  -s [ --skel] arg                 Skeletonization: only keep certain voxels during thinning.
                                   Options: ulti,end, 1isthmus, isthmusulti: delete all voxels except those that change topology.
                                   end: keep voxels with only one neighbor.
                                   1isthmus: keep voxels that are one-isthmus (using LookUpTables) [faster]
                                   isthmus: keep voxels that are one-isthmus or two-isthmus (using LookUpTables) [faster]

  -s [ --skel ] arg                type of skeletonization
  -c [ --select] arg               Select: order in which select voxels in the process.
                                   Options: dmax, first, random
                                   dmax: Use distance map, selecting voxel with max value.
                                   first: Select first pixel (lexicographical order)
                                   random: Select voxel at random.
  -f [ --foreground ] arg (=black) foreground color in binary image
  -m [ --thresholdMin ] arg (=0)   threshold min (excluded) to define binary 
                                   shape
  -M [ --thresholdMax ] arg (=255) threshold max (included) to define binary 
                                   shape
  -p [ --persistence ] arg (=0)    persistence value, implies use of 
                                   persistence algorithm if p>=1
  --profile                        profile algorithm
  -v [ --verbose ]                 verbose output
  -o [ --exportImage ] arg         Export the resulting set of points to a 
                                   image compatible with GenericWriter.
  -e [ --exportSDP ] arg           Export the resulting set of points in a 
                                   simple (sequence of discrete point (sdp)).
  -t [ --visualize ]               visualize result in viewer


 @endcode

 @b Example:

@code
 $  criticalKernelsThinning3D --input ${DGtal}/examples/samples/Al.100.vol --select dmax --skel 1isthmus --persistence 1 -t
 @endcode

 You should obtain such a result:
 @image html resCriticalKernelsThinning3D_select-dmax_skel-1isthmus_persistence-1.png "Resulting visualization."

 @see
 @ref criticalKernelsThinning3D.cpp
 */

#include <iostream>
#include <chrono>
#include <unordered_map>

#include <DGtal/io/viewers/Viewer3D.h>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/SimpleThresholdForegroundPredicate.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"

#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/topology/CubicalComplex.h>
#include <DGtal/topology/CubicalComplexFunctions.h>

#include <DGtal/topology/VoxelComplex.h>
#include <DGtal/topology/VoxelComplexFunctions.h>
#include "DGtal/topology/NeighborhoodConfigurations.h"
#include "DGtal/topology/tables/NeighborhoodTables.h"

// boost::program_options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
// Distance Map
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/images/SimpleThresholdForegroundPredicate.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/VoronoiMap.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
using namespace DGtal;
using namespace std;
using namespace DGtal::Z3i;
namespace po = boost::program_options;

int main(int argc, char* const argv[]){

  /*-------------- Parse command line -----------------------------*/
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
    ( "help,h", "display this message." )
    ( "input,i", po::value<string>()->required(), "Input vol file." )
    ( "skel,s",  po::value<string>()->required(), "type of skeletonization" )
    ( "select,c",  po::value<string>()->required(), "select method for skeletonization" )
    ( "foreground,f",  po::value<string>()->default_value("black"), "foreground color in binary image" )
    ( "thresholdMin,m",  po::value<int>()->default_value(0), "threshold min (excluded) to define binary shape" )
    ( "thresholdMax,M",  po::value<int>()->default_value(255), "threshold max (included) to define binary shape" )
    ( "persistence,p",  po::value<int>()->default_value(0), "persistence value, implies use of persistence algorithm if p>=1" )
    ( "profile",  po::bool_switch()->default_value(false), "profile algorithm" )
    ( "verbose,v",  po::bool_switch()->default_value(false), "verbose output" )
    ( "exportImage,o", po::value<std::string>(), "Export the resulting set of points to a image compatible with GenericWriter.")
    ( "exportSDP,e", po::value<std::string>(), "Export the resulting set of points in a simple (sequence of discrete point (sdp)).")
    ( "visualize,t",  po::bool_switch()->default_value(false), "visualize result in viewer" );
  bool parseOK=true;
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
    po::notify ( vm );
  } catch(const exception& ex) {
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }

  if (!parseOK || vm.count ( "help" ) || !vm.count("input"))
  {
    trace.info() <<
    "Compute the thinning of a volume using the CriticalKernels framework"<< std::endl
    << std::endl << "Basic usage: "<< std::endl
    << "criticalKernelsThinning3D --input <volFileName> --skel <ulti,end, 1isthmus, isthmus> --select "
    " [ -f <white,black> -m <minlevel> -M <maxlevel> -v ] "
    " [--persistence <value> ]" << std::endl
    << "options for --skel {ulti end 1isthmus isthmus}" << std::endl
    << "options for --select = {dmax random first}" << std::endl
    << general_opt << "\n"
    << " Example: \n"
      << "criticalKernelsThinning3D --input ${DGtal}/examples/samples/Al.100.vol --select dmax --skel 1isthmus --persistence 1 --visualize --verbose --outputImage ./Al100_dmax_1isthmus_p1.vol \n";
    return 0;
  }
  //Parse options
  string filename = vm["input"].as<string>();
  bool verbose = vm["verbose"].as<bool>();
  bool visualize = vm["visualize"].as<bool>();
  bool profile = vm["profile"].as<bool>();
  int thresholdMin = vm["thresholdMin"].as<int>();
  int thresholdMax = vm["thresholdMax"].as<int>();
  int persistence = vm["persistence"].as<int>();
  if (vm.count("persistence") && persistence < 0 )
    throw po::validation_error(po::validation_error::invalid_option_value, "persistence");
  string foreground = vm["foreground"].as<string>();
  if (vm.count("foreground") && (!(foreground == "white" || foreground == "black")))
    throw po::validation_error(po::validation_error::invalid_option_value, "foreground");

  string sk_string = vm["skel"].as<string>();
  if (vm.count("skel") &&
     (!( sk_string == "ulti" || sk_string == "end" ||
         sk_string == "isthmus" || sk_string == "1isthmus"))
     )
     throw po::validation_error(po::validation_error::invalid_option_value, "skel");
  string select_string = vm["select"].as<string>();
  if (vm.count("select") &&
     (!( select_string == "random" || select_string == "dmax" ||
         select_string == "first" ))
     )
     throw po::validation_error(po::validation_error::invalid_option_value, "select");
  /*-------------- End of parse -----------------------------*/

  if(verbose){
    std::cout << "Skel: " << sk_string << std::endl;
    std::cout << "Select: " << select_string << std::endl;
    std::cout << "Persistence: " << persistence << std::endl;
    std::cout << "Input: " << filename << std::endl;
  }
  trace.beginBlock("Reading input");
  using Domain = Z3i::Domain ;
  using Image = ImageSelector < Z3i::Domain, unsigned char>::Type ;
  Image image = GenericReader<Image>::import(filename);
  trace.endBlock();

  DigitalSet image_set (image.domain());
  SetFromImage<Z3i::DigitalSet>::append<Image>(
      image_set, image,
      thresholdMin, thresholdMax);


  // Create a VoxelComplex from the set

  using DigitalTopology = DT26_6;
  using DigitalSet =
    DGtal::DigitalSetByAssociativeContainer<Domain ,
      std::unordered_set< typename Domain::Point> >;
  using Object =
    DGtal::Object<DigitalTopology, DigitalSet>;
  using Complex =
    DGtal::VoxelComplex<KSpace, Object>;

  auto & sk = sk_string;
  KSpace ks;
  KSpace::Point d1( KSpace::Point::diagonal( 1 ) );
  ks.init(image.domain().lowerBound() - d1 ,
      image.domain().upperBound() + d1 , true);

  DigitalTopology::ForegroundAdjacency adjF;
  DigitalTopology::BackgroundAdjacency adjB;
  DigitalTopology topo(adjF, adjB, DGtal::DigitalTopologyProperties::JORDAN_DT);
  Object obj(topo,image_set);

  trace.beginBlock("construct with table");
  Complex vc(ks);
  vc.construct(obj.pointSet(), functions::loadTable(simplicity::tableSimple26_6 ));
  trace.endBlock();
  trace.beginBlock("load isthmus table");
  boost::dynamic_bitset<> isthmus_table;
  if (sk == "isthmus")
    isthmus_table = *functions::loadTable(isthmusicity::tableIsthmus);
  else if (sk == "1isthmus")
    isthmus_table = *functions::loadTable(isthmusicity::tableOneIsthmus);
  auto pointMap = *functions::mapZeroPointNeighborhoodToConfigurationMask<Point>();

  trace.endBlock();
  using namespace DGtal::functions ;
  // SKEL FUNCTION:
  std::function< bool(const Complex&, const Cell&) > Skel ;
  if (sk == "ulti") Skel = skelUltimate<Complex>;
  else if (sk == "end") Skel = skelEnd<Complex>;
  else if (sk == "isthmus" || sk == "1isthmus")
      Skel = [&isthmus_table, &pointMap](const Complex & fc,
               const Complex::Cell & c){
        return skelWithTable(isthmus_table, pointMap, fc, c);
      };
  else throw std::runtime_error("Invalid skel string");
  auto start = std::chrono::system_clock::now();

  // SELECT FUNCTION
  /*
   * Calculate distance map even if not requested:
   */
  trace.beginBlock("Create Distance Map");
  using Predicate = Z3i::DigitalSet;
  using L3Metric = ExactPredicateLpSeparableMetric<Z3i::Space, 3>;
  using DT       = DistanceTransformation<Z3i::Space, Predicate, L3Metric>;
  L3Metric l3;
  DT dt(obj.domain(),obj.pointSet(), l3);
  trace.endBlock();

  std::function< std::pair<typename Complex::Cell, typename Complex::Data>(const Complex::Clique&) > Select ;
  auto & sel = select_string;
  if (sel == "random") Select = selectRandom<Complex>;
  else if (sel == "first") Select = selectFirst<Complex>;
  else if (sel == "dmax"){
    Select =
      [&dt](const Complex::Clique & clique){
        return selectMaxValue<DT, Complex>(dt,clique);
      };
  } else throw std::runtime_error("Invalid skel string");

  trace.beginBlock("Thinning");
  Complex vc_new(ks);
  if (persistence == 0)
    vc_new = asymetricThinningScheme< Complex >(
        vc, Select,  Skel, verbose);
  else
    vc_new = persistenceAsymetricThinningScheme< Complex >(
        vc, Select,  Skel, persistence, verbose);
  trace.endBlock();

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds> (end - start) ;
  if (profile) std::cout <<"Time elapsed: " << elapsed.count() << std::endl;

  const auto & thin_set = vc_new.objectSet();
  const auto & all_set = obj.pointSet();

  if (vm.count("exportSDP"))
  {
    std::ofstream out;
    out.open(vm["exportSDP"].as<std::string>().c_str());
    for (auto &p : thin_set)
    {
      out << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
  }

  if (vm.count("exportImage"))
  {
    auto outputFilename = vm["exportImage"].as<std::string>();
    if(verbose)
      std::cout << "outputFilename" << outputFilename << std::endl;

    unsigned int foreground_value = 255;
    auto thin_image = ImageFromSet<Image>::create(thin_set, foreground_value);
    thin_image >> outputFilename;
    // // ITK output
    // typedef itk::ImageFileWriter<Image::ITKImage> ITKImageWriter;
    // typename ITKImageWriter::Pointer writer = ITKImageWriter::New();
    // writer->SetFileName(outputFilename.c_str());
    // writer->SetInput(thin_image.getITKImagePointer());
    // writer->Update();
  }

  if(visualize)
  {
    int argc(1);
    char** argv(nullptr);
    QApplication app(argc, argv);
    Viewer3D<> viewer;
    viewer.setWindowTitle("criticalKernelsThinning3D");
    viewer.show();

    viewer.setFillColor(Color(255, 255, 255, 255));
    viewer << thin_set;

    // All kspace voxels
    viewer.setFillColor(Color(40, 200, 55, 10));
    viewer << all_set;

    viewer << Viewer3D<>::updateDisplay;

    app.exec();
  }
}
