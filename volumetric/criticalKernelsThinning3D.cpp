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
 Positionals:
   1 TEXT:FILE REQUIRED                  Input vol file.

 Options:
   -h,--help                             Print this help message and exit
   -i,--input TEXT:FILE REQUIRED         Input vol file.
   -s,--skel TEXT:{ulti,end,isthmus,1isthmus}=1isthmus
                                         Type of skeletonization. Options: 1isthmus, isthmus, end, ulti.
                                         Options: ulti,end, 1isthmus, isthmusulti: delete all voxels except those that change topology.
                                         end: keep voxels with only one neighbor.
                                         1isthmus: keep voxels that are one-isthmus (using LookUpTables) [faster]
                                         isthmus: keep voxels that are one-isthmus or two-isthmus (using LookUpTables) [faster]
   -c,--select TEXT:{random,dmax,first}=dmax
                                         Select the ordering for skeletonization. Options: dmax, random, first
                                          Options: dmax, first, random
                                          dmax: Use distance map, selecting voxel with max value.
                                          first: Select first pixel (lexicographical order)
                                          random: Select voxel at random.
   -f,--foreground TEXT:{white,black}=black
                                         Foreground color in binary image
   -m,--thresholdMin INT=0               Threshold min (excluded) to define binary shape
   -M,--thresholdMax INT=255             Threshold max (included) to define binary shape
   -p,--persistence INT:POSITIVE=0       Persistence value, implies use of persistence algorithm if p>=1
   --profile                             Profile algorithm
   -v,--verbose                          Verbose output
   -o,--exportImage TEXT                 Export the resulting set of points to a image compatible with GenericWriter.
   -e,--exportSDP TEXT                   Export the resulting set of points in a simple (sequence of discrete point (sdp)).
   -t,--visualize                        Visualize result in viewer



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

// Distance Map
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/VoronoiMap.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"

#include "CLI11.hpp"


using namespace DGtal;
using namespace std;
using namespace DGtal::Z3i;

int main(int argc, char* const argv[]){
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFileName;
  std::string sk_string {"1isthmus"};
  string select_string {"dmax"};
  string foreground {"black"};
  string outputFilenameImg;
  string outputFilenameSDP;

  int thresholdMin {0};
  int thresholdMax {255};
  int persistence {0};
  bool profile {false};
  bool verbose  {false};
  bool visualize {false};
  
  app.description("Compute the thinning of a volume using the CriticalKernels framework\nBasic usage: criticalKernelsThinning3D --input <volFileName> --skel <ulti,end, 1isthmus, isthmus> --select [ -f <white,black> -m <minlevel> -M <maxlevel> -v ] [--persistence <value> ] --persistence <value> ] \n options for --skel {ulti end 1isthmus isthmus} \n options for --select = {dmax random first} \n Example: \n criticalKernelsThinning3D --input ${DGtal}/examples/samples/Al.100.vol --select dmax --skel 1isthmus --persistence 1 --visualize --verbose --exportImage ./Al100_dmax_1isthmus_p1.vol \n");
  app.add_option("-i,--input,1", inputFileName, "Input vol file." )
  ->required()
  ->check(CLI::ExistingFile);

  app.add_option("--skel,-s", sk_string,"Type of skeletonization. Options: 1isthmus, isthmus, end, ulti.", true )
   -> check(CLI::IsMember({"ulti", "end","isthmus", "1isthmus"}));
  app.add_option("--select,-c", select_string, "Select the ordering for skeletonization. Options: dmax, random, first", true)
   -> check(CLI::IsMember({"random", "dmax", "first"}));
  app.add_option("--foreground,-f",foreground, "Foreground color in binary image", true )
   -> check(CLI::IsMember({"white", "black"}));
  app.add_option("--thresholdMin,-m", thresholdMin, "Threshold min (excluded) to define binary shape", true );
  app.add_option("--thresholdMax,-M", thresholdMax, "Threshold max (included) to define binary shape", true );
  app.add_option("--persistence,-p",persistence,"Persistence value, implies use of persistence algorithm if p>=1", true )
  ->check(CLI::PositiveNumber);
  app.add_flag("--profile", profile, "Profile algorithm");
  app.add_flag("--verbose,-v",verbose, "Verbose output");
  app.add_option("--exportImage,-o",outputFilenameImg, "Export the resulting set of points to a image compatible with GenericWriter.");
  app.add_option("--exportSDP,-e",outputFilenameSDP, "Export the resulting set of points in a simple (sequence of discrete point (sdp))." );
  app.add_flag("--visualize,-t", visualize, "Visualize result in viewer");
    
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  
  if(verbose){
    std::cout << "Skel: " << sk_string << std::endl;
    std::cout << "Select: " << select_string << std::endl;
    std::cout << "Persistence: " << persistence << std::endl;
    std::cout << "Input: " << inputFileName << std::endl;
  }
  trace.beginBlock("Reading input");
  using Domain = Z3i::Domain ;
  using Image = ImageSelector < Z3i::Domain, unsigned char>::Type ;
  Image image = GenericReader<Image>::import(inputFileName);
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
  using Complex = DGtal::VoxelComplex<KSpace>;

  auto & sk = sk_string;
  KSpace ks;
  KSpace::Point d1( KSpace::Point::diagonal( 1 ) );
  ks.init(image.domain().lowerBound() - d1 ,
      image.domain().upperBound() + d1 , true);

  trace.beginBlock("construct with table");
  Complex vc(ks);
  vc.construct(image_set, functions::loadTable(simplicity::tableSimple26_6 ));
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
  using L3Metric = ExactPredicateLpSeparableMetric<Z3i::Space, 3>;
  using DT       = DistanceTransformation<Z3i::Space, DigitalSet, L3Metric>;
  L3Metric l3;
  DT dt(image.domain(), image_set, l3);
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


  DigitalSet thin_set(image.domain());
  vc_new.dumpVoxels(thin_set);
  const auto & all_set = image_set;

  if (outputFilenameSDP != "")
  {
    std::ofstream out;
    out.open(outputFilenameSDP.c_str());
    for (auto &p : thin_set)
    {
      out << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
  }

  if (outputFilenameImg != "")
  {
    if(verbose)
      std::cout << "outputFilename" << outputFilenameImg << std::endl;

    unsigned int foreground_value = 255;
    auto thin_image = ImageFromSet<Image>::create(thin_set, foreground_value);
    thin_image >> outputFilenameImg;
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
