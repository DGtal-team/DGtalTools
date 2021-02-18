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
 * @file volEuclideanDT.cpp
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2021/01/25
 *
 *
 * This file is part of the DGtal library.
 */

#include <iostream>
#include <algorithm>
#include <functional>
#include <DGtal/base/Common.h>
#include <DGtal/io/readers/VolReader.h>
#include <DGtal/io/writers/VolWriter.h>
#include <DGtal/io/writers/LongvolWriter.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/geometry/volumes/distance/DistanceTransformation.h>
#include <DGtal/images/SimpleThresholdForegroundPredicate.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include "DGtal/geometry/volumes/distance/PowerMap.h"
#include "DGtal/geometry/volumes/distance/ReducedMedialAxis.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpPowerSeparableMetric.h"

#include <CLI11.hpp>


using namespace std;
using namespace DGtal;
using namespace Z3i;


/**
 @page volEuclideanDT volEuclideanDT
 
 @brief Tool to compute and output Euclidean distance related quantitites (the distance transform, the Voronoi map, the RDMA...).


 @b Usage: 	./volumetric/volEuclideanDT [OPTIONS] 1 [2]


 @b Allowed @b options @b are : 
 @code
 
 @endcode

 @b Example: 

 @see
 @ref volEuclideanDT.cpp

 */

template<typename T>
struct functorCast{
  
  unsigned char operator()(const T &v) {return v%256;}
};


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
  CLI::App app;
  
  app.description("Brutally sub sample a vol file (division by 2 in each direction).\n Basic usage: \n \tvolSubSample --input <volFileName> --o <volOutputFileName> ");

  std::string inputFileName;
  app.add_option("-i,--input,1", inputFileName, "Input vol file." )->required()->check(CLI::ExistingFile);
  std::string outputFileName;
  app.add_option("-o,--output,2", outputFileName, "Output filename.",true);
  
  std::string mode="edt";
  app.add_option("-m,--mode", mode, "Export mode for the distance: {edt (remappred distances to the [0:255] range),sedt (exact square of distances as longvol),voronoi (Voronoi map using hash value per cell), rdma (centre of maximal balls)} (default:edt)", true)-> check(CLI::IsMember({"edt", "sedt", "voronoi", "rdma"}));
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);

  
  trace.beginBlock("Loading file");
  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  Image;
  Image image = VolReader< Image >::importVol ( inputFileName );
  trace.info()<<image<<std::endl;
  trace.endBlock();

  
  
  
  trace.beginBlock("DT");
  typedef functors::SimpleThresholdForegroundPredicate<Image> Predicate;
  Predicate aPredicate(image,0);
  DistanceTransformation<Z3i::Space, Predicate, Z3i::L2Metric> dt(image.domain(), aPredicate, Z3i::l2Metric);
  trace.info()<<dt<<std::endl;
  trace.endBlock();
  
  double valMax = 0.0;
  double valMin = std::numeric_limits<double>::max();
  for(auto v: dt.constRange())
  {
    valMax=std::max(valMax,v);
    valMin=std::min(valMin,v);
  }
  
  if (mode == "edt")
  {
    Image output(image.domain());
    for(auto &voxel: image.domain())
    {
      auto val = dt(voxel);
      unsigned char valC = static_cast<unsigned char>( 255*(val+valMin)/(valMax-valMin) );
      output.setValue(voxel, valC);
    }
    VolWriter<Image>::exportVol(outputFileName, output);
  }
  else
    if (mode == "sedt")
    {
      ImageContainerBySTLVector<Z3i::Domain, uint64_t>   output(image.domain());
      for(auto &voxel: image.domain())
      {
        auto val = l2Metric.rawDistance(voxel, dt.getVoronoiVector(voxel));
        output.setValue(voxel, val);
      }
      size_t lastindex = outputFileName.find_last_of(".");
      std::string rawname = outputFileName.substr(0, lastindex);
      LongvolWriter< ImageContainerBySTLVector<Z3i::Domain, uint64_t>>::exportLongvol(rawname+".longvol", output);
    }
  else
    if (mode == "voronoi")
    {
      auto myhash=[](const Z3i::Point &p){return (((p[0])*43 + p[1]*1777)*123 + p[2]) % 256;};
      Image output(image.domain());
      for(auto &voxel: image.domain())
      {
        auto site = dt.getVoronoiVector(voxel);
        unsigned char v=0;
        if (site == voxel)
          v = 0;
        else
          v = myhash(site);
        output.setValue(voxel, v);
      }
      VolWriter<Image>::exportVol(outputFileName, output);
    }
 else
   if (mode == "rdma")
   {
     typedef ImageContainerBySTLVector<Z3i::Domain, uint64_t> ImageLong;
     ImageLong rawDT(image.domain());
     for(auto &voxel: image.domain())
     {
       auto val = l2Metric.rawDistance(voxel, dt.getVoronoiVector(voxel));
       rawDT.setValue(voxel, val);
     }
     PowerMap<ImageLong, Z3i::L2PowerMetric> powermap(rawDT.domain(), rawDT, l2PowerMetric);
     ReducedMedialAxis<PowerMap<ImageLong, Z3i::L2PowerMetric> >::Type  rdma = ReducedMedialAxis< PowerMap<ImageLong, Z3i::L2PowerMetric> >::getReducedMedialAxisFromPowerMap(powermap);
     auto cpt2=0;
     for(auto v: dt.constRange())
       if (v!=0) cpt2++;
     auto cpt=0;
    
     Image out(image.domain());
     std::vector<std::pair<Point,double>> ma;
     for(auto &voxel: rdma.domain())
     {
       out.setValue(voxel, rdma(voxel) % 256);
       if (rdma(voxel) != 0)
       ma.push_back( std::pair<Point,double>(voxel ,std::sqrt(rdma(voxel))) );
     }
     VolWriter<Image>::exportVol(outputFileName, out);
     
     std::ofstream ofs ("trace.json", std::ofstream::out);
     ofs.precision(8);
     ofs<<std::fixed;
     ofs <<"scale(radius) union() {"<<std::endl;
     auto maxMA=100;
     cpt=0;
     std::sort(ma.begin(), ma.end(), [=](std::pair<Point, double>& a, std::pair<Point, double>& b)
               {return a.second > b.second;});
     
     //Extent
     double scale= (rdma.domain().upperBound()-rdma.domain().lowerBound()).normInfinity();
     
     for(cpt=0; cpt< ma.size(); cpt++)
     {
       Point voxel;
       double d;
       std::tie(voxel,d) = ma[ cpt ];
       RealPoint p(voxel);
       if (d < 2.0) {std::cout<<"+"<<std::flush; continue;}
       p = (p-rdma.domain().lowerBound()) / (scale);
       ofs<< "translate(["<<p[0]<<","<<p[1]<<","<<p[2]<<"]) sphere("<<d / scale<<");"<<std::endl;
     }
     ofs<<"}"<<std::endl;
     ofs.close();
     trace.info()<<"Found "<<cpt<<" spheres out of "<<cpt2<<std::endl;
     
   }
  

  return 0;
}
