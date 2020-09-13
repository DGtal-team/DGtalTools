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
 * @file shapeGenerator.cpp
 * @ingroup Tools
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr)
 * LIRIS (CNRS, UMR 5205),
 *
 * @date 2011/01/04
 *
 * DGtal shape generator
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <string>

#include "CLI11.hpp"

#include "DGtal/base/Common.h"

#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

#include "DGtal/io/writers/PGMWriter.h"
#include "DGtal/io/writers/RawWriter.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/boards/Board2D.h"


using namespace DGtal;


/**
 @page shapeGenerator shapeGenerator
 @brief  Generates shapes using DGtal library.
 


 @b Usage:  shapeGenerator [options] --shape <shapeName> --output <outputBasename>

 @b Allowed @b options @b are:

 @code
  -h,--help                             Print this help message and exit
  -l,--list                             List all available shapes
  -s,--shape TEXT                       Shape name
  -R,--radius FLOAT                     Radius of the shape
  -A,--axis1 FLOAT                      Half big axis of the shape (ellipse)
  -a,--axis2 FLOAT                      Half small axis of the shape (ellipse)
  -r,--smallradius FLOAT=5              Small radius of the shape (default 5)
  -v,--varsmallradius FLOAT=5           Variable small radius of the shape (default 5)
  -k UINT=3                             Number of branches or corners the shape (default 3)
  --phi FLOAT=0                         Phase of the shape (in radian, default 0.0)
  -w,--width FLOAT=10                   Width of the shape (default 10.0)
  -p,--power FLOAT=2                    Power of the metric (default 2.0)
  -o,--output TEXT                      Basename of the output file
  --signature                           Display to the standard output the signature (normal, curvature) at each point of the specified shape contour (middle point of each contour linel)
  -f,--format TEXT                      Output format:
                                            Bitmap {pgm, raw}
                                            Vector {svg} (+ {png,pdf} if libCairo installed) (default pgm)
 @endcode
 You can list the potential shapes:
 @code
 $ contourGenerator --list
 2D Shapes:
	ball	Ball for the Euclidean metric.
		Required parameter(s): --radius [-R]
	square	square (no signature).
		Required parameter(s): --width [-w]
	lpball	Ball for the l_power metric (no signature).
		Required parameter(s): --radius [-R], --power [-p]
	flower	Flower with k petals with radius ranging from R+/-v.
		Required parameter(s): --radius [-R], --varsmallradius [-v], --k [-k], --phi
	ngon	Regular k-gon.
		Required parameter(s): --radius [-R], --k [-k], --phi
	accflower	Accelerated Flower with k petals.
		Required parameter(s): --radius [-R], --varsmallradius [-v], --k [-k], --phi
	ellipse	Ellipse.
		Required parameter(s): --axis1 [-A], --axis2 [-a], --phi

 @endcode

 @b Example:
 @code
 # generate an accflower shape with 6 petals of maximal radius 40 and small radius 10:
 shapeGenerator -s accflower -R 40 -v 10 -k 6   -f pgm  -o test2
 @endcode

 You should obtain such a resulting image:
 @image html resShapeGenerator.png "resulting visualisation of generated shape."


 @see
 @ref shapeGenerator.cpp
 @ref contourGenerator
 */

/**
 * Global vectors to describe the available shapes and their
 * parameters.
 */
std::vector<std::string> shapes2D;
std::vector<std::string> shapesDesc;
std::vector<std::string> shapesParam1;
std::vector<std::string> shapesParam2;
std::vector<std::string> shapesParam3;
std::vector<std::string> shapesParam4;


/**
 * Create the static list of shapes.
 *
 */
void createList()
{
  shapes2D.push_back("ball");
  shapesDesc.push_back("Ball for the Euclidean metric.");
  shapesParam1.push_back("--radius [-R]");
  shapesParam2.push_back("");
  shapesParam3.push_back("");
  shapesParam4.push_back("");

  shapes2D.push_back("square");
  shapesDesc.push_back("square (no signature).");
  shapesParam1.push_back("--width [-w]");
  shapesParam2.push_back("");
  shapesParam3.push_back("");
  shapesParam4.push_back("");

  shapes2D.push_back("lpball");
  shapesDesc.push_back("Ball for the l_power metric (no signature).");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--power [-p]");
  shapesParam3.push_back("");
  shapesParam4.push_back("");

  shapes2D.push_back("flower");
  shapesDesc.push_back("Flower with k petals.");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--varsmallradius [-v],");
  shapesParam3.push_back("--k [-k],");
  shapesParam4.push_back("--phi");

  shapes2D.push_back("ngon");
  shapesDesc.push_back("Regular k-gon.");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--k [-k],");
  shapesParam3.push_back("--phi");
  shapesParam4.push_back("");

  shapes2D.push_back("accflower");
  shapesDesc.push_back("Accelerated Flower with k petals.");
  shapesParam1.push_back("--radius [-R],");
  shapesParam2.push_back("--varsmallradius [-v],");
  shapesParam3.push_back("--k [-k],");
  shapesParam4.push_back("--phi");

  shapes2D.push_back("ellipse");
  shapesDesc.push_back("Ellipse.");
  shapesParam1.push_back("--axis1 [-A],");
  shapesParam2.push_back("--axis2 [-a],");
  shapesParam3.push_back("--phi");
  shapesParam4.push_back("");


}

/**
 * Display the shape list with parameters.
 *
 */
void displayList()
{
  trace.emphase()<<"2D Shapes:"<<std::endl;
  for(unsigned int i=0; i<shapes2D.size(); ++i)
    trace.info()<<"\t"<<shapes2D[i]<<"\t"
		<<shapesDesc[i]<<std::endl
		<<"\t\tRequired parameter(s): "
		<< shapesParam1[i]<<" "
		<< shapesParam2[i]<<" "
		<< shapesParam3[i]<<" "
		<< shapesParam4[i]<<std::endl;

}


/**
 * Check if a given shape is available. If not, we exist with an error.
 * If it is, we return the corresponding index in the global vectors.
 *
 * @param shapeName name of the shape to search.
 *
 * @return index of the shape in the shape vectors.
 */
unsigned int checkAndRetrunIndex(const std::string &shapeName)
{
  unsigned int pos=0;

  while ((pos < shapes2D.size()) && (shapes2D[pos] != shapeName))
    pos++;

  if (pos == shapes2D.size())
    {
      trace.error() << "The specified shape has not found.";
      trace.info()<<std::endl;
      exit(1);
    }

  return pos;
}


/**
 * Functor to export a given shape into an image file
 * (pgm,raw,pdf,svg,...) and to extract its signature.
 *
 * @tparam Set type of the input Set
 * @tparam Image type of the input Image.
 */
template <typename Set, typename Image>
struct Exporter
{

  /**
   * Export a given Set into an image file.
   *
   * @param aSet input set.
   * @param outputName output file name.
   * @param outputFormat output file format.
   *
   */
  static
  void save(const Set &aSet,
	    const std::string outputName,
	    const std::string outputFormat)
  {

    Image  image = ImageFromSet<Image>::template create<Set>(aSet, 255, true);

    if  (outputFormat == "pgm")
      PGMWriter<Image>::exportPGM(outputName+"."+outputFormat,image);
    else
      if (outputFormat == "raw")
	RawWriter<Image>::exportRaw8(outputName+"."+outputFormat,image);
      else
	if (outputFormat == "svg")
	  {
	    Board2D board;
	    board << aSet;
	    board.saveSVG((outputName+"."+outputFormat).c_str());
	  }
	else
#ifdef WITH_CAIRO
	  if (outputFormat == "pdf")
	    {
	      Board2D board;
	      board << aSet;
	      board.saveCairo((outputName+"."+outputFormat).c_str(), Board2D::CairoPDF);

	    }
	  else
	    if (outputFormat == "png")
	      {
		Board2D board;
		board << aSet;
		board.saveCairo((outputName+"."+outputFormat).c_str(), Board2D::CairoPNG);
	      }
	    else
#endif
	      {
		trace.error()<< "Output format: "<<outputFormat<< " not recognized."<<std::endl;
		exit(1);
	      }
  }






  /**
   * Compute and export (std::cout) the boundary of the set and export the signature (normal
   * vector, curvature) at each point of the 2D contour.
   *
   * @param aShape the shape
   * @param aSet input set corresponding to the shape
   * @param aDomain the domain used to construct the set.
   */
  template <typename Shape>
  static
  void exportSignature(const Shape & aShape, Set &aSet, const Z2i::Domain &aDomain)
  {
    trace.beginBlock("Extracting the boundary");
    Z2i::KSpace ks;
    bool space_ok = ks.init( aDomain.lowerBound(),aDomain.upperBound(), true );
    SurfelAdjacency<2> sAdj( true );

    ASSERT(space_ok);
    trace.info() << aSet << std::endl;
    trace.info() << ks
		 << ( space_ok ? " Successfully instantiated" : " Error" )
		 << std::endl;

    std::vector< std::vector< Z2i::Point >  >  vectContoursBdryPointels;
    Surfaces<Z2i::KSpace>::extractAllPointContours4C( vectContoursBdryPointels,
						      ks, aSet, sAdj );
    trace.endBlock();

    ///Export
    std::cout<<"## shapeGenerator signature export"<<std::endl;
    std::cout<<"## shape: "<<aShape<<std::endl;
    std::cout<<"## x\ty\tdx\tdy\tddx\tddy"<<std::endl;
    for(unsigned int i=0; i<vectContoursBdryPointels.size(); i++)
      for(unsigned int j=0 ; j< vectContoursBdryPointels.at(i).size() - 1; j++)
	{
	  Z2i::Space::Point point = (vectContoursBdryPointels.at(i).at(j)
				     + vectContoursBdryPointels.at(i).at(j+1));
	  Z2i::Space::RealPoint midpoint (point[0]/2.0,point[1]/2.0);

	  Z2i::Space::RealPoint xp,xpp;
	  double t = aShape.parameter(midpoint);
	  xp = aShape.xp( t );
	  xpp = aShape.xpp( t );

	  std::cout<< midpoint[0]<<"\t"<<midpoint[1]<<"\t"
		   << xp[0]<<"\t"<<xp[1]<<"\t"
		   << xpp[0]<<"\t"<<xpp[1]<<std::endl;

	}

  }
};

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam(std::string param)
{
  trace.error() <<" Parameter: "<<param<<" is required..";
  trace.info()<<std::endl;
  exit(1);
}

///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
  // parse command line CLI ----------------------------------------------
  CLI::App app;
  std::string shapeName;
  std::string outputName;
  std::string outputFormat {"pgm"};
  double radius;
  double power {2.0};
  double smallradius {5};
  double varsmallradius {5};
  unsigned int k {3};
  double phi {0.0};
  double width {10.0};
  double axis1, axis2;
  
  app.description("Generates shapes using DGtal library.\n Typical use example:\n \t shapeGenerator [options] --shape <shapeName> --output <outputBasename>\n");
  auto listOpt = app.add_flag("--list,-l","List all available shapes");
  auto shapeNameOpt = app.add_option("--shape,-s", shapeName, "Shape name");
  auto radiusOpt = app.add_option("--radius,-R", radius, "Radius of the shape" );
  auto axis1Opt = app.add_option("--axis1,-A", axis1, "Half big axis of the shape (ellipse)" );
  auto axis2Opt = app.add_option("--axis2,-a", axis2, "Half small axis of the shape (ellipse)" );
  auto smallradiusOpt = app.add_option("--smallradius,-r", smallradius, "Small radius of the shape (default 5)", true);
  auto varsmallradiusOpt = app.add_option("--varsmallradius,-v", varsmallradius, "Variable small radius of the shape (default 5)", true );
  auto kOpt = app.add_option("-k", k, "Number of branches or corners the shape (default 3)", true );
  auto phiOpt = app.add_option("--phi", phi, "Phase of the shape (in radian, default 0.0)", true );
  auto widthOpt = app.add_option("--width,-w", width, "Width of the shape (default 10.0)", true );
  auto powerOpt = app.add_option("--power,-p", power, "Power of the metric (default 2.0)", true );
  auto outputNameOpt = app.add_option("--output,-o", outputName, "Basename of the output file");
  auto signatureOpt = app.add_flag("--signature", "Display to the standard output the signature (normal, curvature) at each point of the specified shape contour (middle point of each contour linel)");
  app.add_option("--format,-f", outputFormat, "Output format:\n\t  Bitmap {pgm, raw}\n\t  Vector {svg} (+ {png,pdf} if libCairo installed) (default pgm)" );

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------  

  //List creation
  createList();

  if ( listOpt->count() > 0 )
  {
    displayList();
    return 0;
  }

  if(shapeNameOpt->count()==0) missingParam("--shape");
  if(outputNameOpt->count()==0) missingParam("--output");

  //We check that the shape is known
  unsigned int id = checkAndRetrunIndex(shapeName);
  typedef ImageContainerBySTLVector<Z2i::Domain,unsigned char> Image;

  if (id ==0)
    {
      if (radiusOpt->count()==0) missingParam("--radius");
      Ball2D<Z2i::Space> ball(Z2i::Point(0,0), radius);
      Z2i::Domain domain(ball.getLowerBound(), ball.getUpperBound());
      Z2i::DigitalSet aSet(domain);

      Shapes<Z2i::Domain>::euclideanShaper(aSet, ball);
      Exporter<Z2i::DigitalSet,Image>::save(aSet,outputName,outputFormat);

      if (signatureOpt->count()>0)
        Exporter<Z2i::DigitalSet,Image>::exportSignature(ball,aSet,domain);

      return 0;
    }
  else
    if (id ==1)
      {
      	//if (widthOpt->count()==0) missingParam("--width");
      	ImplicitHyperCube<Z2i::Space> object(Z2i::Point(0,0), width/2.0);
      	Z2i::Domain domain(object.getLowerBound(), object.getUpperBound());
      	Z2i::DigitalSet aSet(domain);

      	Shapes<Z2i::Domain>::euclideanShaper(aSet, object);
      	Exporter<Z2i::DigitalSet,Image>::save(aSet,outputName,outputFormat);

      	if (signatureOpt->count()>0)
      	  {
      	    trace.error()<< "No signature export for this shape.";
      	    trace.info()<<std::endl;
      	  }

      	return 0;
      }
    else
      if (id ==2)
      	{
      	  //if (powerOpt->count()==0) missingParam("--power");
          if (radiusOpt->count()==0) missingParam("--radius");
      	  ImplicitRoundedHyperCube<Z2i::Space> ball(Z2i::Point(0,0), radius, power);
      	  Z2i::Domain domain(ball.getLowerBound(), ball.getUpperBound());
      	  Z2i::DigitalSet aSet(domain);

      	  Shapes<Z2i::Domain>::euclideanShaper(aSet, ball);
      	  Exporter<Z2i::DigitalSet,Image>::save(aSet,outputName,outputFormat);

      	  if (signatureOpt->count()>0)
      	    {
      	      trace.error()<< "No signature export for this shape.";
      	      trace.info()<<std::endl;
      	    }

      	  return 0;
      	}
      else
	if (id ==3)
	  {
	    //if (varsmallradiusOpt->count()==0) missingParam("--varsmallradius");
      if (radiusOpt->count()==0) missingParam("--radius");
      //if (kOpt->count()==0) missingParam("--k");
      //if (phiOpt->count()==0) missingParam("--phi");
	    Flower2D<Z2i::Space> flower(Z2i::Point(0,0), radius, varsmallradius,k,phi);
	    Z2i::Domain domain(flower.getLowerBound(), flower.getUpperBound());
	    Z2i::DigitalSet aSet(domain);

	    Shapes<Z2i::Domain>::euclideanShaper(aSet, flower);
	    Exporter<Z2i::DigitalSet,Image>::save(aSet,outputName,outputFormat);

	    if (signatureOpt->count()>0)
	      Exporter<Z2i::DigitalSet,Image>::exportSignature(flower,aSet,domain);

	    return 0;
	  }
	else
	  if (id ==4)
	    {
	      if (radiusOpt->count()==0) missingParam("--radius");
        //if (kOpt->count()==0) missingParam("--k");
        //if (phiOpt->count()==0) missingParam("--phi");
	      NGon2D<Z2i::Space> object(Z2i::Point(0,0), radius,k,phi);
	      Z2i::Domain domain(object.getLowerBound(), object.getUpperBound());
	      Z2i::DigitalSet aSet(domain);

	      Shapes<Z2i::Domain>::euclideanShaper(aSet, object);
	      Exporter<Z2i::DigitalSet,Image>::save(aSet,outputName,outputFormat);

	      if (signatureOpt->count()>0)
          Exporter<Z2i::DigitalSet,Image>::exportSignature(object,aSet,domain);

	      return 0;
	    }
	  else
	    if (id ==5)
	      {
      		//if (varsmallradiusOpt->count()==0) missingParam("--varsmallradius");
          if (radiusOpt->count()==0) missingParam("--radius");
          //if (kOpt->count()==0) missingParam("--k");
          //if (phiOpt->count()==0) missingParam("--phi");
      		AccFlower2D<Z2i::Space> flower(Z2i::Point(0,0), radius, varsmallradius,k,phi);
      		Z2i::Domain domain(flower.getLowerBound(), flower.getUpperBound());
      		Z2i::DigitalSet aSet(domain);

      		Shapes<Z2i::Domain>::euclideanShaper(aSet, flower);
      		Exporter<Z2i::DigitalSet,Image>::save(aSet,outputName,outputFormat);

      		if (signatureOpt->count()>0)
      		  Exporter<Z2i::DigitalSet,Image>::exportSignature(flower,aSet,domain);

      		return 0;
	      }
	    else
	      //if (id ==6)
	      {
      		if (axis1Opt->count()==0) missingParam("--axis1");
          if (axis2Opt->count()==0) missingParam("--axis2");
          //if (phiOpt->count()==0) missingParam("--phi");
      		Ellipse2D<Z2i::Space> ell(Z2i::Point(0,0), axis1, axis2,phi);
      		Z2i::Domain domain(ell.getLowerBound(), ell.getUpperBound());
      		Z2i::DigitalSet aSet(domain);

      		Shapes<Z2i::Domain>::euclideanShaper(aSet, ell);
      		Exporter<Z2i::DigitalSet,Image>::save(aSet,outputName,outputFormat);

      		if (signatureOpt->count()>0)
      		  Exporter<Z2i::DigitalSet,Image>::exportSignature(ell,aSet,domain);

      		return 0;
	      }
}
