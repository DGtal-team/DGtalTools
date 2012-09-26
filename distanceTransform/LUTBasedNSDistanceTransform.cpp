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
 * @file displayContours.cpp
 * @ingroup Tools
 * @author Nicolas Normand (\c kerautre@loria.fr)
 * IRCCyN (CNRS, UMR 6597), University of Nantes, France
 *
 * @date 2012/09/04
 *
 * DGtal convert grey scales image to fremann contour. 
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "DGtal/base/Common.h"

#include "DGtal/helpers/StdDefs.h"

#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/topology/helpers/Surfaces.h"

//image
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/Color.h"

#ifdef WITH_MAGICK
#include "DGtal/io/readers/MagickReader.h"
#endif


//contour
#include "DGtal/geometry/curves/FreemanChain.h"

//processing
#include "DGtal/geometry/curves/ArithmeticalDSS.h"
#include "DGtal/geometry/curves/GreedyDecomposition.h"
#include "DGtal/geometry/curves/MaximalSegments.h"
#include "DGtal/geometry/curves/FP.h"

//boost
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

//STL
#include <vector>
#include <string>

//Path-based distances
#define WITH_PNG

#include "ImageReader.h"

#include "BaseDistanceDT.h"
#include "D4DistanceDT.h"
#include "D8DistanceDT.h"
#include "RatioNSDistanceDT.h"
#include "PeriodicNSDistanceDT.h"

#include "ImageWriter.h"

using namespace DGtal;

////////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

#define PBM_FILE_FORMAT	1
#define PNG_FILE_FORMAT	2

int main( int argc, char** argv )
{
    // parse command line ----------------------------------------------------//
    po::options_description general_opt("Options and arguments: ");
    general_opt.add_options()
	("help,h", "Display this message")
	("city-block,4", "Use the city block distance")
	("chessboard,8", "Use the chessboard distance")
	("sequence,s", po::value<std::string>(),
	 "One period of the sequence of neighborhoods given as a list of 1 "
	 "and 2 separated by \" \" or \",\". Space characters must be escaped "
	 "from the shell.")
	("ratio,r", po::value<std::string>(),
	 "Ratio of neighborhood 2 given as the rational number num/den "
	 "(with den >= num >= 0 and den > 0).")
	("center,c", "Center the distance transform (the default is an asymmetric "
	 "distance transform)")
	("input,i", po::value<std::string>(), "Read from file \"arg\" instead of "
	 "stdin.")
	("drawContourPoint", po::value<double>(), "<size> display contour points as disk of radius <size>")    
	("fillContour", "fill the contours with default color")
	("lineWidth", po::value<double>()->default_value(1.0), "Define the linewidth of the contour (SDP format)") 
	("withProcessing", po::value<std::string>(), "Processing (used only with --FreemanChain):\n\t DSS segmentation {DSS}\n\t  Maximal segments {MS}\n\t Faithful Polygon {FP}\n\t Minimum Length Polygon {MLP}")   
	("outputEPS", po::value<std::string>(), " <filename> specify eps format (default format output.eps)")
	("outputSVG", po::value<std::string>(), " <filename> specify svg format.")
	("outputFIG", po::value<std::string>(), " <filename> specify fig format.")
#ifdef WITH_CAIRO
	("outputPDF", po::value<std::string>(), "outputPDF <filename> specify pdf format. ")
	("outputPNG", po::value<std::string>(), "outputPNG <filename> specify png format.")
	("invertYaxis", " invertYaxis invert the Y axis for display contours (used only with --SDP)")
#endif
#ifdef WITH_MAGICK
	("backgroundImage", po::value<std::string>(), "backgroundImage <filename> <alpha> : display image as background with transparency alpha (defaut 1) (transparency works only if cairo is available)")
	("alphaBG", po::value<double>(), "alphaBG <value> 0-1.0 to display the background image in transparency (default 1.0)")
#endif
	("scale", po::value<double>(), "scale <value> 1: normal; >1 : larger ; <1 lower resolutions  )");
    //------------------------------------------------------------------------//

    bool parseOK = true;
    po::variables_map vm;
    try {
	po::store(po::parse_command_line(argc, argv, general_opt), vm);  
    }
    catch (const std::exception& ex) {
	parseOK = false;
	trace.info() << "Error checking program options: " << ex.what() << endl;
    }

    po::notify(vm);    
    if(!parseOK || vm.count("help") || argc <= 1 ||
       //(vm.count("chessboard") + vm.count("city-block") + vm.count("ratio") + vm.count("sequence") != 1)
       (vm.count("chessboard") + vm.count("city-block") != 1)
       ) {
	trace.info() <<
	    "Compute the 2D translated neighborhood-sequence "
	    "distance transform of a binary image" << std::endl <<
	    "Basic usage: " << std::endl
	<< "\tLUTBasedNSDistanceTransform [-f filename] [-c] (-4|-8|-r <num/den>|-s <sequence>) [-t (pgm|png)]" <<std::endl
	<< general_opt << "\n";
	return 0;
    }

    FILE *input = stdin;

    Board2D aBoard;
    double scale = 1.0;
    aBoard.setUnit (0.05*scale, LibBoard::Board::UCentimeter);

    int inputFormat = 0;
    //sourceOfDT source = undefined;
    int translateFlag = 0;
    char *myName = argv[0];
    char *outputFormat = NULL;
    bool lineBuffered = false;
    BaseDistance *dist = NULL;

    if (vm.count("city-block")){
	dist = new D4Distance();
    }
    else if (vm.count("chessboard")){
	dist = new D8Distance();
    }

    ImageConsumer<GrayscalePixelType> *output = createImageWriter("-", outputFormat, lineBuffered);
    if (output == NULL) {
	fprintf(stderr, "Unable to create image output stream (wrong format?)\n");
    }

    if (translateFlag) {
	output = dist->newDistanceTransformUntranslator(output);
    }
    BaseDistanceTransform *dt = dist->newTranslatedDistanceTransform(output);

#ifdef WITH_NETPBM
    if (inputFormat == 0) {
	char c = fgetc(input);
	ungetc(c, input);

	if (c == 'P') {
	    inputFormat = PBM_FILE_FORMAT;
	}
    }

    if (inputFormat == PBM_FILE_FORMAT) {
	PBMImageReader producer(dt, input);
	dt = NULL;

	while (!feof(input)) {
	    producer.produceAllRows();

	    int c;
	    while ((c = fgetc(input)) == '\n');
	    if (c != EOF)
		ungetc(c, input);
	}
    }
#endif

#ifdef WITH_PNG
    unsigned char signature[8];
    int readBytes = 0;

    if (inputFormat == 0) {
	readBytes = fread(signature, 1, 8, input);
    }

    if (inputFormat == PNG_FILE_FORMAT ||
	inputFormat == 0 && readBytes == 8 && png_check_sig(signature, 8)) {
	
	inputFormat = PNG_FILE_FORMAT;
	PNGImageReader producer(dt, input);
	dt = NULL;

	do {
	    // Assumes following images, if any, are png
	    producer.produceAllRows(8);
	}
	while (fread(signature, 1, 8, input) == 8);
    }
#endif

    if (inputFormat == 0) {
	fprintf(stderr, "Input image format not recognized\n");
    }

    if (dt != NULL) {
	delete dt;
    }

    fclose(input);
    return 0;
}
