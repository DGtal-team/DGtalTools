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
 * @file LUTBasedNSDistanceTransform.cpp
 * @ingroup Tools
 * @author Nicolas Normand (\c Nicolas.Normand@polytech.univ-nantes.fr)
 * LUNAM Université, Université de Nantes, IRCCyN UMR CNRS 6597
 *
 * @date 2012/09/28
 *
 * LUTBasedNSDistanceTransform computes the 2D translated neighborhood-sequence
 * distance transform of a binary image. It reads the input images from its
 * standard input and writes the result to its standard output.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////


#if !defined(WITH_NETPBM) && !defined(WITH_PNG)
#pragma message ( "Neither libpng nor libnetpbm is available, image i/o is disabled," )
#pragma message ( "use cmake with option -DWITH_PNG or -DWITH_NETPBM" )
#endif

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
#include "DGtal/geometry/curves/FP.h"

//boost
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

//STL
#include <vector>
#include <string>

//Path-based distances
#include "ImageReader.h"

#include "NeighborhoodSequenceDistance.h"
#include <boost/tokenizer.hpp>

#include "ImageWriter.h"

using namespace DGtal;

////////////////////////////////////////////////////////////////////////////////

std::vector<int> parseSequence(std::string string) {
    std::vector<int> args;
    boost::tokenizer<> tok(string);
    for (boost::tokenizer<>::iterator beg = tok.begin(); beg != tok.end(); ++beg){
	//args.push_back(boost::lexical_cast<short>(*beg));
	args.push_back(boost::lexical_cast<short>(*beg));
    }
    return args;
}

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
	("output,o", po::value<std::string>(), "Output file name, optionally "
	 "prefixed with the file format and ':'")
	("outputFormat,t", po::value<std::string>(), "Output file format")
	("lineBuffered,l", "Flush output after each produced row.")
	("input,i", po::value<std::string>(), "Read from file \"arg\" instead of "
	 "stdin.");
    //------------------------------------------------------------------------//

    bool parseOK = true;
    po::variables_map vm;
    try {
	po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }
    catch (const std::exception& ex) {
	parseOK = false;
	trace.info() << "Error checking program options: " << ex.what() << std::endl;
    }

    po::notify(vm);
    if(!parseOK || vm.count("help") || argc <= 1 ||
       (vm.count("chessboard") +
	vm.count("city-block") +
	vm.count("ratio") +
	vm.count("sequence") != 1))
    {
	trace.info() <<
	    "Compute the 2D translated neighborhood-sequence "
	    "distance transform of a binary image" << std::endl <<
	    "Basic usage: " << std::endl
	<< "\tLUTBasedNSDistanceTransform [-i filename] [-c] (-4|-8|-r <num/den>|-s <sequence>) [-t (pgm|png)]" <<std::endl
	<< general_opt << "\n";
	return 0;
    }

    //sourceOfDT source = undefined;
    NeighborhoodSequenceDistance *dist = NULL;

    // Distance selection ----------------------------------------------------//
    if (vm.count("city-block")) {
	dist = NeighborhoodSequenceDistance::newD4Instance();
    }
    else if (vm.count("chessboard")) {
	dist = NeighborhoodSequenceDistance::newD8Instance();
    }
    else if (vm.count("ratio")) {
      	boost::rational<int> ratio;
	std::istringstream iss(vm["ratio"].as<std::string>());
	iss >> ratio;
	//trace.info() << ratio << std::endl;
/*	if (!parseRate(vm["ratio"].as<string>(), &num, &den)) {
	    trace.error() << "Unable to parse num/den from \"%s\"\n"
			  << vm["ratio"].as<string>();
	    exit(-1);
	}*/
	if (ratio < 0 || ratio > 1) {
	    std::cerr <<
		"Invalid ratio " << ratio << std::endl <<
		"correct ratios num/den are between 0 and 1 inclusive" << std::endl;
	    exit(-1);
	}
	dist = NeighborhoodSequenceDistance::newInstance(ratio);
    }
    else if (vm.count("sequence")) {
      std::vector<int> sequence = parseSequence(vm["sequence"].as<std::string>());
	dist = NeighborhoodSequenceDistance::newInstance(sequence);
    }
    //------------------------------------------------------------------------//

    // Output ----------------------------------------------------------------//
    ImageConsumer<GrayscalePixelType> *output;
    {
	std::string outputFile("-");
	std::string outputFormat("");
	bool lineBuffered = false;

	if (vm.count("output")) {
          outputFile = vm["output"].as<std::string>();
	}
	if (vm.count("outputFormat")) {
          outputFormat = vm["outputFormat"].as<std::string>();
	}
	if (vm.count("lineBuffered")) {
	    lineBuffered = true;
	}

	output = createImageWriter(outputFile, outputFormat, lineBuffered);

	if (output == NULL) {
	    std::cerr << "Unable to create image output stream (unrecognized format?)" << std::endl;
	}

	if (vm.count("center")) {
	    output = dist->newDistanceTransformUntranslator(output);
	}
    }
    NeighborhoodSequenceDistanceTransform *dt = dist->newTranslatedDistanceTransform(output);
    //------------------------------------------------------------------------//

    // Input -----------------------------------------------------------------//
    FILE *input;
    {
	std::string inputFile("-");
	if (vm.count("input")) {
          inputFile = vm["input"].as<std::string>();
	}
	if (inputFile == "-") {
	    input = stdin;
	}
	else {
	    input = fopen(inputFile.c_str(), "r");
	    // FIXME: where is fclose?
	    if (input == NULL)
		std::cerr << "Unable to open input stream";
	}
    }
    int inputFormat = 0;

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
	(inputFormat == 0 && readBytes == 8 && png_check_sig(signature, 8))) {

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
	std::cerr << "Input image format not recognized" << std::endl;
    }
    //------------------------------------------------------------------------//

    if (dt != NULL) {
	delete dt;
    }

    fclose(input);
    return 0;
}
