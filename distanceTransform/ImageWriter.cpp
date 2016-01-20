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
 * @file ImageWriter.cpp
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

#include <boost/algorithm/string/predicate.hpp>
#include <stdio.h>
#include "ImageWriter.h"

/**
 * Creates an ImageWriter.
 *
 * The output image format is determined, in that order, by:
 * - the **format** parameter if not NULL,
 * - a prefix ended by ':' in the file format (*e.g.* 'png:filename'),
 * - the file extension,
 * If one of these methods specifies a format that is not available, no
 * ImageWriter is created and the function return NULL.
 * Il no format is speficied at all, the default format is used in the last
 * resort.
 */
ImageConsumer<GrayscalePixelType> *createImageWriter(std::string filename, std::string format, bool lineBuffered) {
    FILE *output = NULL;

    // Format wasn't specified in arguments, check if there is a prefix for it.
    if (format == "") {
	size_t n = filename.find(':');
	if (n != std::string::npos) {
	    format = filename.substr(0, n);
	    filename = filename.substr(n+1);
	}
    }

    // Format wasn't specified in arguments nor in the filename prefix, check if
    // there the file has an extension
    if (format == "") {
	size_t n = filename.rfind('.');
	if (n != std::string::npos)
	    format = filename.substr(n+1);
    }

    if (filename == "-") {
	output = stdout;
    }
    else {
	output = fopen(filename.c_str(), "w");
	// FIXME: where is fclose?
	if (output == NULL)
	    return NULL;
    }

#ifdef WITH_NETPBM
    if (boost::iequals(format, "pgm")) {
	return new PGMImageWriter(output, lineBuffered);
    }
#endif
#ifdef WITH_PNG
    if (boost::iequals(format, "png")) {
	return new PNGImageWriter(output, lineBuffered);
    }
#endif

    if (format != "") {
	// We can't obey the specified format, give up
	return NULL;
    }

    // No format specified, use default
#ifdef WITH_NETPBM
    return new PGMImageWriter(output);
#else
#   ifdef WITH_PNG
    return new PNGImageWriter(output);
#   endif
#endif

    return NULL;
}
