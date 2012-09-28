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

//#include "MMChainConfig.h"
#include "ImageWriter.h"

bool checkFormat(char const *fromname, char const *format) {
    if (fromname == NULL)
	return false;

    int l = strlen(format);
    return strncasecmp(fromname, format, l) == 0 &&
	   (fromname[l] == '\0' || fromname[l] == ':');
}

ImageConsumer<GrayscalePixelType> *createImageWriter(char const *filename, char const *format, bool lineBuffered) {
    FILE *output = NULL;

    if (format == NULL) {
	format = index(filename, ':');
	if (format != NULL) {
	    char const *t = format;
	    format = filename;
	    filename = t + 1;
	}
    }

    if (strcmp(filename, "-") == 0) {
	output = stdout;
    }
    else {
	output = fopen(filename, "w");
	// FIXME: where is fclose?
	if (output == NULL)
	    return NULL;
    }
			 
    if (format == NULL && filename != NULL) {
	// Guess format from file extension
	format = rindex(filename, '.');
	if (format != NULL)
	    format++;
    }

#ifdef WITH_NETPBM
    if (checkFormat(format, "pgm")) {
	return new PGMImageWriter(output, 1);
    }
#endif
#ifdef WITH_PNG
    if (checkFormat(format, "png")) {
	return new PNGImageWriter(output, lineBuffered);
    }
#endif

    // We can't obey the specified format, give up
    if (format != NULL) {
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
