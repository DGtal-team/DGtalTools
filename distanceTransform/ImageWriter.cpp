/*
 *  ImageWriter.cpp
 *  LUTBasedNSDistanceTransform
 *  $Id: ImageWriter.cpp 108 2012-07-11 06:09:48Z Nicolas.Normand $
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
