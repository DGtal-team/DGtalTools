/*
 *  ImageWriter.h
 *  LUTBasedNSDistanceTransform
 *  $Id: ImageWriter.h 112 2012-07-11 07:58:28Z Nicolas.Normand $
 */

#include "LUTBasedNSDistanceTransformConfig.h"

#include "ImageFilter.h"
#ifdef WITH_NETPBM
#   include "PGMImageWriter.h"
#endif
#ifdef WITH_PNG
#   include "PNGImageWriter.h"
#endif

ImageConsumer<GrayscalePixelType>* createImageWriter(char const *filename = NULL, char const *format = NULL, bool lineBuffered=false);
