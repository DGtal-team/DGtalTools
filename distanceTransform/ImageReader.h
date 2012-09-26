/*
 *  ImageReader.h
 *  LUTBasedNSDistanceTransform
 *  $Id: ImageReader.h 112 2012-07-11 07:58:28Z Nicolas.Normand $
 */

#include "LUTBasedNSDistanceTransformConfig.h"

#ifdef WITH_NETPBM
#   include "PBMImageReader.h"
#endif
#ifdef WITH_PNG
#   include "PNGImageReader.h"
#endif
