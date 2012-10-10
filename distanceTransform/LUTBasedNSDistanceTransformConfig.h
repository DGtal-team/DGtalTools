
#ifndef LUBBasedNSDistanceTransformConfig_H
#define LUBBasedNSDistanceTransformConfig_H

#include <limits>

typedef unsigned char  BinaryPixelType;
typedef unsigned short GrayscalePixelType;

const static int GRAYSCALE_MAX = std::numeric_limits<GrayscalePixelType>::max();
#define BINARY_WHITE_PIXEL 0
#define BINARY_BLACK_PIXEL 1

#endif
