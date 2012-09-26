/*
 *  PNGImageWriter.h
 *  LUTBasedNSDistanceTransform
 *  $Id: PNGImageWriter.h 108 2012-07-11 06:09:48Z Nicolas.Normand $
 */

#include "ImageFilter.h"
#include <png.h>

class PNGImageWriter : public ImageConsumer<GrayscalePixelType> {
public:
    PNGImageWriter(FILE* output, bool lineBuffered = false);

    void beginOfImage(int cols, int rows);

    void endOfImage();

    void processRow(const GrayscalePixelType* inputRow);

protected:
    png_structp _png_ptr;
    png_infop _info_ptr;

    FILE* _output;
    bool _lineBuffered;
};
