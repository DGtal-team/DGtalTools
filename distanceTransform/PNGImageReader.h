/*
 *  PNGImageReader.h
 *  LUTBasedNSDistanceTransform
 *  $Id: PNGImageReader.h 124 2012-07-16 10:07:55Z Nicolas.Normand $
 */

#include "ImageFilter.h"
#include "png.h"

class PNGImageReader: public RowImageProducer<BinaryPixelType> {
public:
    PNGImageReader(ImageConsumer<BinaryPixelType>* consumer, FILE *input);

    void produceAllRows(size_t readBytes);

private:
    void startImage(size_t readBytes = 0);

    png_structp _png_ptr;
    png_infop _info_ptr;

    png_uint_32 _cols, _rows;
    FILE *_input;
    BinaryPixelType *_inputRow;
};
