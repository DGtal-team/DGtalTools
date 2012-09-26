// $Id: PGMImageWriter.h 106 2012-07-10 20:40:36Z Nicolas.Normand $

extern "C" {
#include <pam.h>
#include <assert.h>
}
//#include <algorithm>
//#include <iostream>

#include "ImageFilter.h"

class PGMImageWriter : public ImageConsumer<GrayscalePixelType> {
public:
    PGMImageWriter(FILE* output, int format = 0);

    void beginOfImage(int cols, int rows);
    void processRow(const GrayscalePixelType* inputRow);
    void endOfImage();

protected:
    int _cols;
    int _format;
    FILE* _output;
    struct pam _outpam;
    tuple * _tuplerow;
};
