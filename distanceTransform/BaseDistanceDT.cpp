// $Id: BaseDistanceDT.cpp 106 2012-07-10 20:40:36Z Nicolas.Normand $

#include <assert.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <unistd.h>

#include <algorithm>

#include "BaseDistanceDT.h"

void BaseDistanceTransform::rotate() {
    GrayscalePixelType *t = dtLines[2];
    dtLines[2] = dtLines[1];
    dtLines[1] = dtLines[0];
    dtLines[0] = t;
}

void BaseDistanceTransform::beginOfImage(int cols, int rows) {
    assert(!_inited);
    assert(_cols == 0);
    assert(dtLines[0] == NULL);
    assert(dtLines[1] == NULL);
    assert(dtLines[2] == NULL);

    _cols = cols;
    dtLines[0] = (GrayscalePixelType *) malloc((2 + cols + 1) * sizeof(GrayscalePixelType));
    memset(dtLines[0], 0, (2 + cols + 1) * sizeof(GrayscalePixelType));
    dtLines[1] = (GrayscalePixelType *) malloc((2 + cols + 1) * sizeof(GrayscalePixelType));
    memset(dtLines[1], 0, (2 + cols + 1) * sizeof(GrayscalePixelType));
    dtLines[2] =  (GrayscalePixelType *) malloc((2 + cols + 1) * sizeof(GrayscalePixelType));
    memset(dtLines[2], 0, (2 + cols + 1) * sizeof(GrayscalePixelType));

    _consumer->beginOfImage(cols, rows);

    _inited = true;
}

void BaseDistanceTransform::endOfImage() {
    _consumer->endOfImage();

    free(dtLines[0]);
    free(dtLines[1]);
    free(dtLines[2]);

    dtLines[0] = NULL;
    dtLines[1] = NULL;
    dtLines[2] = NULL;
    _cols = 0;
    _inited = false;
}

BaseDistanceTransform::BaseDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) :
super(consumer),
_inited(false),
_cols(0) {

    dtLines[0] = NULL;
    dtLines[1] = NULL;
    dtLines[2] = NULL;
}

BaseDistanceTransform::~BaseDistanceTransform() {
}
