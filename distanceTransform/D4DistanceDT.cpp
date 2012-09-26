// $Id: D4DistanceDT.cpp 124 2012-07-16 10:07:55Z Nicolas.Normand $

#include <assert.h>

#include <algorithm>

#include "D4DistanceDT.h"

BaseDistanceTransform* D4Distance::newTranslatedDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) const {
    return new D4DistanceTransform(consumer);
}

DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType>* D4Distance::newDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer) const {
    // FIXME: distance max set to 0
    return new D4DistanceTransformUntranslator(consumer, 0);
}

void D4DistanceTransform::processRow(const BinaryPixelType *imageRow) {
    int col;
#define N1_COUNT  4
    static vect n1[N1_COUNT]   = {{-1, 1}, {0, 1}, {1, 1}, {0, 2}};

    for (col = 0; col < _cols; col++) {
	if (imageRow[col] == 0)
	    dtLines[0][col + 2] = 0;
	else {
	    GrayscalePixelType val;
	    int k;

	    val = GRAYSCALE_MAX;
	    for (k = 0; k < N1_COUNT; k++) {
		assert(n1[k].y >= 0);
		assert(n1[k].y <= 2);
		assert(col + 2 - n1[k].x >= 0);
		assert(col + 2 - n1[k].x < _cols + 3);
		val = std::min(val, dtLines[n1[k].y][col + 2 - n1[k].x]);
	    }

	    dtLines[0][col + 2] = val + 1;
	}
    }
    _consumer->processRow(dtLines[0]+2);
    rotate();
}

D4DistanceTransform::D4DistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) : BaseDistanceTransform(consumer) {
}

D4DistanceTransformUntranslator::D4DistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer, int dMax) :
super(consumer, 0),
_dMax(dMax),
_imageDMax(0) {
}

D4DistanceTransformUntranslator::~D4DistanceTransformUntranslator() {
}

void D4DistanceTransformUntranslator::beginOfImage(int cols, int rows) {
    assert(_imageDMax == 0);
    _imageDMax = (_dMax == 0) ? INT_MAX : _dMax;
    _imageDMax = std::min(_imageDMax, (cols + 1) / 2);
    _imageDMax = std::min(_imageDMax, (rows + 1) / 2);
    super::beginOfImage(cols, rows, _imageDMax + 1);
}

void D4DistanceTransformUntranslator::endOfImage() {
    super::endOfImage();
    _imageDMax = 0;
}

// untranslate is called once for each row of the input image, plus one extra time
// with null-valued translated DT to flush all DT values
void D4DistanceTransformUntranslator::processRow(const GrayscalePixelType* inputRow) {
    int dtmax = 1;  // Not 0 to avoid outputing the extra row

    super::processRow(inputRow);

    for (int col = 0; col < _cols; col++) {

	int dtn, dtp;
	int dy;

	dtp = _tdtRows[1][col];
	if (_tdtRows[0][col] == 0) {
	    assert(_outputRows[(_curRow + _dtRowCount) % _dtRowCount][col] == (GrayscalePixelType) -1);
	    _outputRows[(_curRow + _dtRowCount) % _dtRowCount][col] = 0;
	}

	dtn = _tdtRows[0][col];
	dtmax = std::max(dtmax, dtn);
	dtn = std::max(0, dtn - 1);

	for (int r = dtn + 1;
	     r <= dtp;
	     r++) {

	    dy = r - 1;
	    assert(_curRow - 1 - dy >= 0);
	    assert(_outputRows[(_curRow - 1 - dy) % _dtRowCount][col] == (GrayscalePixelType) -1);
	    _outputRows[(_curRow - 1 - dy) % _dtRowCount][col] = r;
	    //printf("%d ", r);
	}
	//printf("\n");
    }
    _curRow++;
    //printf("max: %d\n", dtmax);
    for (; _outRow < _curRow - dtmax; _outRow++) {
	_consumer->processRow(_outputRows[_outRow % _dtRowCount]);
#ifndef NDEBUG
	for (int col = 0; col < _cols; col++) {
	    assert(_outputRows[_outRow % _dtRowCount][col] != (GrayscalePixelType) -1);
	    _outputRows[_outRow % _dtRowCount][col] = -1;
	}
#endif
    }
}
