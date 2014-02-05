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
 * @file D8DistanceDT.cpp
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

#include <assert.h>
#include <limits.h>

#include <algorithm>

#include "D8DistanceDT.h"

NeighborhoodSequenceDistanceTransform* D8Distance::newTranslatedDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) const {
    return new D8DistanceTransform(consumer);
}

DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType>* D8Distance::newDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer) const {
    // FIXME: distance max set to 0
    return new D8DistanceTransformUntranslator(consumer, 0);
}

void D8DistanceTransform::processRow(const BinaryPixelType *imageRow) {
    int col;
#define N2_COUNT  8
    static vect n2[N2_COUNT] = {	    vect(1, 0), vect(2, 0),
				vect(0, 1), vect(1, 1), vect(2, 1),
				vect(0, 2), vect(1, 2), vect(2, 2)};

    for (col = 0; col < _cols; col++) {
	if (imageRow[col] == 0)
	    dtLines[0][col + 2] = 0;
	else {
	    GrayscalePixelType val;
	    int k;

	    val = GRAYSCALE_MAX;
	    for (k = 0; k < N2_COUNT; k++) {
		assert(n2[k][1] >= 0);
		assert(n2[k][1] <= 2);
		assert(col + 2 - n2[k][0] >= 0);
		assert(col + 2 - n2[k][0] < _cols + 3);
		val = std::min(val, dtLines[n2[k][1]][col + 2 - n2[k][0]]);
	    }

	    dtLines[0][col + 2] = val + 1;
	}
    }
    _consumer->processRow(dtLines[0]+2);
    rotate();
}

D8DistanceTransform::D8DistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) :
NeighborhoodSequenceDistanceTransform(consumer) {
}

D8DistanceTransformUntranslator::D8DistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer, int dMax) :
super(consumer, 1),
_dMax(dMax),
_imageDMax(0) {
}

D8DistanceTransformUntranslator::~D8DistanceTransformUntranslator() {
}

void D8DistanceTransformUntranslator::beginOfImage(int cols, int rows) {
    assert(_imageDMax == 0);
    _imageDMax = (_dMax == 0) ? INT_MAX : _dMax;
    _imageDMax = std::min(_imageDMax, (cols + 1) / 2);
    _imageDMax = std::min(_imageDMax, (rows + 1) / 2);
    super::beginOfImage(cols, rows, _imageDMax + 1);
}

void D8DistanceTransformUntranslator::endOfImage() {
    super::endOfImage();
    _imageDMax = 0;
}

// untranslate is called once for each row of the input image, plus one extra time
// with null-valued translated DT to flush all DT values
void D8DistanceTransformUntranslator::processRow(const GrayscalePixelType* inputRow) {
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

	dtn = _tdtRows[0][col + 1];
	dtn = std::max(0, dtn - 1);

	for (int r = dtn + 1, dx = r - 1;
	     r <= dtp;
	     r++) {

	    dy = r - 1;

	    assert(_curRow - 1 - dy >= 0);
	    assert(_outputRows[(_curRow - 1 - dy) % _dtRowCount][col - dx] == (GrayscalePixelType) -1);
	    _outputRows[(_curRow - 1 - dy) % _dtRowCount][col - dx] = r;
	    dx++;
	}
    }
    _curRow++;
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
