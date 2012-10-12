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
 * @file RatioNSDistanceDT.cpp
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

#include "RatioNSDistanceDT.h"

#ifndef NDEBUG
#define RBS(num, den, dir, n) (((n) * (num) - (dir) + (den)) / (den) - 1)
#define RBSi(num, den, dir, n) RBS(den, num, !(dir), n)
#define MATHBF2(num, den, n) RBS(num, den, 0, n)
#define MATHBF2i(num, den, n) RBSi(num, den, 0, n)
#define MATHBF1(num, den, n) (1 + RBS((den - num), den, 1, n))
#define MATHBF1i(num, den, n) RBSi((den - num), den, 1, (n - 1))
#define C1(num, den, n) (MATHBF1i(num, den, MATHBF1(num, den, n) + 1) + 1)
#define C2(num, den, n) (MATHBF2i(num, den, MATHBF2(num, den, n) + 1) + 1)
#endif

RatioNSDistance::RatioNSDistance(boost::rational<int> ratio) :
_ratio(ratio),
mbf1(RationalBeattySequence(1 - ratio, ratio.denominator() - 1)),
mbf2(RationalBeattySequence(ratio, 0)),
mbf1i(RationalBeattySequence(1 - ratio, ratio.denominator() - 1).invert()),
mbf2i(RationalBeattySequence(ratio, 0).invert()) {
#ifndef NDEBUG
    RationalBeattySequence B(ratio + 1, 0);
    for (int n = 1; n < 1000; n++) {
	assert(mbf1(n) + mbf2(n) == n);
	assert(mbf1(n) - mbf1(n - 1) == (B(n)-B(n-1) == 1));
	assert(mbf2(n) - mbf2(n - 1) == (B(n)-B(n-1) == 2));
	assert(MATHBF1(ratio.numerator(), ratio.denominator(), n) == mbf1(n));
	assert(MATHBF2(ratio.numerator(), ratio.denominator(), n) == mbf2(n));
	assert(MATHBF1i(ratio.numerator(), ratio.denominator(), n) == mbf1i(n));
	assert(MATHBF2i(ratio.numerator(), ratio.denominator(), n) == mbf2i(n));
    }
#endif
}

NeighborhoodSequenceDistanceTransform* RatioNSDistance::newTranslatedDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) const {
    return new RatioNSDistanceTransform(consumer, _ratio);
}

DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType>* RatioNSDistance::newDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer) const {
    // FIXME: distance max set to 0
    return new RatioNSDistanceTransformUntranslator(consumer, 0, _ratio);
}

void RatioNSDistanceTransform::processRow(const BinaryPixelType *imageRow) {
    int col;
#define N1_SETMINUS_N2_COUNT  1
#define N2_SETMINUS_N1_COUNT  5
#define N1_CAP_N2_COUNT 3
    static vect n1[N1_SETMINUS_N2_COUNT] = {vect(-1, 1)};
    static vect n2[N2_SETMINUS_N1_COUNT] = {			     vect(1, 0), vect(2, 0),
										 vect(2, 1),
								     vect(1, 2), vect(2, 2)};
    static vect n12[N1_CAP_N2_COUNT]	 = {		 vect(0, 1), vect(1, 1),
							 vect(0, 2)};

    for (col = 0; col < _cols; col++) {
	if (imageRow[col] == 0)
	    dtLines[0][col + 2] = 0;
	else {
	    GrayscalePixelType val;
	    GrayscalePixelType dt;
	    int k;

	    val = GRAYSCALE_MAX;
	    for (k = 0; k < N1_SETMINUS_N2_COUNT; k++) {
		assert(n1[k][1] >= 0);
		assert(n1[k][1] <= 2);
		assert(col + 2 - n1[k][0] >= 0);
		assert(col + 2 - n1[k][0] < _cols + 3);
		val = std::min(val, dtLines[n1[k][1]][col + 2 - n1[k][0]]);
	    }
	    //assert(C1(d.num, d.den, (int) val) == d.mbf1i(d.mbf1(val)+1)+1);
	    dt = d.mbf1i(d.mbf1(val)+1)+1;

	    val = GRAYSCALE_MAX;
	    for (k = 0; k < N2_SETMINUS_N1_COUNT; k++) {
		assert(n2[k][1] >= 0);
		assert(n2[k][1] <= 2);
		assert(col + 2 - n2[k][0] >= 0);
		assert(col + 2 - n2[k][0] < _cols + 3);
		val = std::min(val, dtLines[n2[k][1]][col + 2 - n2[k][0]]);
	    }
	    //assert(C2(d.num, d.den, (int) val) == d.mbf2i(d.mbf2(val)+1)+1);
	    dt = std::min((int) dt, d.mbf2i(d.mbf2(val)+1)+1);

	    val = GRAYSCALE_MAX;
	    for (k = 0; k < N1_CAP_N2_COUNT; k++) {
		assert(n12[k][1] >= 0);
		assert(n12[k][1] <= 2);
		assert(col + 2 - n12[k][0] >= 0);
		assert(col + 2 - n12[k][0] < _cols + 3);
		val = std::min(val, dtLines[n12[k][1]][col + 2 - n12[k][0]]);
	    }
	    dt = std::min((int) dt, val + 1);

	    dtLines[0][col + 2] = dt;
	}
    }
    _consumer->processRow(dtLines[0]+2);
    rotate();
}

RatioNSDistanceTransform::RatioNSDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer, boost::rational<int> ratio) :
    NeighborhoodSequenceDistanceTransform(consumer),
    d(ratio) {
}

RatioNSDistanceTransformUntranslator::RatioNSDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer, int dMax, boost::rational<int> ratio) :
    super(consumer, marginRight),
    _dMax(dMax),
    _imageDMax(0),
    d(ratio) {
}

RatioNSDistanceTransformUntranslator::~RatioNSDistanceTransformUntranslator() {
}

void RatioNSDistanceTransformUntranslator::beginOfImage(int cols, int rows) {
    assert(_imageDMax == 0);
    _imageDMax = (_dMax == 0) ? INT_MAX : _dMax;
    _imageDMax = std::min(_imageDMax, (cols + 1) / 2);
    _imageDMax = std::min(_imageDMax, (rows + 1) / 2);
    super::beginOfImage(cols, rows, _imageDMax + 1);
}

void RatioNSDistanceTransformUntranslator::endOfImage() {
    super::endOfImage();
    _imageDMax = 0;
}

// untranslate is called once for each row of the input image, plus one extra time
// with null-valued translated DT to flush all DT values
void RatioNSDistanceTransformUntranslator::processRow(const GrayscalePixelType* inputRow) {
    int dtmax = 1;  // Not 0 to avoid outputing the extra row

    super::processRow(inputRow);

    for (int col = 0; col < _cols; col++) {

	int dtn, dtp;

	dtp = _tdtRows[1][col];
	if (_tdtRows[0][col] == 0) {
	    assert(_outputRows[(_curRow + _dtRowCount) % _dtRowCount][col] == (GrayscalePixelType) -1);
	    _outputRows[(_curRow + _dtRowCount) % _dtRowCount][col] = 0;
	}

	dtn = _tdtRows[0][col];
	dtmax = std::max(dtmax, dtn);
	dtn = std::max(0, dtn - 1);
	//assert(C1(d.num, d.den, dtn) == d.mbf1i(d.mbf1(dtn)+1)+1);
	//assert(MATHBF2(d.num, d.den, d.mbf1i(d.mbf1(dtn)+1)) == d.mbf2(d.mbf1i(d.mbf1(dtn)+1)));
	for (int r = d.mbf1i(d.mbf1(dtn)+1)+1, dx = d.mbf2(r-1);
	     r <= dtp;
	     r = d.mbf1i(d.mbf1(r)+1)+1) {

	    int dy = r - 1;
	    assert (dx == d.mbf2(r-1));
	    assert(_curRow - 1 - dy >= 0);
	    assert(_outputRows[(_curRow - 1 - dy) % _dtRowCount][col - dx] == (GrayscalePixelType) -1);
	    _outputRows[(_curRow - 1 - dy) % _dtRowCount][col - dx] = r;
	    // Let s be the next radius where neighborhood 1 is used
	    // s = r + c1[r % period]
	    // between r and s, neighborhood 2 is used s - r - 1 times
	    // so dx is raised by s - r - 1 = c1[r % period] - 1
	    //assert(dx + c1[r % period] - 1 == mathbf2(d, r + c1[r % period] - 1));
	    //assert(dx + C1(d.num, den, r) - r - 1 == MATHBF2(d.num, den, C1(d.num, den, r) - 1));
	    //dx += C1(d.num, den, r) - r - 1;
	    //assert(MATHBF2(d.num, d.den, C1(d.num, d.den, r) - 1) == d.mbf2(d.mbf1i(d.mbf1(r)+1)));
	    dx = d.mbf2(d.mbf1i(d.mbf1(r)+1));
	    //assert(C1(d.num, d.den, r) == d.mbf1i(d.mbf1(r)+1)+1);
	}

	dtn = _tdtRows[0][col + 1];
	dtn = std::max(0, dtn - 1);
	//assert(C2(d.num, d.den, dtn) == d.mbf2i(d.mbf2(dtn)+1)+1);
	//assert(MATHBF2(d.num, d.den, d.mbf2i(d.mbf2(dtn)+1)) == d.mbf2(d.mbf2i(d.mbf2(dtn)+1)));
	for (int r = d.mbf2i(d.mbf2(dtn)+1)+1, dx = d.mbf2(d.mbf2i(d.mbf2(dtn)+1));
	     r <= dtp;
	     r = d.mbf2i(d.mbf2(r)+1)+1) {

	    int dy = r - 1;
	    //assert(dx == MATHBF2(d.num, d.den, r - 1));

	    assert(_curRow - 1 - dy >= 0);
	    assert(_outputRows[(_curRow - 1 - dy) % _dtRowCount][col - dx] == (GrayscalePixelType) -1);
	    _outputRows[(_curRow - 1 - dy) % _dtRowCount][col - dx] = r;
	    // Next time we use neighborhood 2, dx is increased by one
	    //assert(dx + 1 == MATHBF2(d.num, d.den, C2(d.num, d.den, r) - 1));
	    dx++;
	    //assert(C2(d.num, d.den, r) == d.mbf2i(d.mbf2(r)+1)+1);
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
