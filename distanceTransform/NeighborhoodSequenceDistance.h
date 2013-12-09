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
 * @file NeighborhoodSequenceDistance.h
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

#ifndef BASE_DISTANCE_DT_H
#define BASE_DISTANCE_DT_H

#include "ImageFilter.h"

#include "string.h"

#include <vector>
#include <boost/rational.hpp>
#include <DGtal/kernel/PointVector.h>

typedef DGtal::PointVector<2, int> vect;

template <typename inputPixelType, typename outputPixelType>
class DistanceTransformUntranslator: public ImageFilter<inputPixelType, outputPixelType> {
private:
    typedef ImageFilter<inputPixelType, outputPixelType> super;
protected:
    DistanceTransformUntranslator(ImageConsumer<outputPixelType>* consumer, int rightMargin) :
    super(consumer),
    _inited(false),
    _cols(0),
    _rightMargin(rightMargin),
    _curRow(0),
    _dtRowCount(0),
    _outputRows(NULL) {
	_tdtRows[0] = NULL;
	_tdtRows[1] = NULL;
    }

    void beginOfImage(int cols, int rows, int dtRowCount) {
	assert(!_inited);
	assert(_cols == 0);
	assert(_curRow == 0);
	assert(_dtRowCount == 0);
	assert(_outputRows == NULL);
	assert(_tdtRows[0] == NULL);
	assert(_tdtRows[1] == NULL);

	_dtRowCount = dtRowCount;
	_outputRows = (GrayscalePixelType **)malloc(_dtRowCount * sizeof(GrayscalePixelType *));
	assert(_outputRows);
	_tdtRows[0] = (inputPixelType *) malloc((cols + 1) * sizeof(inputPixelType));
	memset(_tdtRows[0], 0, (cols + 1) * sizeof(GrayscalePixelType));
	assert(_tdtRows[0]);
	_tdtRows[1] = (inputPixelType *) malloc((cols + 1) * sizeof(inputPixelType));
	assert(_tdtRows[1]);
	memset(_tdtRows[1], 0, (cols + 1) * sizeof(GrayscalePixelType));

	_curRow = 1;	// Start at 1 to avoid modulo of negative problem
	_outRow = _curRow;

	_cols = cols;

	for (int row = 0; row < _dtRowCount; row++) {
	    _outputRows[row] = (outputPixelType *) malloc(cols * sizeof(outputPixelType));
    #ifndef NDEBUG
	    // Set all values to -1 to check later that each pixel is assigned a
	    // value exactly once.
	    // - assert pixel is -1 before setting a value
	    // - assert pixel is not -1 before outputting and resetting it
	    for (int col = 0; col < cols; col++) {
		_outputRows[row][col] = -1;
	    }
    #endif
	}

	super::beginOfImage(cols, rows);

	_inited = true;
    }

    void processRow(const inputPixelType* inputRow) {
	inputPixelType *t = _tdtRows[1];
	_tdtRows[1] = _tdtRows[0];
	_tdtRows[0] = t;
	if (inputRow == NULL)
	    memset(_tdtRows[0], 0, _cols * sizeof(inputPixelType));
	else
	    memcpy(_tdtRows[0], inputRow, _cols * sizeof(inputPixelType));
    }

    void endOfImage() {
	this->processRow(NULL);

	_inited = false;
	_cols = 0;
	_curRow = 0;
	_dtRowCount = 0;
	for (int row = 0; row < _dtRowCount; row++) {
	    free(_outputRows[row]);
	}
	free(_outputRows);
	_outputRows = NULL;
	free(_tdtRows[0]);
	free(_tdtRows[1]);
	_tdtRows[0] = NULL;
	_tdtRows[1] = NULL;

	super::endOfImage();
    }

    bool _inited;
    int _cols;
    const int _rightMargin;
    int _curRow;
    int _outRow;
    int _dtRowCount;
    outputPixelType **_outputRows;
    inputPixelType *_tdtRows[2];
};

class NeighborhoodSequenceDistanceTransform: public ImageFilter<BinaryPixelType, GrayscalePixelType> {
private:
    typedef ImageFilter<BinaryPixelType, GrayscalePixelType> super;
public:
    NeighborhoodSequenceDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer);
    ~NeighborhoodSequenceDistanceTransform();

    void beginOfImage(int cols, int rows);
    //void processRow(const bit* inputRow);
    void endOfImage();

protected:
    void rotate();

    bool _inited;
    int _cols;
    GrayscalePixelType* dtLines[3];
};

/**
 * \brief Abstract factory for the instanciation of neighborhood sequence
 * distances.
 *
 * Instances of concrete subclasses provides factory methods for the creation
 * of image filters aimed at:
 * - computing a translated image transform,
 * - recentering the translated image transform.
 */
class NeighborhoodSequenceDistance {
public:
    virtual NeighborhoodSequenceDistanceTransform* newTranslatedDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) const = 0;
    virtual DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType>* newDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer) const = 0;

    static NeighborhoodSequenceDistance* newD4Instance();
    static NeighborhoodSequenceDistance* newD8Instance();
    static NeighborhoodSequenceDistance* newInstance(boost::rational<int> ratio);
    static NeighborhoodSequenceDistance* newInstance(std::vector<int> sequence);
};

#endif
