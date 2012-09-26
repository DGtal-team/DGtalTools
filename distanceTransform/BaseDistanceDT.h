// $Id: BaseDistanceDT.h 106 2012-07-10 20:40:36Z Nicolas.Normand $

#ifndef BASE_DISTANCE_DT_H
#define BASE_DISTANCE_DT_H

#include "ImageFilter.h"

typedef struct {
    int x;
    int y;
} vect;

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

	_curRow = 1;	// Start at 1 to avoid modulus of negative problem
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

class BaseDistanceTransform: public ImageFilter<BinaryPixelType, GrayscalePixelType> {
private:
    typedef ImageFilter<BinaryPixelType, GrayscalePixelType> super;
public:
    BaseDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer);
    ~BaseDistanceTransform();

    void beginOfImage(int cols, int rows);
    //void processRow(const bit* inputRow);
    void endOfImage();

protected:
    void rotate();

    bool _inited;
    int _cols;
    GrayscalePixelType* dtLines[3];
};

class BaseDistance {
public:
    virtual BaseDistanceTransform* newTranslatedDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) const = 0;
    virtual DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType>* newDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer) const = 0;
};

#endif
