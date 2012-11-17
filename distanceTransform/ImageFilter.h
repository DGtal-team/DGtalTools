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
 * @file ImageFilter.h
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

#ifndef IMAGE_FILTER_H
#define IMAGE_FILTER_H

#include <sys/types.h>

#include <assert.h>
#include <algorithm>
#include <iostream>

#include "LUTBasedNSDistanceTransformConfig.h"

void freerow(void *);

/**
 * \brief Template class for image consumers. An ImageConsumer accepts image
 * data one row at a time from top to bottom.
 *
 * A call to beginOfImage() initialises the ImageConsumer,
 * followed by as many calls to processRow() as the number of rows in the image
 * and a final call to endOfImage().
 */
template <typename inputPixelType>
class ImageConsumer {
public:
    virtual ~ImageConsumer() { }

    virtual void processRow(const inputPixelType* inputRow) = 0;
    virtual void beginOfImage(int cols, int rows) = 0;
    virtual void endOfImage() = 0;
};

/**
 * \brief Template class for image producers. A RowImageProducer generates image
 * data row by row from top to bottom.
 *
 * The data is pushed to the ImageConsumer attached to it.
 */
template <typename outputPixelType>
class RowImageProducer {
public:
    RowImageProducer(ImageConsumer<outputPixelType>* consumer) : _consumer(consumer) {
	if (_consumer == NULL) {
	    exit(1);
	}
    }
    bool hasMoreRows();
    void produceRow();
    void produceAllRows();
    virtual ~RowImageProducer() { delete _consumer; }
protected:
    ImageConsumer<outputPixelType>* _consumer;
};

/**
 * \brief Template class for image filters. An ImageFilter is an image consumer
 * and producer.
 *
 * It redirects calls to beginOfImage(), processRow() and endOfImage() to the
 * next ImageConsumer (optionally modifying the size and content of the image).
 */
template <typename inputPixelType, typename outputPixelType>
class ImageFilter: public ImageConsumer<inputPixelType> {
public:
    ImageFilter(ImageConsumer<outputPixelType>* consumer) : _consumer(consumer) {
	if (_consumer == NULL) {
	    exit(1);
	}
    }
    virtual ~ImageFilter() { delete _consumer; }
    void beginOfImage(int cols, int rows) {_consumer->beginOfImage(cols, rows);}
    void endOfImage() {_consumer->endOfImage();}
protected:
    ImageConsumer<outputPixelType>* _consumer;
};

/**
 * \brief A TeeImageFilter forwards the image data to two ImageConsumer intances
 * instead of one.
 */
template <typename inputPixelType, typename outputPixelType>
class TeeImageFilter: public ImageFilter<inputPixelType, outputPixelType> {
public:
    typedef ImageFilter<inputPixelType, outputPixelType> super;

    TeeImageFilter(ImageConsumer<inputPixelType>* consumer, ImageConsumer<inputPixelType>* consumer2) :
	ImageFilter<inputPixelType, outputPixelType>(consumer),
	_consumer2(consumer2) {}
    ~TeeImageFilter() { delete _consumer2; }

    void beginOfImage(int cols, int rows) {
	super::beginOfImage(cols, rows);

	if (_consumer2) {
	    _consumer2->beginOfImage(cols, rows);
	}
    }

    void processRow(const outputPixelType* inputRow);
    /*
    void processRow(const outputPixelType* inputRow) {
	_consumer->processRow(inputRow);

	if (_consumer2) {
	    _consumer2->processRow(inputRow);
	}
    }*/

    void endOfImage() {
	super::endOfImage();

	if (_consumer2) {
	    _consumer2->endOfImage();
	}
    }
protected:
    ImageConsumer<outputPixelType>* _consumer2;
};

#endif
