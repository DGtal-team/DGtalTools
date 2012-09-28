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
 * @file BaseDistanceDT.cpp
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
