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
 * @file NeighborhoodSequenceDistance.cpp
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


#include <algorithm>

#include "NeighborhoodSequenceDistance.h"
#include "PeriodicNSDistanceDT.h"
#include "RatioNSDistanceDT.h"
#include "D4DistanceDT.h"
#include "D8DistanceDT.h"

void NeighborhoodSequenceDistanceTransform::rotate() {
    GrayscalePixelType *t = dtLines[2];
    dtLines[2] = dtLines[1];
    dtLines[1] = dtLines[0];
    dtLines[0] = t;
}

void NeighborhoodSequenceDistanceTransform::beginOfImage(int cols, int rows) {
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

void NeighborhoodSequenceDistanceTransform::endOfImage() {
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

NeighborhoodSequenceDistanceTransform::NeighborhoodSequenceDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) :
super(consumer),
_inited(false),
_cols(0) {

    dtLines[0] = NULL;
    dtLines[1] = NULL;
    dtLines[2] = NULL;
}

NeighborhoodSequenceDistanceTransform::~NeighborhoodSequenceDistanceTransform() {
}

NeighborhoodSequenceDistance* NeighborhoodSequenceDistance::newD4Instance() {
    return new D4Distance();
}

NeighborhoodSequenceDistance* NeighborhoodSequenceDistance::newD8Instance() {
    return new D8Distance();
}

NeighborhoodSequenceDistance* NeighborhoodSequenceDistance::newInstance(std::vector<int> sequence) {
    int countOfNeighbors[2] = {0, 0};
    for (std::vector<int>::iterator it = sequence.begin(); it != sequence.end(); it++) {
	BOOST_VERIFY(*it == 1 || *it == 2);
	countOfNeighbors[*it - 1]++;
    }
    if (countOfNeighbors[0] == 0) {
	// d8
	return new D8Distance();
    }
    else if (countOfNeighbors[1] == 0) {
	// d4
	return new D4Distance();
    }
    else {
	return new PeriodicNSDistance(sequence);
    }		
}

NeighborhoodSequenceDistance* NeighborhoodSequenceDistance::newInstance(boost::rational<int> ratio) {
    if (ratio == 0) {
	return new D4Distance();
    }
    else if (ratio == 1) {
	return new D8Distance();
    }
    else {
	return new RatioNSDistance(ratio);
    }	
}
