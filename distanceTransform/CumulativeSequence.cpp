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
 * @file CumulativeSequence.cpp
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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include <algorithm>

#include "CumulativeSequence.h"

int CumulativeOfPeriodicSequence::equals(CumulativeOfPeriodicSequence& seq2) {
    if (_sequence.size() != seq2._sequence.size())
	return 0;
#warning Do this properly
    return memcmp(_sequence.data(), seq2._sequence.data(), sizeof(CumulativeOfPeriodicSequence) + _sequence.size() * sizeof(int)) == 0;
}

int mod(int a, int b) {
    return ((a % b) + b) % b;
}

CumulativeOfPeriodicSequence *CumulativeOfPeriodicSequence::inverse() {
    //CumulativeOfPeriodicSequence *inv = CumulativeOfPeriodicSequenceCreate(_sequence[seq->period - 1], 0, NULL);
    CumulativeOfPeriodicSequence *inv = new CumulativeOfPeriodicSequence(_sequence[_sequence.size() - 1]);
    inv->_offset = 0;

    int xx, yy = 0;
    yy = mod(_offset, inv->_sequence.size());
    //yy = CumulativeOfPeriodicSequenceValueAtIndex(seq, 1);
    //yy %= inv->_sequence.size();
    for (xx = 0; xx < _sequence.size(); xx++) {
	//yy = CumulativeOfPeriodicSequenceValueAtIndex(seq, xx + 1);
	//yy %= inv->_sequence.size();
	inv->_sequence[yy]++;
	yy += _sequence[xx] - (xx > 0 ? _sequence[xx-1] : 0);
	yy %= inv->_sequence.size();
    }

    for (xx = 1; xx < inv->_sequence.size(); xx++) {
	inv->_sequence[xx] += inv->_sequence[xx-1];
    }
    
    if (valueAtIndex(1) > 0) {
	yy = valueAtIndex(1);
	// Find first increasing index
	xx = 1;
	while (valueAtIndex(xx) == yy) xx++;
	// First positive term in inverse must be equal to xx - 1
	// i.e., inv->value[yy % inv->_sequence.size()] + (yy / inv->_sequence.size()) * inv->_sequence[inv->_sequence.size() - 1] + inv->_offset == xx - 1
	inv->_offset = xx - 1 - inv->_sequence[yy % inv->_sequence.size()] - (yy / inv->_sequence.size()) * inv->_sequence[inv->_sequence.size() - 1];
	//TODO: assert(CumulativeOfPeriodicSequenceValueAtIndex(inv, yy) == 0);
	//TODO: assert(CumulativeOfPeriodicSequenceValueAtIndex(inv, yy+1) == xx - 1);
    }
    else {
	// Find first positive value
	xx = 0;
	while (valueAtIndex(xx) == 0) xx++;
	// seq(xx) > 0 and seq(xx-1) <= 0 then seqinv(1) = xx-1;
	// i.e., inv->value[0] + inv->_offset == xx-1;
	
	inv->_offset = xx - 1 - inv->_sequence[0];
    }

    return inv;
}

void CumulativeOfPeriodicSequenceFree(CumulativeOfPeriodicSequence *seq) {
    free(seq);
}

int CumulativeOfPeriodicSequence::valueAtIndex(int i) {
    assert(i >= 0);
    if (i == 0) return 0;
    i--;
    return std::max(_sequence[(i % _sequence.size())] +
		    (long int) (i / _sequence.size()) * _sequence[_sequence.size() - 1] +
		    _offset, 0L);
}

void CumulativeOfPeriodicSequence::print() {
    int i;

    //printf("period: %d, offset: %d, values: ", _sequence.size(), seq->offset);
    printf("(%d", _sequence[0]);
    for (i = 1; i < _sequence.size(); i++) {
	//printf("%d ", _sequence[i]);
	printf(",%d", _sequence[i] - (i > 0 ? _sequence[i-1] : 0));
    }
    printf(")");
    if (_offset > 0) {
	printf("(+%d)", _offset);
    }
    else if (_offset < 0) {
	printf("(%d)", _offset);
    }
    printf("\n");
}
