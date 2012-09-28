// $Id: CumulativeSequence.cpp 94 2012-07-04 07:32:53Z Nicolas.Normand $

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
