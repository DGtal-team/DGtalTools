// $Id: CumulativeSequence.cpp 94 2012-07-04 07:32:53Z Nicolas.Normand $

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include <algorithm>

#include "CumulativeSequence.h"

struct CumulativeOfPeriodicSequence {
    int period;
    int offset;
    int values[0];
};

CumulativeOfPeriodicSequence *CumulativeOfPeriodicSequenceCreate(int period, int offset, int *values) {
    CumulativeOfPeriodicSequence *seq = (CumulativeOfPeriodicSequence *)malloc(sizeof(CumulativeOfPeriodicSequence) + period * sizeof(int));
    if (seq != NULL) {
	int i;

	seq->period = period;
	seq->offset = offset;

	for (i = 0; i < period; i++) {
	    if (values != NULL)
		seq->values[i] = values[i] + (i > 0 ? seq->values[i-1] : 0);
	    else
		seq->values[i] = 0;
	}
    }
    return seq;
}

int CumulativeOfPeriodicSequenceEquals(CumulativeOfPeriodicSequence *seq1, CumulativeOfPeriodicSequence *seq2) {
    if (seq1->period != seq2->period)
	return 0;
    return memcmp(seq1, seq2, sizeof(CumulativeOfPeriodicSequence) + seq1->period * sizeof(int)) == 0;
}

int mod(int a, int b) {
    return ((a % b) + b) % b;
}

CumulativeOfPeriodicSequence *CumulativeOfPeriodicSequenceCreateInverse(CumulativeOfPeriodicSequence *seq) {
    CumulativeOfPeriodicSequence *inv = CumulativeOfPeriodicSequenceCreate(seq->values[seq->period - 1], 0, NULL);

    int xx, yy = 0;
    yy = mod(seq->offset, inv->period);
    //yy = CumulativeOfPeriodicSequenceValueAtIndex(seq, 1);
    //yy %= inv->period;
    for (xx = 0; xx < seq->period; xx++) {
	//yy = CumulativeOfPeriodicSequenceValueAtIndex(seq, xx + 1);
	//yy %= inv->period;
	inv->values[yy]++;
	yy += seq->values[xx] - (xx > 0 ? seq->values[xx-1] : 0);
	yy %= inv->period;
    }

    for (xx = 1; xx < inv->period; xx++) {
	inv->values[xx] += inv->values[xx-1];
    }
    
    if (CumulativeOfPeriodicSequenceValueAtIndex(seq, 1) > 0) {
	yy = CumulativeOfPeriodicSequenceValueAtIndex(seq, 1);
	// Find first increasing index
	xx = 1;
	while (CumulativeOfPeriodicSequenceValueAtIndex(seq, xx) == yy) xx++;
	// First positive term in inverse must be equal to xx - 1
	// i.e., inv->value[yy % inv->period] + (yy / inv->period) * inv->values[inv->period - 1] + inv->offset == xx - 1
	inv->offset = xx - 1 - inv->values[yy % inv->period] - (yy / inv->period) * inv->values[inv->period - 1];
	assert(CumulativeOfPeriodicSequenceValueAtIndex(inv, yy) == 0);
	assert(CumulativeOfPeriodicSequenceValueAtIndex(inv, yy+1) == xx - 1);
    }
    else {
	// Find first positive value
	xx = 0;
	while (CumulativeOfPeriodicSequenceValueAtIndex(seq, xx) == 0) xx++;
	// seq(xx) > 0 and seq(xx-1) <= 0 then seqinv(1) = xx-1;
	// i.e., inv->value[0] + inv->offset == xx-1;
	
	inv->offset = xx - 1 - inv->values[0];
    }

    return inv;
}

void CumulativeOfPeriodicSequenceFree(CumulativeOfPeriodicSequence *seq) {
    free(seq);
}

int CumulativeOfPeriodicSequenceValueAtIndex(CumulativeOfPeriodicSequence *seq, int i) {
    assert(i >= 0);
    if (i == 0) return 0;
    i--;
    return std::max(seq->values[(i % seq->period)] +
		    (i / seq->period) * seq->values[seq->period - 1] +
		    seq->offset, 0);
}

void CumulativeOfPeriodicSequencePrint(CumulativeOfPeriodicSequence *seq) {
    int i;

    //printf("period: %d, offset: %d, values: ", seq->period, seq->offset);
    printf("(%d", seq->values[0]);
    for (i = 1; i < seq->period; i++) {
	//printf("%d ", seq->values[i]);
	printf(",%d", seq->values[i] - (i > 0 ? seq->values[i-1] : 0));
    }
    printf(")");
    if (seq->offset > 0) {
	printf("(+%d)", seq->offset);
    }
    else if (seq->offset < 0) {
	printf("(%d)", seq->offset);
    }
    printf("\n");
}
