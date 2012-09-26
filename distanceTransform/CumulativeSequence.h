// $Id: CumulativeSequence.h 94 2012-07-04 07:32:53Z Nicolas.Normand $

extern "C" {
typedef struct CumulativeOfPeriodicSequence CumulativeOfPeriodicSequence;

CumulativeOfPeriodicSequence *CumulativeOfPeriodicSequenceCreate(int period, int offset, int *values);

CumulativeOfPeriodicSequence *CumulativeOfPeriodicSequenceCreateInverse(CumulativeOfPeriodicSequence *seq);

void CumulativeOfPeriodicSequenceFree(CumulativeOfPeriodicSequence *seq);

int CumulativeOfPeriodicSequenceValueAtIndex(CumulativeOfPeriodicSequence *seq, int i);

int CumulativeOfPeriodicSequenceEquals(CumulativeOfPeriodicSequence *seq1, CumulativeOfPeriodicSequence *seq2);

void CumulativeOfPeriodicSequencePrint(CumulativeOfPeriodicSequence *seq);
}
