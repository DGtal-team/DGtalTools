// $Id: RationalBeattySequence.h 101 2012-07-09 07:41:29Z Nicolas.Normand $

extern "C" {
typedef struct RationalBeattySequence RationalBeattySequence;

RationalBeattySequence *RationalBeattySequenceCreate(int num, int den, int dir);

RationalBeattySequence *RationalBeattySequenceCreateInverse(RationalBeattySequence *seq);

RationalBeattySequence *RationalBeattySequenceCreateComplement(RationalBeattySequence *seq);

void RationalBeattySequenceFree(RationalBeattySequence *seq);

int RationalBeattySequenceValueAtIndex(RationalBeattySequence *seq, int i);

void RationalBeattySequencePrint(RationalBeattySequence *seq);

int RationalBeattySequenceEquals(RationalBeattySequence *seq1, RationalBeattySequence *seq2);
}

class RationalBeattySeq {
private:
    int num;
    int den;
    int offset;

public:
    RationalBeattySeq(int num, int den, int offset) :
    num(num),
    den(den),
    offset(offset) { }

    RationalBeattySeq invert() const {
	return RationalBeattySeq(den, num, -offset - 1);
    }

    RationalBeattySeq complement() const {
	return RationalBeattySeq(num, num - den, -offset - 1);
    }

    int operator()(int n) const {
	//assert(n >= 0);
	// Floor dir: floor(n*tau) -> (n*num)/den
	// Ceil dir: ceil(n*tau - 1) -> (n*num+den-1)/den - 1 -> (n*num-1)/den
	return (num * n + offset) / den;
    }

    void print() const {
	printf("⌊(%d*n",num);
	if (offset!=0)
	    printf("%+d", offset);
	printf(")/%d⌋\n", den);
	//printf("num: %d, den: %d, dir: %d\n", num, den, offset);
    }

    bool equals(RationalBeattySeq &otherSeq) const {
	return num == otherSeq.num &&
	       den == otherSeq.den &&
	       offset == otherSeq.offset;
    }
};
