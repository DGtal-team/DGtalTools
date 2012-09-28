// $Id: CumulativeSequence.h 94 2012-07-04 07:32:53Z Nicolas.Normand $

#include <vector>

class CumulativeOfPeriodicSequence {
public:
    CumulativeOfPeriodicSequence(int length, int offset = 0) :
	_sequence(length),
	_offset(offset) {
    }

    CumulativeOfPeriodicSequence(std::vector<int> sequence, int offset = 0) :
        _sequence(sequence),
	_offset(offset)
    {
	int sum = 0;
	for (std::vector<int>::iterator it = _sequence.begin(); 
	     it != _sequence.end();
	     it++)
	{
	    sum += *it;
	    *it = sum;
	}
    }

    CumulativeOfPeriodicSequence *inverse();

    int valueAtIndex(int i);

    int equals(CumulativeOfPeriodicSequence& seq);

    void print();
private:
    std::vector<int> _sequence;
    int _offset;
};
