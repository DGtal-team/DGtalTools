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

#include <assert.h>
#include <algorithm>

#include "CumulativeSequence.h"

bool CumulativeOfPeriodicSequence::operator ==(const CumulativeOfPeriodicSequence& seq2) const {
    return this->_sequence == seq2._sequence;
}

int mod(int a, int b) {
    return ((a % b) + b) % b;
}

CumulativeOfPeriodicSequence CumulativeOfPeriodicSequence::invert() const {
    //CumulativeOfPeriodicSequence *inv = CumulativeOfPeriodicSequenceCreate(_sequence[seq->period - 1], 0, NULL);
    CumulativeOfPeriodicSequence inv(_sequence[_sequence.size() - 1]);
    inv._offset = 0;

    int xx, yy = 0;
    yy = mod(_offset, inv._sequence.size());
    //yy = CumulativeOfPeriodicSequenceValueAtIndex(seq, 1);
    //yy %= inv._sequence.size();
    for (xx = 0; xx < (int) _sequence.size(); xx++) {
	//yy = CumulativeOfPeriodicSequenceValueAtIndex(seq, xx + 1);
	//yy %= inv._sequence.size();
	inv._sequence[yy]++;
	yy += _sequence[xx] - (xx > 0 ? _sequence[xx-1] : 0);
	yy %= inv._sequence.size();
    }

    for (xx = 1; xx < (int) inv._sequence.size(); xx++) {
	inv._sequence[xx] += inv._sequence[xx-1];
    }
    
    if ((*this)(1) > 0) {
	yy = (*this)(1);
	// Find first increasing index
	xx = 1;
	while ((*this)(xx) == yy) xx++;
	// First positive term in inverse must be equal to xx - 1
	// i.e., inv.value[yy % inv._sequence.size()] + (yy / inv._sequence.size()) * inv._sequence[inv._sequence.size() - 1] + inv._offset == xx - 1
	inv._offset = xx - 1 - inv._sequence[yy % inv._sequence.size()] - (yy / inv._sequence.size()) * inv._sequence[inv._sequence.size() - 1];
	//TODO: assert(CumulativeOfPeriodicSequenceValueAtIndex(inv, yy) == 0);
	//TODO: assert(CumulativeOfPeriodicSequenceValueAtIndex(inv, yy+1) == xx - 1);
    }
    else {
	// Find first positive value
	xx = 0;
	while ((*this)(xx) == 0) xx++;
	// seq(xx) > 0 and seq(xx-1) <= 0 then seqinv(1) = xx-1;
	// i.e., inv.value[0] + inv._offset == xx-1;
	
	inv._offset = xx - 1 - inv._sequence[0];
    }

    return inv;
}

int CumulativeOfPeriodicSequence::operator() (int i) const {
    assert(i >= 0);
    if (i == 0) return 0;
    i--;
    return std::max(_sequence[(i % _sequence.size())] +
		    (long int) (i / _sequence.size()) * _sequence[_sequence.size() - 1] +
		    _offset, 0L);
}

std::ostream &operator<<(std::ostream &out, const CumulativeOfPeriodicSequence &seq) {
  
    out << '(';

    std::vector<int>::const_iterator it = seq._sequence.begin();
    out << *it;
    int prev = *it;

    for (it++; it != seq._sequence.end(); it++) {
	out << ',' << *it - prev;
	prev = *it;
    }

    out <<  ')';

    if (seq._offset > 0) {
	out << "(+" <<  seq._offset << ')';
    }
    else if (seq._offset < 0) {
	out << '(' << seq._offset << ')';
    }

    return out;
}
