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
 * @file RationalBeattySequence.h
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

#include <boost/rational.hpp>

class RationalBeattySequence {
private:
    boost::rational<int> _ratio;
    int _offset;

public:
    RationalBeattySequence(boost::rational<int> ratio, int offset) :
	_ratio(ratio),
	_offset(offset) { }

    RationalBeattySequence invert() const {
	boost::rational<int> r(_ratio.denominator(), _ratio.numerator());
	return RationalBeattySequence(r, -_offset - 1);
    }

    RationalBeattySequence complement() const {
	boost::rational<int> r(_ratio.numerator(), _ratio.numerator() - _ratio.denominator());
	return RationalBeattySequence(r, -_offset - 1);
    }

    int operator()(int n) const {
	//assert(n >= 0);
	// Floor dir: floor(n*tau) -> (n*num)/den
	// Ceil dir: ceil(n*tau - 1) -> (n*num+den-1)/den - 1 -> (n*num-1)/den
	return (_ratio.numerator() * n + _offset) / _ratio.denominator();
    }

    void print() const {
	printf("⌊(%d*n", _ratio.numerator());
	if (_offset!=0)
	    printf("%+d", _offset);
	printf(")/%d⌋\n", _ratio.denominator());
	//printf("num: %d, den: %d, dir: %d\n", num, den, offset);
    }

    bool equals(RationalBeattySequence &otherSeq) const {
	return _ratio == otherSeq._ratio &&
	       _offset == otherSeq._offset;
    }
};
