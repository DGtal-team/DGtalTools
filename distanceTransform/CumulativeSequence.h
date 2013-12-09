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
 * @file CumulativeSequence.h
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

#include <vector>
#include <ostream>

/**
 * \brief This class represents non-decreasing sequences of integers produced
 * by the cumulative sum of periodic sequences.
 */
class CumulativeOfPeriodicSequence {
public:
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

    CumulativeOfPeriodicSequence invert() const;

    int operator()(int i) const;

    bool operator==(const CumulativeOfPeriodicSequence& seq) const;

    friend std::ostream &operator<<(std::ostream &out, const CumulativeOfPeriodicSequence &seq);

protected:
    CumulativeOfPeriodicSequence(int length, int offset = 0) :
    _sequence(length),
    _offset(offset) {
    }

private:
    std::vector<int> _sequence;
    int _offset;
};
