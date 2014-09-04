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
 * @file sequenceTest.h
 * @ingroup Tools
 * @author Nicolas Normand (\c Nicolas.Normand@polytech.univ-nantes.fr)
 * LUNAM Université, Université de Nantes, IRCCyN UMR CNRS 6597
 *
 * @date 2012/10/11
 *
 * LUTBasedNSDistanceTransform computes the 2D translated neighborhood-sequence
 * distance transform of a binary image. It reads the input images from its
 * standard input and writes the result to its standard output.
 *
 * This file is part of the DGtal library.
 */

#include <boost/assert.hpp>

/**
 * Checks if two sequences are Lambek-Moser inverse sequences, *i.e.* if the
 * sequences obtained by adding the rank of each term
 * (s1(1)+1,s1(2)+2,...,s1(k)+k,... and s2(1)+1,s2(2)+2,...,s2(k)+k,...) are
 * complementary sequences \cite lambek1954amm.
 *
 * The test uses the extra assumption that the sequences are non decreasing.
 * It has the structure of a merge algorithm for sorted lists.
 * Practically, the test is performed from **n=1** to a finite limit **N**
 * (defaults to 10000).
 */
template<typename sequence>
bool testLambekMoserInverseSequences(sequence &s1, sequence &s2, unsigned int N = 10000) {
    unsigned int i1 = 1;
    unsigned int i2 = 1;
    for (unsigned int n = 1; n <= N; n++) {
	// Each positive integer n is expected to be equal to a sequence term
	// plus its rank.
	if (s1(i1) + i1 == n)
	    i1++;
	else if (s2(i2) + i2 == n)
	    i2++;
	else
	    // ... otherwise the test fails
	    return false;
    }
    // The test succeded (so far...)
    return true;
}

/**
 * Checks if two sequences are complementary sequences, *i.e.* if each positive
 * integer **n** appears once in exactly one of the sequences.
 *
 * The test uses the extra assumption that the sequences are non decreasing.
 * It has the structure of a merge algorithm for sorted lists.
 * Practically, the test is performed from **n=1** to a finite limit **N**
 * (defaults to 10000).
 */
template<typename sequence>
bool testComplementarySequences(sequence &s1, sequence &s2, unsigned int N = 10000) {
    int i1 = 1;
    int i2 = 1;
    for (int n = 1; n <= N; n++) {
	// Each positive integer n is expected to be found in one sequence.
      if (s1(i1) == n)
	    i1++;
	else if (s2(i2) == n)
	    i2++;
	else
	    // ... otherwise the test fails
	    return false;
    }
    // The test succeded (so far...)
    return true;
}
