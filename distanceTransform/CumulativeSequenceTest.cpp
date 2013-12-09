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
 * @file CumulativeSequenceTest.cpp
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

#include <iostream>

#include "sequenceTest.h"
#include "CumulativeSequence.h"

using namespace std;

void testSequence(vector<int> sequence, int offset) {
    int i;

    cout << "Sequence: ";
    CumulativeOfPeriodicSequence cs(sequence, offset);
    for (i = 1; i < 15; i++) {
	cout << cs(i) << ' ';
    }
    cout << cs << endl;

    cout << "Inverse of sequence: ";
    CumulativeOfPeriodicSequence csi(cs.invert());
    BOOST_VERIFY(testLambekMoserInverseSequences(cs, csi));
    for (i = 1; i < 15; i++) {
	cout << csi(i) << ' ';
    }
    cout << csi << endl;

    cout << "Inverse of inverse of sequence: ";
    CumulativeOfPeriodicSequence csii(csi.invert());
    BOOST_VERIFY(testLambekMoserInverseSequences(csi, csii));
    for (i = 1; i < 15; i++) {
	cout << csii(i) << ' ';
    }
    cout << csii << endl;

    BOOST_VERIFY(csii == cs);
}

int main(int argc, char** argv) {
    //int period;
    //int seq2[4] = {2, 0, 0, 1};
    vector<int> seq2;
    seq2.push_back(2);
    seq2.push_back(0);
    seq2.push_back(0);
    seq2.push_back(1);

    testSequence(seq2, -4);

    testSequence(seq2, -1);
    testSequence(seq2,  0);
    testSequence(seq2,  1);
    testSequence(seq2,  10);
    testSequence(seq2,  100);
    testSequence(seq2, -10);
    testSequence(seq2, -100);

    //int sequence[4] = {1, 2, 0, 3};	// cf. Table 1
    vector<int> sequence;
    sequence.push_back(1);
    sequence.push_back(2);
    sequence.push_back(0);
    sequence.push_back(3);
    
    testSequence(sequence, 0);

    return 0;
}
