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
 * @file RationalBeattySequenceTest.cpp
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
#include "RationalBeattySequence.h"

using namespace std;

void testSequence(boost::rational<int> ratio, int dir) {
    int i;

    RationalBeattySequence bs(ratio, dir);
    cout << "Sequence: " << bs << endl;
    for (i = 1; i < 15; i++) {
	cout << bs(i) << ' ';
    }
    cout << endl << endl;

    RationalBeattySequence bsi = bs.invert();
    BOOST_VERIFY(testLambekMoserInverseSequences(bs, bsi));
    cout << "Inverse of sequence: " << bsi << endl;
    for (i = 1; i < 15; i++) {
	cout << bsi(i) << ' ';
    }
    cout << endl << endl;

    RationalBeattySequence bsii = bsi.invert();
    BOOST_VERIFY(testLambekMoserInverseSequences(bsi, bsii));
    cout << "Inverse of inverse of sequence: " << bsii << endl;
    for (i = 1; i < 15; i++) {
	cout << bsii(i) << ' ';
    }
    cout << endl << endl;

    RationalBeattySequence bsc = bs.complement();
    BOOST_VERIFY(testComplementarySequences(bs, bsc));
    cout << "Complement of sequence: " << bsc << endl;
    for (i = 1; i < 15; i++) {
	cout << bsc(i) << ' ';
    }
    cout << endl << endl;

    RationalBeattySequence bscc = bsc.complement();
    BOOST_VERIFY(testComplementarySequences(bsc, bscc));
    cout << "Complement of complement of sequence: " << bscc << endl;
    for (i = 1; i < 15; i++) {
	cout << bscc(i) << ' ';
    }
    cout << endl << endl;

    BOOST_VERIFY(bs == bsii);
    BOOST_VERIFY(bs == bscc);
}

int main(int argc, char** argv) {
    try {
	cout << "Rational Beatty sequence with rate 3/2" << endl;
	cout << "Test should proceed without errors" << endl;
	testSequence(boost::rational<int>(3, 2), 0);
    }
    catch (exception& e) {
	cout << e.what() << endl;
    }

    cout << endl << endl;

    try {
	cout << "Rational Beatty sequence with rate 1/2" << endl;
	cout << "No complementary sequence (rate <= 1), should throw an exception" << endl;
	testSequence(boost::rational<int>(1, 2), 0);
    }
    catch (exception& e) {
	cout << e.what() << endl;
    }

    cout << endl << endl;

    try {
	cout << "Rational Beatty sequence with rate 1" << endl;
	cout << "No complementary sequence (rate <= 1), should throw an exception" << endl;
	testSequence(boost::rational<int>(1, 1), 0);
    }
    catch (exception& e) {
	cout << e.what() << endl;
    }

    try {
	cout << "Rational Beatty sequence with rate 0" << endl;
	cout << "No inverse sequence (rate <= 0), should throw an exception" << endl;
	testSequence(boost::rational<int>(0, 1), 0);
    }
    catch (exception& e) {
	cout << e.what() << endl;
    }

    return 0;
}
