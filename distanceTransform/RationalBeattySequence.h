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

#include "DGtal/io/Color.h"

#include <boost/rational.hpp>

/**
 * \brief this class represents rational Beatty sequences.
 */
class RationalBeattySequence {
private:
    boost::rational<int> _ratio;
    int _offset;

public:
    /**
     * Constructor. Creates a rational Beatty sequence with parameter **ratio**
     * and offset **offset**.
     */
    RationalBeattySequence(boost::rational<int> ratio, int offset = 0) :
  _ratio(ratio),
  _offset(offset) { }

    /**
     * Creates the Lambek-Moser inverse of the rational Beatty sequence, *i.e.*
     * the sequence such that:
     * ∀(n,m) ∈ ℕ², this(n) < m or inverse(n) < m but not both \cite lambek1954amm.
     *
     * \throws DGtal::InfiniteNumberException if the sequence isn't a
     * non-decreasing integer sequence and doesn't have an inverse
     * (*i.e.* if #_ratio ≤ 0).
     */
    RationalBeattySequence invert() const {
  DGtal::InfiniteNumberException dgtalinfty;
  if (_ratio <= 0) {
      throw dgtalinfty;
  }
  boost::rational<int> r(_ratio.denominator(), _ratio.numerator());
  return RationalBeattySequence(r, -_offset - 1);
    }

    /**
     * Creates the complementary sequence of this rational Beatty sequence. Each
     * integer is produced once by the initial or the complementary sequence.
     *
     * \throws DGtal::InfiniteNumberException if the sequence hasn't a
     * complementary sequence (*i.e.* if #_ratio ≤ 1).
     */
    RationalBeattySequence complement() const {
  DGtal::InfiniteNumberException dgtalinfty;
  if (_ratio <= 1) {
      throw dgtalinfty;
  }
  boost::rational<int> r(_ratio.numerator(), _ratio.numerator() - _ratio.denominator());
  return RationalBeattySequence(r, -_offset - 1);
    }

    /**
     * Computes the value of the sequence for index **n**.
     */
    int operator()(unsigned int n) const {
  //assert(n >= 0);
  // Floor dir: floor(n*tau) -> (n*num)/den
  // Ceil dir: ceil(n*tau - 1) -> (n*num+den-1)/den - 1 -> (n*num-1)/den
  return (_ratio.numerator() * n + _offset) / _ratio.denominator();
    }

    friend std::ostream &operator<<(std::ostream &out, const RationalBeattySequence &seq) {
  out << "⌊";
  if (seq._offset != 0)
      out << '(';
  if (seq._ratio.numerator() == -1)
      out << '-';
  else if (seq._ratio.numerator() != 1)
      out << seq._ratio.numerator();
  out << 'n';
  if (seq._offset > 0)
      out << '+';
  if (seq._offset != 0)
      out << seq._offset << ')';
  out << '/' << seq._ratio.denominator() << "⌋";
  return out;
    }

    bool operator==(RationalBeattySequence &otherSeq) const {
  return _ratio == otherSeq._ratio &&
  _offset == otherSeq._offset;
    }
};
