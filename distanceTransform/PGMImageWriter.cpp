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
 * @file PGMImageWriter.cpp
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

#include "PGMImageWriter.h"

#include <string.h>

PGMImageWriter::PGMImageWriter(FILE* output, int format) :
_cols(0),
_format(format),
_output(output) {
    _outpam.size        = sizeof(struct pam);
    _outpam.len         = sizeof(struct pam);
    _outpam.file        = output;
    _outpam.format      = PGM_FORMAT;
    _outpam.plainformat = 1;
}

void
PGMImageWriter::beginOfImage(int cols, int rows) {
    _cols = cols;
    _outpam.height	         = rows;
    _outpam.width	         = cols;
    _outpam.depth            = 1;
    _outpam.maxval           = 255;
    _outpam.bytes_per_sample = 1;
    strncpy(_outpam.tuple_type, PAM_PGM_TUPLETYPE, sizeof(_outpam.tuple_type));
#ifdef PAM_HAVE_ALLOCATION_DEPTH
    _outpam.allocation_depth = sizeof(GrayscalePixelType);
#endif
#ifdef PAM_HAVE_COMMENT_P
    _outpam.comment_p        = NULL;
#endif
    _tuplerow = pnm_allocpamrow(&_outpam);
    pnm_writepaminit(&_outpam);
    //pgm_writepgminit(_output, cols, rows, 255, _format);
}

void
PGMImageWriter::endOfImage() {
    pnm_freepamrow(_tuplerow);
    _tuplerow = NULL;
}

void
PGMImageWriter::processRow(const GrayscalePixelType* inputRow) {
    // FIXME: type mismatch, change to pgm
    for (int column = 0; column < _outpam.width; ++column) {
	_tuplerow[column][0] = inputRow[column];
    }
    pnm_writepamrow(&_outpam, _tuplerow);
    //pgm_writepgmrow(_output, inputRow, _cols, 255, _format);
}
