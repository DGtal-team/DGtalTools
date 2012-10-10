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
 * @file PNGImageReader.cpp
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

#include "PNGImageReader.h"

PNGImageReader::PNGImageReader(ImageConsumer<BinaryPixelType>* consumer, FILE *input) :
RowImageProducer<BinaryPixelType>(consumer),
_input(input) {
}

void
PNGImageReader::startImage(size_t readBytes)
{
    _png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!_png_ptr) {
	return;
    }

    _info_ptr = png_create_info_struct(_png_ptr);
    if (!_info_ptr) {
	png_destroy_read_struct(&_png_ptr, (png_infopp)NULL, (png_infopp)NULL);
	return;
    }

    if (setjmp(png_jmpbuf(_png_ptr))) {
	png_destroy_read_struct(&_png_ptr, &_info_ptr, (png_infopp)NULL);
	std::cerr << "Error during PNGImageReader::PNGImageReader" << std::endl;
	exit(1);
    }

    png_init_io(_png_ptr, _input);

    // Tell libpng we've already read 'readBytes' bytes
    png_set_sig_bytes (_png_ptr, readBytes);

    png_read_info(_png_ptr, _info_ptr);
    int bit_depth, color_type, interlace_type;
    png_get_IHDR(_png_ptr, _info_ptr, &_cols, &_rows, &bit_depth, &color_type,
		 &interlace_type, NULL, NULL);

    /* Convert to grayscale */
    if (color_type == PNG_COLOR_TYPE_RGB ||
        color_type == PNG_COLOR_TYPE_RGB_ALPHA)
	png_set_rgb_to_gray_fixed(_png_ptr, 1 /* error_action */,
				  -1 /* int red_weight */, -1 /* int green_weight */);

    /* Drop alpha channel */
    if (color_type & PNG_COLOR_MASK_ALPHA)
        png_set_strip_alpha(_png_ptr);

    /* Downscale 16 bits/pixel grayscale images to 8 bits/pixel */
    if (bit_depth == 16)
        png_set_strip_16(_png_ptr);

    /* Expand grayscale images to the full 8 bits from 1, 2, or 4 bits/pixel */
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
	png_set_expand_gray_1_2_4_to_8(_png_ptr);

    //png_set_invert_mono(_png_ptr);

    png_read_update_info(_png_ptr, _info_ptr);
    png_get_IHDR(_png_ptr, _info_ptr, &_cols, &_rows, &bit_depth, &color_type,
		 &interlace_type, NULL, NULL);

    // Now, everything should be grayscale, 8 bit/pixel
    assert(bit_depth == 8);
    assert(color_type == PNG_COLOR_TYPE_GRAY);

    //int number_of_passes = png_set_interlace_handling(_png_ptr);

    _inputRow = (BinaryPixelType *)malloc(png_get_rowbytes(_png_ptr, _info_ptr));
    _consumer->beginOfImage(_cols, _rows);
}

void PNGImageReader::produceAllRows(size_t readBytes) {
    startImage(readBytes);
    if (setjmp(png_jmpbuf(_png_ptr))) {
	png_destroy_read_struct(&_png_ptr, &_info_ptr, (png_infopp)NULL);
	std::cerr << "Error during PNGImageReader::produceAllRows" << std::endl;
	exit(1);
    }

    while (_rows-- > 0) {
	png_read_row(_png_ptr, (png_bytep)_inputRow, NULL);
	for (unsigned int col = 0; col < _cols; col++) {
	    _inputRow[col] = _inputRow[col] > 127 ? BINARY_WHITE_PIXEL : BINARY_BLACK_PIXEL;
	}
	_consumer->processRow(_inputRow);
    }
    _consumer->endOfImage();
    free(_inputRow);
    _inputRow = NULL;
    png_read_end(_png_ptr, NULL);
}
