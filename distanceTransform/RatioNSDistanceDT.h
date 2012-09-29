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
 * @file RatioNSDistanceDT.h
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

#include "BaseDistanceDT.h"
#include "RationalBeattySequence.h"

class RatioNSDistance: public BaseDistance {
public:
    RatioNSDistance(boost::rational<int> ratio);
    BaseDistanceTransform* newTranslatedDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) const;
    DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType>* newDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer) const;

    boost::rational<int> _ratio;
    RationalBeattySequence mbf1;
    RationalBeattySequence mbf2;
    RationalBeattySequence mbf1i;
    RationalBeattySequence mbf2i;
};

class RatioNSDistanceTransform : public BaseDistanceTransform {
protected:
    const RatioNSDistance d;

public:
    RatioNSDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer, boost::rational<int> ratio);

    void processRow(const BinaryPixelType *imageRow);
};

class RatioNSDistanceTransformUntranslator: public DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType> {
public:
    RatioNSDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer, int dMax, boost::rational<int> ratio);
    ~RatioNSDistanceTransformUntranslator();

    void beginOfImage(int cols, int rows);
    void processRow(const GrayscalePixelType* inputRow);
    void endOfImage();

protected:
    typedef DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType> super;
    void rotate();

    static const int marginRight = 1;
    int _dMax;
    int _imageDMax;
    const RatioNSDistance d;
};
