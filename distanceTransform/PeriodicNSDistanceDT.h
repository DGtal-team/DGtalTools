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
 * @file PeriodicNSDistanceDT.h
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

#include <vector>

class PeriodicNSDistance: public BaseDistance {
public:
    PeriodicNSDistance(std::vector<int> Bvalues);
    ~PeriodicNSDistance();
    BaseDistanceTransform* newTranslatedDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) const;
    DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType>* newDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer) const;

    int mathbf2(int r) const;
    GrayscalePixelType C1(int r) const;
    GrayscalePixelType C2(int r) const;

    int period;
    int *c1;

protected:
    int *c2;
    int *mathbf1d;
    int *mathbf2d;
};

class PeriodicNSDistanceTransform : public BaseDistanceTransform {
public:
    PeriodicNSDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer, const PeriodicNSDistance *d);
    ~PeriodicNSDistanceTransform();

    void processRow(const BinaryPixelType *imageRow);
    void untranslate(int cols);
protected:
    const PeriodicNSDistance *_d;
};

class PeriodicNSDistanceTransformUntranslator: public DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType> {
public:
    PeriodicNSDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer, int dMax, const PeriodicNSDistance *d);
    ~PeriodicNSDistanceTransformUntranslator();

    void beginOfImage(int cols, int rows);
    void processRow(const GrayscalePixelType* inputRow);
    void endOfImage();

protected:
    typedef DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType> super;
    void rotate();

    static const int marginRight = 1;
    const int _dMax;
    const PeriodicNSDistance *_d;

    int _imageDMax;
};
