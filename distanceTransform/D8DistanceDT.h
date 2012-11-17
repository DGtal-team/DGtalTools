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
 * @file D8DistanceDT.h
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

#include "NeighborhoodSequenceDistance.h"

/**
 * \brief This class represents the chessboard distance. It provides
 * factory methods to create translated distance transform filters and
 * distance transform untranslator filters.
 */
class D8Distance: public NeighborhoodSequenceDistance {
public:
    NeighborhoodSequenceDistanceTransform* newTranslatedDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) const;
    DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType>* newDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer) const;
};

/**
 * \brief Implements a single scan translated chessboard distance transform.
 */
class D8DistanceTransform: public NeighborhoodSequenceDistanceTransform {
public:
    D8DistanceTransform(ImageConsumer<GrayscalePixelType>* consumer);

    void processRow(const BinaryPixelType *imageRow);
};

/**
 * \brief Implements a recentering algorithm for the translated chessboard
 * distance transform.
 */
class D8DistanceTransformUntranslator: public DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType> {
private:
    typedef DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType> super;
public:
    D8DistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer, int dMax);
    ~D8DistanceTransformUntranslator();

    void beginOfImage(int cols, int rows);
    void processRow(const GrayscalePixelType* inputRow);
    void endOfImage();

protected:
    int _dMax;
    int _imageDMax;
};
