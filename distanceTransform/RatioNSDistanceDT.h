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

#include "NeighborhoodSequenceDistance.h"
#include "RationalBeattySequence.h"

/**
 * \brief This class represents neighborhood sequence distances defined by
 * a ratio of neighborhoods. It provides factory methods to create translated
 * distance transform filters and distance transform untranslator filters.
 */
class RatioNSDistance: public NeighborhoodSequenceDistance {
public:
    /**
     * Constructor.  Creates a neighborhood sequence distance defined by a
     * a ratio of neighborhoods.
     *
     * \param ratio ratio of appearance of the 2-neighborhood in the sequence.
     * **ratio** should be greater than 0 and less than 1. For instance,
     * **ratio=1/2** creates the octagonal distance, *i.e.* the neighborhood
     * sequence distance with the periodic sequence of neighborhoods 1,2,1,2...
     * (strict alternance of neighborhoods).
     * The sequence of neighborhoods is the first difference of a rational
     * Beatty sequence with parameter **ratio+1**. In the octagonal case, the
     * sequence is B(r)=⌊3r/2⌋-⌊3(r-1)/2⌋.
     .
     * \latexonly B(r)=\lfloor3r/2\rfloor\endlatexonly
     */
    RatioNSDistance(boost::rational<int> ratio);

    NeighborhoodSequenceDistanceTransform* newTranslatedDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) const;

    DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType>* newDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer) const;

    friend class RatioNSDistanceTransform;
    friend class RatioNSDistanceTransformUntranslator;

protected:
    /**
     * ratio of appearance of the 2-neighborhood in the sequence.
     */
    boost::rational<int> _ratio;
    /**
     * counts the occurrences of neighorhood 1 in the sequence of neighborhoods.
     * #mbf1 is a Rational Beatty sequence with parameter 1 - #_ratio.
     */
    RationalBeattySequence mbf1;
    /**
     * counts the occurrences of neighorhood 2 in the sequence of neighborhoods.
     * #mbf2 is a Rational Beatty sequence with parameter #_ratio.
     */
    RationalBeattySequence mbf2;
    //! Lambek-Moser inverse of #mbf1.
    RationalBeattySequence mbf1i;
    //! Lambek-Moser inverse of #mbf2.
    RationalBeattySequence mbf2i;
};

/**
 * \brief Implements a single scan translated distance transform for distances
 * defined by a ratio of neighborhoods.
 */
class RatioNSDistanceTransform : public NeighborhoodSequenceDistanceTransform {
protected:
    const RatioNSDistance d;

public:
    RatioNSDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer, boost::rational<int> ratio);

    void processRow(const BinaryPixelType *imageRow);
};

/**
 * \brief Implements a recentering algorithm for the translated distance
 * transforms defined by a ratio of neighborhoods.
 */
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
