// $Id: RatioNSDistanceDT.h 106 2012-07-10 20:40:36Z Nicolas.Normand $

#include "BaseDistanceDT.h"
#include "RationalBeattySequence.h"

class RatioNSDistance: public BaseDistance {
public:
    RatioNSDistance(int num, int den);
    BaseDistanceTransform* newTranslatedDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) const;
    DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType>* newDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer) const;

    int num;
    int den;
    RationalBeattySeq mbf1;
    RationalBeattySeq mbf2;
    RationalBeattySeq mbf1i;
    RationalBeattySeq mbf2i;
};

class RatioNSDistanceTransform : public BaseDistanceTransform {
protected:
    const RatioNSDistance d;

public:
    RatioNSDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer, int num, int den);

    void processRow(const BinaryPixelType *imageRow);
};

class RatioNSDistanceTransformUntranslator: public DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType> {
public:
    RatioNSDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer, int dMax, int num, int den);
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
