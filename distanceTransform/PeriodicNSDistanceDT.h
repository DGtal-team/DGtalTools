// $Id: PeriodicNSDistanceDT.h 106 2012-07-10 20:40:36Z Nicolas.Normand $

#include "BaseDistanceDT.h"

class PeriodicNSDistance: public BaseDistance {
public:
    PeriodicNSDistance(int period, int *Bvalues);
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
