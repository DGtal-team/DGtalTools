// $Id: D4DistanceDT.h 106 2012-07-10 20:40:36Z Nicolas.Normand $

#include "BaseDistanceDT.h"

class D4Distance: public BaseDistance {
public:
    BaseDistanceTransform* newTranslatedDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) const;
    DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType>* newDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer) const;
};

class D4DistanceTransform: public BaseDistanceTransform {
public:
    D4DistanceTransform(ImageConsumer<GrayscalePixelType>* consumer);

    void processRow(const BinaryPixelType *imageRow);

private:
    typedef BaseDistanceTransform super;
};

class D4DistanceTransformUntranslator: public DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType> {
private:
    typedef DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType> super;
public:
    D4DistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer, int dMax);
    ~D4DistanceTransformUntranslator();

    void beginOfImage(int cols, int rows);
    void processRow(const GrayscalePixelType* inputRow);
    void endOfImage();
    
protected:
    int _dMax;
    int _imageDMax;
};
