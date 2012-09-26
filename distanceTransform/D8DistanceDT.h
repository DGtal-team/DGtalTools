// $Id: D8DistanceDT.h 106 2012-07-10 20:40:36Z Nicolas.Normand $

#include "BaseDistanceDT.h"

class D8Distance: public BaseDistance {
public:
    BaseDistanceTransform* newTranslatedDistanceTransform(ImageConsumer<GrayscalePixelType>* consumer) const;
    DistanceTransformUntranslator<GrayscalePixelType, GrayscalePixelType>* newDistanceTransformUntranslator(ImageConsumer<GrayscalePixelType>* consumer) const;
};

class D8DistanceTransform: public BaseDistanceTransform {
public:
    D8DistanceTransform(ImageConsumer<GrayscalePixelType>* consumer);
    
    void processRow(const BinaryPixelType *imageRow);
};

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
