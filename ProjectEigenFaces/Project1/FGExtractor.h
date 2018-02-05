#ifndef __FGExtractor_h_
#define __FGExtractor_h_
#include "CImg.h"
#include <vector>

template<class T>
class FGExtractor
{
public:
    typedef cimg_library::CImg<T>   ImageType;

    ImageType estimageForeground(ImageType image);
    std::vector<std::pair<int, int>> getForegroundPixelsPositions(ImageType image);

    FGExtractor();
    ~FGExtractor();
private:
    ImageType computeAverage();
};

#endif //__FGExtractor_h_
