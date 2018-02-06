#ifndef __FGExtractor_h_
#define __FGExtractor_h_
#include "CImg.h"
#include <vector>

template<class T>
class FGExtractor
{
public:
    typedef cimg_library::CImg<T>   ImageType;

    ImageType estimageForeground(ImageType);
    std::vector<std::pair<int, int>> getForegroundPixelsPositions(ImageType, int);

    FGExtractor();
    ~FGExtractor();
    ImageType computeAverage(ImageType&, int);
    ImageType computeVariance(ImageType&, int);
private:


    struct Neighborhood
    {
        size_t mini, maxi, minj, maxj, size;
        Neighborhood() {};
        Neighborhood* place(std::pair<int,int> center, int size, std::pair<int,int> maxValues)
        {
            size_t offset = size % 2,
                   sz = size - offset;
            mini = center.first - sz / 2; if(mini < 0) mini = 0;
            maxi = center.first + sz / 2; if(maxi > maxValues.first) maxi = maxValues.first;
            minj = center.second - sz / 2; if(minj < 0) minj = 0;
            maxj = center.second + sz / 2; if(maxj > maxValues.second) maxj = maxValues.second;
            this->size = (maxj - minj + 1) * (maxi - mini + 1);
            return this;
        }
        template<class R, class F, class I>
        R sumApply(I& image, F&& f, R initial_value) {
            R bin = initial_value;
            for(size_t i = mini; i <= maxi; ++i) {
                for(size_t j = minj; j <= maxj; ++j) {
                    bin += f(image(i,j), i, j);
                }
            }
            return bin;
        }
        template<class R, class F, class I>
        R sumApplyNomalized(I& image, F&& f, R initial_value) {
            R val = sumApply(image, std::forward<F>(f), initial_value) / R(size);
            return val;
        }
    };
};

#endif //__FGExtractor_h_

#if not FGExtractor_Manual_Instanciation
#define FGExtractor_Manual_Instanciation 1
#include "FGExtractor.cpp"
#endif