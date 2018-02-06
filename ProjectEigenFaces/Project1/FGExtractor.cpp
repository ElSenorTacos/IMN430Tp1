#include "FGExtractor.h"


template <class T>
inline FGExtractor<T>::FGExtractor()
{
}

template <class T>
inline FGExtractor<T>::~FGExtractor()
{
}

template <class T>
typename FGExtractor<T>::ImageType FGExtractor<T>::estimageForeground(ImageType image)
{
    //ImageType
}

template <class T>
std::vector<std::pair<int, int>> FGExtractor<T>::getForegroundPixelsPositions(ImageType input, int windowSize)
{
    cimg_library::CImg<double> average(computeAverage(input, windowSize));
    std::vector<std::pair<int, int>> interests;
    std::pair<int, int> limits{ input.width(), input.height() };
    Neighborhood neighborhood;
    cimg_forXY(input, x, y)
    {
        double variance = neighborhood.place({ x, y }, double(windowSize), limits)->sumApplyNomalized(input, [&](T pixel, size_t i, size_t j){
            return std::pow(double(pixel) - average(i,j), 2.f);
        }, 0.f);
        if (abs((input(x,y) - average(x,y)) / variance) > 0.5)
        {
            interests.emplace_back( x, y );
        }
    }
    return std::forward<std::vector<std::pair<int, int>>>(interests);
}

template<class T>
inline typename FGExtractor<T>::ImageType FGExtractor<T>::computeAverage(ImageType& input, int windowSize)
{
    cimg_library::CImg<double> output(input.width(), input.height(), input.depth(), input.spectrum());
    std::pair<int, int> limits{ input.width(), input.height() };
    Neighborhood neighborhood;
    cimg_forXY(input, x, y)
    {
        output(x,y) = neighborhood.place({ x, y }, double(windowSize), limits)->sumApplyNomalized(input, [&](T pixel, size_t i, size_t j){
            return double(pixel);
        }, 0.f);
    }
    return std::forward<ImageType>(output);
}

template<class T>
inline typename FGExtractor<T>::ImageType FGExtractor<T>::computeVariance(ImageType& input, int windowSize)
{
    cimg_library::CImg<double> output(input.width(), input.height(), input.depth(), input.spectrum()),
                               average = computeAverage(input, windowSize);
    std::pair<int, int> limits{ input.width(), input.height() };
    Neighborhood neighborhood;
    cimg_forXY(input, x, y)
    {
        output(x,y) = neighborhood.place({ x, y }, windowSize, limits)->sumApplyNomalized(input, [&](T pixel, size_t i, size_t j){
            return std::pow(double(pixel) - average(i,j), 2.f);
        }, 0.f);
    }
    return std::forward<ImageType>(output);
}
