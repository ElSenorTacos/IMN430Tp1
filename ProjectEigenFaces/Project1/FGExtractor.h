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
	void erosionImage(std::vector<std::pair<int, int>>& interests, ImageType& input);
	void dilatationImage(std::vector<std::pair<int, int>>& interests, ImageType& input);
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
            maxi = center.first + sz / 2; if(maxi >= maxValues.first) maxi = maxValues.first - 1;
            minj = center.second - sz / 2; if(minj < 0) minj = 0;
            maxj = center.second + sz / 2; if(maxj >= maxValues.second) maxj = maxValues.second - 1;
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
		template<class F, class I>
		void affect(I& image, F&& f) {
			for (size_t i = mini; i <= maxi; ++i) {
				for (size_t j = minj; j <= maxj; ++j) {
					image(i,j) = f(image, i, j);
				}
			}
		}
    };
};

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
void FGExtractor<T>::erosionImage(std::vector<std::pair<int, int>>& interests, ImageType& input)
{
	std::vector<std::pair<int, int>> keptList;
	size_t pointListSize = interests.size();

	for (int i = 0; i < pointListSize; ++i)
	{
		int x = interests[i].first;
		int y = interests[i].second;

		if (input(x + 1, y, 0) != 0 &&
			input(x - 1, y, 0) != 0 &&
			input(x, y + 1, 0) != 0 &&
			input(x, y - 1, 0) != 0 &&
			input(x + 1, y + 1, 0) != 0 &&
			input(x - 1, y - 1, 0) != 0 &&
			input(x + 1, y - 1, 0) != 0 &&
			input(x - 1, y + 1, 0) != 0)
		{
			keptList.emplace_back(x, y);
		}
	}
	interests = keptList;
}


template <class T>
void FGExtractor<T>::dilatationImage(std::vector<std::pair<int, int>>& interests, ImageType& input)
{
	std::vector<std::pair<int, int>> keptList;
	size_t pointListSize = interests.size();

	Neighborhood neigborhood;
	std::pair<int, int> limits{ input.width(), input.height() };

	for (int i = 0; i < pointListSize; ++i)
	{
		neigborhood.place(interests[i], 3, limits)->affect(input, [](ImageType&, int, int) { return 255.f; });
	}


	interests.clear();

	cimg_forXY(input, x, y)
	{
		if(input(x, y, 0) != 0)
		{
			interests.emplace_back(x, y);
		}
	}

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
		double variance = neighborhood.place({ x, y }, double(windowSize), limits)->sumApplyNomalized(input, [&](T pixel, size_t i, size_t j) {
			return std::pow(double(pixel) - average(i, j), 2.f);
		}, 0.f);
		if ( abs((input(x, y) - average(x, y)) / variance) > 0.1 && (x > 2 && y >2 ) )
		{
			interests.emplace_back(x, y);
		}
	}

	ImageType pointMap(input.width(), input.height(), input.depth(), input.spectrum());
	pointMap.fill(0);
	for (int i = 0; i <interests.size(); ++i)
	{
		for (int chanel = 0; chanel < input.spectrum(); ++chanel)
		{
			pointMap(interests[i].first, interests[i].second, 0, chanel) = 255;
		}
	}

	dilatationImage(interests, pointMap);
	erosionImage(interests, pointMap);

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
		output(x, y) = neighborhood.place({ x, y }, double(windowSize), limits)->sumApplyNomalized(input, [&](T pixel, size_t i, size_t j) {
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
		output(x, y) = neighborhood.place({ x, y }, windowSize, limits)->sumApplyNomalized(input, [&](T pixel, size_t i, size_t j) {
			return std::pow(double(pixel) - average(i, j), 2.f);
		}, 0.f);
	}
	return std::forward<ImageType>(output);
}

#endif //__FGExtractor_h_