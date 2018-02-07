#ifndef __ImageVector_h_
#define __ImageVector_h_

#include <vector>
#include "Eigen/Cholesky"
#include "CImg.h"

template<class T>
class ImageVector
{
public:
    typedef     Eigen::VectorXd                      VectorizedComponentType;
    typedef     std::vector<VectorizedComponentType>    VectorizedImageType;
    typedef     VectorizedImageType::iterator           VectorizedImageIteratorType;
    typedef     VectorizedImageType::const_iterator     VectorizedImageConstIteratorType;

    typedef     cimg_library::CImg<T>                   ImageType;

	ImageVector() {}
    explicit ImageVector(ImageType image) { vectorize(image); };
    ~ImageVector();

    VectorizedImageConstIteratorType begin() { return imageComponents.begin(); }
    VectorizedImageConstIteratorType end() { return imageComponents.end(); }

    void clear() { imageComponents.clear(); }
    void resize(const size_t size) { clear(); imageComponents.resize(size); }
    void setComponent(const size_t, VectorizedComponentType);
	void setComponentFree(const size_t, VectorizedComponentType);
    VectorizedComponentType getComponent(size_t i){ return imageComponents.at(i); }
    
	void setPixelSize(int size) { for (int i = 0; i < componentsCount(); ++i) { imageComponents[i].resize(size); } }

    size_t componentsCount() const { return imageComponents.size(); }
    size_t pixelCount() const { return imageComponents.front().size(); }

    void save(std::string name, size_t sizeX, size_t sizeY);

private:

    void vectorize(ImageType);
    void initialize(const int, const int);

    VectorizedImageType     imageComponents;
};

template<class T>
ImageVector<T>::~ImageVector()
{
}

template <class T>
void ImageVector<T>::setComponent(const size_t componentIndex, VectorizedComponentType component)
{
	assert(imageComponents.at(componentIndex).size() == component.size());
	imageComponents.at(componentIndex) = component;
}

template <class T>
void ImageVector<T>::setComponentFree(const size_t componentIndex, VectorizedComponentType component)
{
	imageComponents.at(componentIndex) = component;
}

template<class T>
void ImageVector<T>::vectorize(ImageType image)
{
	using namespace cimg_library;

	int width = image.width();

	initialize(image.spectrum(), width * image.height());

	cimg_forXY(image, x, y)
	{
		for (int i = 0; i < image.spectrum(); ++i)
		{
			imageComponents[i][x + y * width] = image(x, y, i);
		}
	}
}

template <class T>
void ImageVector<T>::initialize(const int nbChannels, const int size)
{
	imageComponents.clear();
	imageComponents.resize(nbChannels);
	for(int i = 0; i < nbChannels; ++i)
	{
		imageComponents[i].resize(size);
	}
}

template <class T>
void ImageVector<T>::save(std::string name, size_t sizeX, size_t sizeY)
{
	ImageType output(sizeX, sizeY, 1, imageComponents.size());
	cimg_forXY(output, x, y)
	{
		for (size_t i = 0; i < imageComponents.size(); ++i)
		{
			output(x, y, 1, i) = imageComponents.at(i)(x + y * sizeX);
		}
	}
}

#endif //__ImageVector_h_