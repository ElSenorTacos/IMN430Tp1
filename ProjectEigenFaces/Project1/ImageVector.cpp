#include "ImageVector.h"




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
    std::for_each(begin(), end(), [](VectorizedImageIteratorType componentVector)
    {
        componentVector->resize(size);
    });
}
