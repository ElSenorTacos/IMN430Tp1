#ifndef __ImageVector_h_
#define __ImageVector_h_

#include <vector>
#include "Eigen/Cholesky"
#include "CImg.h"

template<class T>
class ImageVector
{
public:
    typedef     Eigen::VectorXd                         VectorizedComponentType;
    typedef     std::vector<VectorizedComponentType>    VectorizedImageType;
    typedef     VectorizedImageType::iterator           VectorizedImageIteratorType;
    typedef     VectorizedImageType::const_iterator     VectorizedImageConstIteratorType;

    typedef     cimg_library::CImg<T>                   ImageType;

    explicit ImageVector(ImageType image) { vectorize(image); };
    ~ImageVector();

    VectorizedImageConstIteratorType begin() { return imageComponents.begin(); }
    VectorizedImageConstIteratorType end() { return imageComponents.end(); }

    void clear() { imageComponents.clear(); }
    void resize(const size_t size) { clear(); imageComponents.resize(size); }
    void setComponent(const size_t, VectorizedComponentType);
    VectorizedComponentType getComponent(size_t i){ return imageComponents.at(i); }
    
    void setPixelSize(int size) { std::for_each(begin(), end(), [&](VectorizedComponentType& component) { component.resize(size); }); }

    size_t componentsCount() const { return imageComponents.size(); }
    size_t pixelCount() const { return imageComponents.front().size(); }

    void save(std::string name, size_t sizeX, size_t sizeY);

private:

    void vectorize(ImageType);
    void initialize(const int, const int);

    VectorizedImageType     imageComponents;
};

#endif //__ImageVector_h_