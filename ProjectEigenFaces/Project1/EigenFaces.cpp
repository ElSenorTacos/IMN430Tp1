#include "EigenFaces.h"

template<class T>
EigenFaces<T>::EigenFaces(std::string folder, std::string file)
    : parser{ folder, file }
{
}

template<class T>
EigenFaces<T>::~EigenFaces()
{
}

template<class T>
EigenFaces<T>* EigenFaces<T>::apply(const size_t nbComponents)
{
    LinImgPackType vectorizedImages;

    this->image0 = parser.next();
    vectorizedImages.emplace_back({ this->image0 });

    std::for_each(parser.begin() + 1, parser.end(),[]()
    {
        vectorizedImages.emplace_back({ this->realign(parser.next()) });
    });

    this->pca(vectorizedImages, nbComponents);

    return this;
}

template<class T>
typename EigenFaces<T>::ImageType EigenFaces<T>::reconstruct(std::string fileName)
{
    return ImageType();
}

template <class T>
typename EigenFaces<T>::ImageType EigenFaces<T>::realign(ImageType image)
{
    return LinImgPackType();
}

template<class T>
void EigenFaces<T>::pca(LinImgPackType dataVectors, size_t nbComponents)
{
    this->eigenVecImages.resize(dataVectors.size());
    std::for_each(this->eigenVecImages.begin(), this->eigenVecImages.end(), [&](LinImgPackType& vectorImage)
    {
        vectorImage.resize(dataVectors.front().componentsCount());
        vectorImage.setPixelSize(nbComponents);
    });

    Eigen::MatrixXd pcaMatrix(dataVectors.size(), dataVectors.front().pixelCount());
    int imageIndex = 0;
    for (int component = 0; component < dataVectors.front().componentsCount(); ++component)
    {
        std::for_each(dataVectors.begin(), dataVectors.end(), [&](LinearImageType imageVector)
        {
            pcaMatrix.row(imageIndex) = imageVector;
            ++imageIndex;
        });
        std::vector<EigenLinearImageType> eigenImages = pca(pcaMatrix, nbComponents);
        for (int i = 0; i < eigenImages.size(); ++i)
        {
            this->eigenVecImages.at(i).setComponent(component, eigenImages.at(i));
        }
    }
}
