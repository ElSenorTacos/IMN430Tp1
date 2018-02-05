#ifndef __EigenFaces_h_
#define __EigenFaces_h_
#include "CImg.h"
#include <string>
#include <vector>
#include "Eigen/Cholesky"
#include "ImageVector.h"
#include "ImageParser.h"

template<class T>
class EigenFaces
{
    ImageParser<T> parser;
public:

    typedef cimg_library::CImg<T>               ImageType;
    typedef ImageVector<T>                      LinearImageType;
    typedef std::vector<LinearImageType>        LinImgPackType;
    typedef Eigen::VectorXd                     EigenLinearImageType;

    EigenFaces() = delete;
    explicit EigenFaces(std::string, std::string);
    ~EigenFaces();

    EigenFaces<T>* apply(size_t);
    ImageType reconstruct(std::string);
protected:
    ImageType       image0;
    LinImgPackType  eigenVecImages;

    ImageType realign(const ImageType&, const ImageType&) const;
    void pca(LinImgPackType, size_t);
	std::vector<EigenLinearImageType> pca(const Eigen::MatrixXd&, size_t&);
};

#endif //__EigenFaces_h_