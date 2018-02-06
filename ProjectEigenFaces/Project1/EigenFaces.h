#ifndef __EigenFaces_h_
#define __EigenFaces_h_
#include "CImg.h"
#include <string>
#include <vector>
#include "Eigen/Dense"
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
	typedef std::vector<Eigen::VectorXd>        EigenLinearImagePackType;

    EigenFaces() = delete;
    explicit EigenFaces(std::string, std::string);
    ~EigenFaces();

    void apply(const size_t);
    ImageType reconstruct(std::string);
protected:
    ImageType       image0;
    LinImgPackType  eigenVecImages;

    ImageType realign(const ImageType&, const ImageType&);
    void pca(LinImgPackType, size_t);
	EigenLinearImagePackType pca(Eigen::MatrixXd&, size_t&);
};

#endif //__EigenFaces_h_

#if not EigenFaces_Manual_Instanciation
#define EigenFaces_Manual_Instanciation 1
#include "EigenFaces.cpp"
#endif