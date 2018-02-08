#ifndef __EigenFaces_h_
#define __EigenFaces_h_
#include "CImg.h"
#include <string>
#include <vector>
#include "Eigen/Dense"
#include "ImageVector.h"
#include "ImageParser.h"
#include "FGExtractor.h"
#include <string>

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

    void apply(const size_t, int);
    ImageType reconstruct(std::string);
protected:
    ImageType       image0;
    LinImgPackType  eigenVecImages;

    ImageType realign(const ImageType&, const ImageType&);
    void pca(LinImgPackType&, size_t, int);
	EigenLinearImagePackType pca(Eigen::MatrixXd&, size_t&);
};

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
void EigenFaces<T>::apply(const size_t nbComponents, int imgSize)
{
	LinImgPackType vectorizedImages;

	this->image0 = parser.next();
	vectorizedImages.emplace_back(this->image0.get_resize(imgSize, imgSize, 1, this->image0.spectrum()));

	std::for_each(parser.begin() + 1, parser.end(), [&](std::string a)
	{
		vectorizedImages.emplace_back(this->realign(this->image0, parser.next()).resize(imgSize, imgSize, 1, this->image0.spectrum()));
	});

	this->pca(vectorizedImages, nbComponents, imgSize);
}

template<class T>
typename EigenFaces<T>::ImageType EigenFaces<T>::reconstruct(std::string fileName)
{
	return ImageType();
}

template <class T>
typename EigenFaces<T>::ImageType EigenFaces<T>::realign(const ImageType& model, const ImageType& imageToAlign)
{
	std::vector<std::pair<int, int>> modelInterests, imageInterests;
	FGExtractor<unsigned char> extractor;
	modelInterests = extractor.getForegroundPixelsPositions(model.get_RGBtoYCbCr().get_channel(0).equalize(255, 0, 255), 5);
	imageInterests = extractor.getForegroundPixelsPositions(imageToAlign.get_RGBtoYCbCr().get_channel(0).equalize(255, 0, 255), 5);

	int N1 = modelInterests.size();
	int N2 = imageInterests.size();

	Eigen::MatrixXd F1(2, N1);
	Eigen::MatrixXd F2(2, N2);

	for (int i = 0; i < N1; ++i)
	{
		F1(0, i) = modelInterests[i].first;
		F1(1, i) = modelInterests[i].second;
	}
	for (int i = 0; i < N2; ++i)
	{
		F2(0, i) = imageInterests[i].first;
		F2(1, i) = imageInterests[i].second;
	}

	size_t modelPcaMaxIndex = 0;
	size_t imagePcaMaxIndex = 0;
	std::vector< Eigen::VectorXd > pcaModel = this->pca(F1, modelPcaMaxIndex);
	std::vector< Eigen::VectorXd > pcaImage = this->pca(F2, imagePcaMaxIndex);

	Eigen::VectorXd modelMainAxis = pcaModel[modelPcaMaxIndex];
	Eigen::VectorXd imageMainAxis = pcaImage[imagePcaMaxIndex];

	double dotprod(imageMainAxis.dot(modelMainAxis));
	double angle(acos(dotprod)*180.0 / 3.14159);


	ImageType rotatedImage = imageToAlign;
	rotatedImage.rotate(-1.0*angle);

	ImageType outImage(imageToAlign.width(), imageToAlign.height(), imageToAlign.depth(), imageToAlign.spectrum());
	if (rotatedImage.width() > imageToAlign.width() &&
		rotatedImage.height() > imageToAlign.height())
	{
		size_t startX = 0;
		size_t startY = 0;

		startX = std::max((rotatedImage.width() / 2) - (outImage.width() / 2), 0);
		startY = std::max((rotatedImage.height() / 2) - (outImage.height() / 2), 0);

		cimg_forXY(outImage, x, y)
		{
			for (size_t chanel = 0; chanel < outImage.spectrum(); ++chanel)
			{
				outImage(x, y, chanel) = rotatedImage(x + startX, y + startY, chanel);
			}
		}
	}
	else
	{
		outImage = rotatedImage;
	}

	return std::move(outImage);
}

template<class T>
void EigenFaces<T>::pca(LinImgPackType& dataVectors, size_t nbComponents, int imgSize)
{
	this->eigenVecImages.resize(nbComponents);
	std::for_each(this->eigenVecImages.begin(), this->eigenVecImages.end(), [&](LinearImageType& vectorImage)
	{
		vectorImage.resize(dataVectors.front().componentsCount());
		vectorImage.setPixelSize(nbComponents);
	});
	size_t sz = dataVectors.front().pixelCount();
	Eigen::MatrixXd pcaMatrix(dataVectors.front().pixelCount(), dataVectors.size());
	int imageIndex = 0;
	for (int component = 0; component < dataVectors.front().componentsCount(); ++component)
	{
		std::for_each(dataVectors.begin(), dataVectors.end(), [&](LinearImageType imageVector)
		{
			pcaMatrix.col(imageIndex) = imageVector.getComponent(component).cast<double>();
			++imageIndex;
		});
		std::vector<EigenLinearImageType> eigenImages = pca(pcaMatrix, sz);
		for (int i = 0; i < nbComponents; ++i)
		{
			this->eigenVecImages.at(i).setComponentFree(component, eigenImages.at(i));
		}
		imageIndex = 0;
	}
	for (int i = 0; i < nbComponents; ++i)
	{
		this->eigenVecImages.at(i).save("eigenV" + std::to_string(i) + ".ppm", imgSize, imgSize);
	}
}

template<class T>
typename EigenFaces<T>::EigenLinearImagePackType EigenFaces<T>::pca(Eigen::MatrixXd& data, size_t& maxIndex)
{
	int N = data.cols();
	if (N == 0)
	{
		return std::vector< EigenLinearImageType >();
	}

	EigenLinearImageType Xbarre = data.rowwise().mean();
	EigenLinearImageType O(N);
	O.setOnes();

	Eigen::MatrixXd Q = (1.0 / (N - 1.0)) * (data - Xbarre * O.transpose()) * (data - Xbarre * O.transpose()).transpose();

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ES(Q);

	ES.eigenvalues().real().maxCoeff(&maxIndex);

	std::vector< EigenLinearImageType > eigenVectors;

	for (int i = 0; i < data.rows(); ++i)
	{
		EigenLinearImageType vp = ES.eigenvectors().col(i).real();
		eigenVectors.emplace_back(vp);
	}
	return eigenVectors;
}

#endif //__EigenFaces_h_