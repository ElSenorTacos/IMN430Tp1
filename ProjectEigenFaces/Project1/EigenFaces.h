#ifndef __EigenFaces_h_
#define __EigenFaces_h_
#include "CImg.h"
#include <string>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Householder"
#include "ImageVector.h"
#include "ImageParser.h"
#include "FGExtractor.h"
#include <string>

struct eigen {
	double value;
	Eigen::VectorXd *vect;
	eigen() : value{ std::numeric_limits<double>::min() }, vect{ nullptr } {}
	eigen(eigen& oth) : value{ oth.value }, vect{ oth.vect } {}
	eigen(eigen&& oth) : value{ oth.value }, vect{ oth.vect } {}
	explicit eigen(double valp, Eigen::VectorXd& vp) : value{ valp }, vect{ &vp } {}
	bool operator<(eigen const &oth) const {
		return value >= oth.value;
	}
	eigen operator=(eigen& oth) {
		return oth;
	}
};
struct myclass {
	bool operator() (eigen& i, eigen& j) { return (i.value>=j.value); }
} bigger;

template<class T>
class EigenFaces
{

    ImageParser<T> parser;
public:

    typedef cimg_library::CImg<T>               ImageType;
    typedef ImageVector<T>                      LinearImageType;
    typedef std::vector<LinearImageType>        LinImgPackType;
	typedef cimg_library::CImg<double>			EigenLinImgType;
	typedef std::vector<ImageVector<double>>    EigenLinImgPackType;
    typedef Eigen::VectorXd                     EigenLinearImageType;
	typedef std::vector<Eigen::VectorXd>        EigenLinearImagePackType;

    EigenFaces() = delete;
    explicit EigenFaces(std::string, std::string);
    ~EigenFaces();

    void apply(const size_t, int);
    ImageType reconstruct(std::string, int);
protected:
    ImageType			 image0;
	EigenLinImgPackType  eigenVecImages;
	EigenLinImgType		 mean;

    ImageType realign(const ImageType&, const ImageType&);
    void pca(LinImgPackType&, size_t, int);
	EigenLinearImagePackType pca(Eigen::MatrixXd&, size_t&, bool = false);
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
	vectorizedImages.emplace_back(this->image0.get_resize(imgSize, imgSize, 1, this->image0.spectrum()).get_RGBtoYCbCr().get_channel(0));

	std::for_each(parser.begin() + 1, parser.end(), [&](std::string a)
	{
		vectorizedImages.emplace_back(/*this->realign(this->image0, */parser.next()/*)*/.resize(imgSize, imgSize, 1, this->image0.spectrum()).get_RGBtoYCbCr().get_channel(0));
	});

	this->pca(vectorizedImages, nbComponents, imgSize);
}

template<class T>
typename EigenFaces<T>::ImageType EigenFaces<T>::reconstruct(std::string fileName, int imgSize)
{
	EigenLinImgType output(imgSize, imgSize, this->image0.depth(), this->image0.spectrum()); output.fill(0);
	EigenLinImgType ref(/*this->realign(image0, */parser.load(fileName)/*)*/.get_resize(imgSize, imgSize, this->image0.depth(), this->image0.spectrum()));
	ImageVector<double> dummy(ref - mean);
	for (int ch = 0; ch < this->eigenVecImages.front().componentsCount(); ++ch) {
		Eigen::VectorXd c(this->eigenVecImages.size());
		Eigen::MatrixXd A(this->eigenVecImages.size(), this->eigenVecImages.front().pixelCount());
		for (int i = 0; i < this->eigenVecImages.size(); ++i) {
			ImageVector<double> vecp((this->eigenVecImages[i].getImage(imgSize, imgSize)));
			A.row(i) = vecp.getComponent(ch);
		}
		c = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(A.transpose()).solve(dummy.getComponent(ch));
		for (int i = 1; i < c.rows(); ++i)
		{
			output.channel(ch) += c[i] * this->eigenVecImages[i].getImage(imgSize, imgSize).channel(ch);
		}
	}
	output += mean;
	return output.abs();
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
	std::for_each(this->eigenVecImages.begin(), this->eigenVecImages.end(), [&](ImageVector<double>& vectorImage)
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
		this->mean = ImageVector<double>::vectToImage(pcaMatrix.rowwise().mean(), imgSize, imgSize);
		std::vector<EigenLinearImageType> eigenImages = pca(pcaMatrix, sz, true);
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
typename EigenFaces<T>::EigenLinearImagePackType EigenFaces<T>::pca(Eigen::MatrixXd& data, size_t& maxIndex, bool sorted)
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

	std::vector< EigenLinearImageType > eigenVectors;
	VectorXd eigenVals{ ES.eigenvalues() };

	for (int i = 0; i < data.rows(); ++i)
	{
		eigenVectors.emplace_back(ES.eigenvectors().col(i).real());
	}

	if (sorted) {
		std::vector<eigen> vpsorted;
		for (int i = 0; i < eigenVectors.size(); ++i) {
			vpsorted.emplace_back(eigenVals[i], eigenVectors[i]);
		}
		std::sort(vpsorted.begin(), vpsorted.end(), bigger);
		std::vector<EigenLinearImageType> eig2;
		std::for_each(vpsorted.rbegin(), vpsorted.rend(), [&](eigen e) {
			eig2.push_back(*(e.vect));
		});
		eigenVectors = eig2;
	}
	else {
		ES.eigenvalues().maxCoeff(&maxIndex);
	}

	return eigenVectors;
}

#endif //__EigenFaces_h_