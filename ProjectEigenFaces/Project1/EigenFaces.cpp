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
typename EigenFaces<T>::ImageType EigenFaces<T>::realign(const ImageType& model, const ImageType& imageToAlign) const
{
	std::vector<int> modelX, modelY;
	std::vector<int> imageX, imageY;

	cimg_forXY(model, x, y)
	{
		bool awesomePixel = false;
		for (int chanel = 0; chanel < model.spectrum(); ++chanel)
		{
			if (model.get_channel(chanel)(x, y) > 50)
			{
				awesomePixel = true;
			}
		}

		if (awesomePixel)
		{
			modelX.push_back(x);
			modelY.push_back(y);
		}
	}

	cimg_forXY(imageToAlign, x, y)
	{
		bool awesomePixel = false;
		for (int chanel = 0; chanel < imageToAlign.spectrum(); ++chanel)
		{
			if (imageToAlign.get_channel(chanel)(x, y) > 50)
			{
				awesomePixel = true;
			}
		}

		if (awesomePixel)
		{
			imageX.push_back(x);
			imageY.push_back(y);
		}
	}

	int N1 = modelX.size();
	int N2 = imageX.size();

	Eigen::MatrixXd F1(2, N1);
	Eigen::MatrixXd F2(2, N2);

	for (int i = 0; i < N1; ++i)
	{
		F1(0, i) = modelX[i];
		F1(1, i) = modelY[i];
	}
	for (int i = 0; i < N2; ++i)
	{
		F2(0, i) = imageX[i];
		F2(1, i) = imageY[i];
	}

	double modelPcaMaxIndex = 0;
	double imagePcaMaxIndex = 0;
	std::vector< Eigen::VectorXd > pcaModel = pca(F1, modelPcaMaxIndex);
	std::vector< Eigen::VectorXd > pcaImage = pca(F2, modelPcaMaxIndex);

	Eigen::VectorXd modelMainAxis = pcaModel[modelPcaMaxIndex];
	Eigen::VectorXd imageMainAxis = pcaImage[imagePcaMaxIndex];

	double dotprod(imageMainAxis.dot(modelMainAxis));
	double angle(acos(dotprod)*180.0 / 3.14159);

	ImageType outImage = imageToAlign;
	outImage.rotate(-1.0*angle);

	return std::move(outImage);
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

template<class T>
std::vector<typename EigenFaces<T>::EigenLinearImageType> pca(const Eigen::MatrixXd& data, size_t& maxIndex)
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

	Eigen::EigenSolver<Eigen::MatrixXd> ES(Q);

	ES.eigenvalues().real().maxCoeff(&maxIndex);

	std::vector< EigenLinearImageType > eigenVectors;

	for (int i = 0; i < data.rows(); ++i)
	{
		EigenLinearImageType vp = ES.eigenvectors().col(i).real();
		eigenVectors.emplace_back(vp);
	}
	return eigenVectors;
}
