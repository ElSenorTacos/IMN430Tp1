#include <iostream>
#include <vector>
#include <cmath>

#include "Eigen/Dense"

#define cimg_display 0
#include "CImg.h"

using namespace Eigen;
using namespace cimg_library;
using namespace std;


vector< VectorXd > pca(MatrixXd data, double& maxIndex)
{
	int N = data.cols();
	if (N == 0)
	{
		return vector< VectorXd >();
	}


	VectorXd Xbarre = data.rowwise().mean();
	VectorXd O(N);
	O.setOnes();

	MatrixXd Q = (1.0 / (N - 1.0)) * (data - Xbarre * O.transpose()) * (data - Xbarre * O.transpose()).transpose();

	EigenSolver<MatrixXd> ES(Q);

	ES.eigenvalues().real().maxCoeff(&maxIndex);

	vector< VectorXd > eigenVectors;
																						
	for (int i = 0; i < data.rows(); ++i)
	{
		VectorXd vp = ES.eigenvectors().col(i).real();

		eigenVectors.emplace_back(vp);
	}
	return eigenVectors;
}



/*
void filtreGaussien5X5(double filtre[][5], double sig)
{
	double r, s = 2.0 * sig * sig;
	double somme = 0.0;

	for (int x = -2; x <= 2; ++x)
	{
		for (int y = -2; y <= 2; ++y)
		{
			r = sqrt(x*x + y * y);
			filtre[x + 2][y + 2] = (exp(-(r*r) / s)) / (3.14159265359 * s);
			somme += filtre[x + 2][y + 2];
		}
	}

	// normalising the Kernel
	for (int i = 0; i < 5; ++i)
	{
		for (int j = 0; j < 5; ++j)
		{
			filtre[i][j] /= somme;
		}
	}
}*/

double average(const CImg<unsigned char>& image, double& sigma, size_t X, size_t Y)
{
	double Xcenter;
	double Ycenter;
	if (X < 2)
	{
		Xcenter = 2;
	}
	else if(X > (image.width() - 3))
	{
		Xcenter = image.width() - 3;
	}
	else
	{
		Xcenter = X;
	}

	if (Y < 2)
	{
		Ycenter = 2;
	}
	else if (Y >(image.height() - 3))
	{
		Ycenter = image.height() - 3;
	}
	else
	{
		Ycenter = Y;
	}


	double somme = 0;
	for (int i = -2; i <= 2; ++i)
	{
		for (int j = -2; j <= 2; ++j)
		{
			somme += image(Xcenter + i, Ycenter + j);
		}
	}
	double mu = somme / 25.0f;

	sigma = 0;
	for (int i = -2; i <= 2; ++i)
	{
		for (int j = -2; j <= 2; ++j)
		{
			sigma += pow((image(Xcenter + i, Ycenter + j) - mu),2);
			
		}
	}

	sigma = sigma / 25.0f;
	//cout << "sig " << sigma << endl;
	return mu;
}

CImg<unsigned char> filterPoints(const CImg<unsigned char>& image, vector<int>& imageX, vector<int>& imageY)
{
	imageX.clear();
	imageY.clear();
	CImg<unsigned char> copy(image);
	
	CImg<unsigned char> imageGrey = copy.get_RGBtoYCbCr().get_channel(0).equalize(255,0,255);
	imageGrey.save("equal.ppm");
	
	double sig = 0;


	cimg_forXY(imageGrey, x, y)
	{

		double mu = average(imageGrey, sig, x, y);

		if (abs((imageGrey(x,y) - mu) / sig) > 0.5)
		{
			imageX.push_back(x);
			imageY.push_back(y);
		}
	}

	CImg<unsigned char> imageFiltered( image.width(), image.height(), image.depth(), image.spectrum());
	imageFiltered.fill(0);
	for (int i = 0; i < imageX.size(); ++i)
	{
		imageFiltered(imageX[i], imageY[i],0) = image(imageX[i], imageY[i],0);
		imageFiltered(imageX[i], imageY[i],1) = image(imageX[i], imageY[i],1);
		imageFiltered(imageX[i], imageY[i],2) = image(imageX[i], imageY[i],2);
	}	


	return imageFiltered;
}

CImg<unsigned char> realign(const CImg<unsigned char>& model, const CImg<unsigned char>& imageToAlign)
{
	vector<int> modelX, modelY;
	vector<int> imageX, imageY;

	CImg<unsigned char> modelFiltered = filterPoints(model, modelX, modelY);
	CImg<unsigned char> imageFiltered = filterPoints(imageToAlign, imageX, imageY);
	modelFiltered.save("BLUR1.ppm");
	imageFiltered.save("BLUR2.ppm");

	int N1 = modelX.size();
	int N2 = imageX.size();

	MatrixXd F1(2, N1);
	MatrixXd F2(2, N2);

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
	vector< VectorXd > pcaModel = pca(F1, modelPcaMaxIndex);
	vector< VectorXd > pcaImage = pca(F2, imagePcaMaxIndex);

	VectorXd modelMainAxis = pcaModel[modelPcaMaxIndex];
	VectorXd imageMainAxis = pcaImage[imagePcaMaxIndex];

	double dotprod(imageMainAxis.dot(modelMainAxis));
	double angle(acos(dotprod)*180.0 / 3.14159);

	CImg<unsigned char> outImage = imageToAlign;
	outImage.rotate(-1.0*angle);

	return outImage;
}

int main(int argc, const char * argv[]) {

	CImg<unsigned char> image1("16.ppm");
	CImg<unsigned char> image2("6.ppm");
	
	vector<int> imageX; vector<int> imageY;

	CImg<unsigned char> image2Aligned = realign(image1, image2);

	image2Aligned.save("output3.ppm");



	return 0;
}