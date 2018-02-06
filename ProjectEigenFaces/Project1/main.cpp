#include <iostream>
#include <vector>
#include <cmath>

#include "Eigen/Dense"

#define cimg_display 0
#include "CImg.h"

using namespace Eigen;
using namespace cimg_library;
using namespace std;

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

int main(int argc, const char * argv[]) {
    /*
	CImg<unsigned char> image1("oval1.ppm");

	CImg<unsigned char> image1b(image1.get_RGBtoYCbCr().get_channel(0));

	vector<int> xC, yC;

	cimg_forXY(image1b, x, y)
	{
		if (image1b(x, y) > 50)
		{
			xC.push_back(x);
			yC.push_back(y);
		}
	}

	int N = xC.size();

	MatrixXd F(2, N);

	for (int i = 0; i < N; ++i)
	{
		F(0, i) = xC[i];
		F(1, i) = yC[i];
	}

	Vector2d Xbarre = F.rowwise().mean();
	VectorXd O(N);
	O.setOnes();

	MatrixXd Q = (1.0 / (N - 1.0)) * (F - Xbarre*O.transpose()) * (F - Xbarre*O.transpose()).transpose();

	EigenSolver<MatrixXd> ES(Q);

	int maxIndex;
	ES.eigenvalues().real().maxCoeff(&maxIndex);

	Vector2d vp1 = ES.eigenvectors().col(maxIndex).real();

	Vector2d yaxis;
	yaxis << 0, -1;

	double dotprod(vp1.dot(yaxis));
	double angle(acos(dotprod)*180.0 / 3.14159);

	image1.rotate(-1.0*angle);

	image1.save("output.ppm");
     */

	CImg<unsigned char> image1("16.ppm");
	CImg<unsigned char> image2("6.ppm");

	vector<int> imageX; vector<int> imageY;

	CImg<unsigned char> image2Aligned = realign(image1, image2);

	image2Aligned.save("output3.ppm");



	return 0;
}