#include <iostream>
#include <vector>
#include <cmath>
//todo remove includes when parse works
#include <fstream>
#include <algorithm>
#include <streambuf>
#include <iterator>
//
#include "Eigen/Dense"

#define cimg_display 0
#include "CImg.h"
#include "ImageParser.h"

using namespace Eigen;
using namespace cimg_library;
using namespace std;

int main(int argc, const char * argv[]) {
	ImageParser<unsigned char> image("resized", "name.txt");
	//CImg<unsigned char>test("1.png");
	ImageParser<unsigned char> testImageParser("test", "testname.txt");
	CImg<unsigned char> testCImg = testImageParser.next();
	testCImg.save("output.ppm");
	testImageParser.load("test", "oval1.ppm").save("teleportedOval");
	testImageParser.setBegin();
	unsigned int index = 0;
	while (true)
	{
		CImg<unsigned char> current = testImageParser.next(); 
		if (current != CImg<unsigned char>())
		{
			current.save(std::string(std::to_string(index) + ".ppm").c_str());
			++index;
		}
		else
		{
			break;
		}
	}

//	CImg<unsigned char> image1("oval1.ppm");

	//CImg<unsigned char> image1b(image1.get_RGBtoYCbCr().get_channel(0));
	/*	vector<int> xC, yC;

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

	image1.save("output.ppm");*/

	return 0;
}