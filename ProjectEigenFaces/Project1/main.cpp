
#include "Eigen/Dense"
#include "EigenFaces.h"
#include <iostream>

	using namespace Eigen;
	using namespace cimg_library;
	using namespace std;

int main(int argc, const char * argv[]) {
    if(argc < 3) return 1;

	/*for (int i = 1; i < 47; i++) {
		std::string fileX = std::to_string(i) + ".ppm";
		std::string fileXOutput = "outpout" + fileX;
		CImg<unsigned char> image1(fileXOutput.c_str());

		image1.resize(8, 8);
		image1.save(fileXOutput.c_str());
	}*/
	size_t sz, nbVects;
	cout << "What size : ";
	cin >> sz;
	cout << "Nb of v_prop : ";
	cin >> nbVects;


    EigenFaces<unsigned char> eigenFaces(argv[1], argv[2]);
    eigenFaces.apply(nbVects, sz);
	string recons, ext;
	while (true) {
		cout << "Image to rebuild (without extension) : ";
		cin >> recons;
		cout << "Extension : ";
		cin >> ext;
		EigenFaces<unsigned char>::ImageType output = eigenFaces.reconstruct(recons + ext, sz);
		string ouputn = recons + "_recons" + ext;
		output.save(ouputn.c_str());
		cout << "Another ? ";
		cin >> recons;
		if (recons == "n" || recons == "N") {
			break;
		}
	}

	return 0;
}
