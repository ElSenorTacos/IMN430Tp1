
#include "Eigen/Dense"
#include "EigenFaces.h"
#include <iostream>

	using namespace Eigen;
	using namespace cimg_library;
	using namespace std;

int main(int argc, const char * argv[]) {
	if (argc < 3) {
		std::cout << "Set argument 1 to the filePath containing the images" << std::endl;
		std::cout << "Set argument 2 to either : - the fileNames containing the images names of the database" << std::endl;
		std::cout << "                           - the fileNames containing the eigenvectors of the database" << std::endl;
		return 1;
	}

	size_t sz, nbVects;



	EigenFaces<unsigned char> eigenFaces(argv[1], argv[2]);
	if (!eigenFaces.existingDB()) {
		cout << "What size : ";
		cin >> sz;
		eigenFaces.apply(sz);
	}
	string recons;
	while (true) {
		cout << "Image to rebuild : ";
		cin >> recons;
		cout << "Nb of v_prop : ";
		cin >> nbVects;

		Eigen::VectorXd reconstructionCoeffs;
		EigenFaces<unsigned char>::ImageType output = eigenFaces.reconstruct(recons, nbVects, reconstructionCoeffs);

		string ouputn = "recons_" + std::to_string(nbVects) + "_" + recons + ".pgm";
		output.save(ouputn.c_str());
		cout << "Another ? ";
		cin >> recons;
		if (recons == "n" || recons == "N") {
			break;
		}
	}

	return 0;
}
