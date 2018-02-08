
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
	cout << "What size bra : ";
	cin >> sz;
	cout << "Nb of v_prop bra : ";
	cin >> nbVects;


    EigenFaces<unsigned char> eigenFaces(argv[1], argv[2]);
    eigenFaces.apply(nbVects, sz);

	return 0;
}
