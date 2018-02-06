
#include "Eigen/Dense"
#include "EigenFaces.h"

using namespace Eigen;

int main(int argc, const char * argv[]) {
    if(argc < 3) return 1;

    EigenFaces<unsigned char> eigenFaces(argv[0], argv[1]);
    eigenFaces.apply(5);

	return 0;
}