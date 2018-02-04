#include "ImageParser.h"
#include <fstream>
#include <algorithm>
#include <iterator>

template<typename T>
ImageParser<T>::ImageParser(const std::string& filePath, const std::string& folderPath) : folderPath{ folderPath }
{
	std::ifstream in{ folderPath + "/" + filePath };
	if (in.is_open())
	{
		imagePaths(std::istream_iterator<std::string>(in), std::istream_iterator<std::string>());
		currentFilePath = imagePaths.begin();
	}
}

template<typename T>
inline std::vector<std::string>::iterator ImageParser<T>::begin()
{
	return imagePaths.begin();
}

template<typename T>
inline std::vector<std::string>::const_iterator ImageParser<T>::begin() const
{
	return imagePaths.cbegin();
}

template<typename T>
inline std::vector<std::string>::iterator ImageParser<T>::end()
{
	return imagePaths.end();
}

template<typename T>
inline std::vector<std::string>::const_iterator ImageParser<T>::end() const
{
	return imagePaths.cend();
}

template<typename T>
cimg_library::CImg<T> ImageParser<T>::load(const std::string& imageName, const std::string& imagePath) const
{
	return cimg_library::CImg<T>(imagePath + "/" + imageName);
}

template<typename T>
cimg_library::CImg<T> ImageParser<T>::next()
{
	if (currentImagePath != imagePaths.end())
	{
		cimg_library::CImg<T> image(folderPath+ "/" + (*currentImagePath));
		++it;
		return image;
	}
	return cimg_library::CImg<T>;
}