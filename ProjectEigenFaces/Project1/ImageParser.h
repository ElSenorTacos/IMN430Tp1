#ifndef IMAGE_PARSER_H_ 
#define IMAGE_PARSER_H_

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>
#include "CImg.h"

template<typename T>
class ImageParser //Not robust: assumes that image name extensions are CImg compatible
{
	std::string folderPath;
	std::vector<std::string>::const_iterator currentImagePath;
	std::vector<std::string> imagePaths;
public:
	ImageParser() = delete;
	ImageParser(const std::string& folderPath, const std::string& filePath) : folderPath{ folderPath }
	{
		std::ifstream in{ folderPath + "/" + filePath };
		if (in.is_open())
		{
			std::copy(std::istream_iterator<std::string>(in), std::istream_iterator<std::string>(), std::back_inserter(imagePaths));
			currentImagePath = imagePaths.begin();
		}
	}
	ImageParser(const ImageParser&) = delete;
	ImageParser(ImageParser&&) = delete;
	~ImageParser() = default;

	inline std::vector<std::string>::iterator begin() { return imagePaths.begin(); }
	inline std::vector<std::string>::const_iterator begin() const { return imagePaths.cbegin(); }
	inline std::vector<std::string>::iterator end() { return imagePaths.end(); }
	inline std::vector<std::string>::const_iterator end() const { return imagePaths.cend(); }
	inline void setBegin() { currentImagePath = imagePaths.begin(); }
	inline const size_t size() const { return imagePaths.size(); }
	cimg_library::CImg<T> load(const std::string& genericFolderPath, const std::string& imageName) const //load image wrapper
	{
		return cimg_library::CImg<T>(std::string(genericFolderPath + "/" + imageName).c_str());
	}

	cimg_library::CImg<T> load(const std::string& imageName) const //load image wrapper using internal folder path
	{
		return load(folderPath, imageName);
	}

	cimg_library::CImg<T> next() //return current image and place internal iterator on next image : return empty CImg if end
	{
		if (currentImagePath != imagePaths.end())
		{
			cimg_library::CImg<T> image(std::string(folderPath + "/" + *currentImagePath).c_str());
			++currentImagePath;
			return image;
		}
		return cimg_library::CImg<T>();
	}
};
#endif // !IMAGE_PARSER_H_ 