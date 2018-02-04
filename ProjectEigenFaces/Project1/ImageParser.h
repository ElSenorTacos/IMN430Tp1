#ifndef IMAGE_PARSER_H_ 
#define IMAGE_PARSER_H_

#include "CImg.h"
#include <string>
#include <vector>

template<typename T>
class ImageParser
{
	std::string folderPath;
	std::vector<std::string>::const_iterator currentImagePath;
	std::vector<std::string> imagePaths;
public:
	ImageParser() = delete;
	ImageParser(const std::string&, const std::string&);
	ImageParser(const ImageParser&) = delete;
	ImageParser(ImageParser&&) = delete;
	~ImageParser() = default;

	std::vector<std::string>::iterator begin();
	std::vector<std::string>::const_iterator begin() const;
	std::vector<std::string>::iterator end();
	std::vector<std::string>::const_iterator end() const;
	cimg_library::CImg<T> load(const std::string& imageName, const std::string& imagePath) const;
	cimg_library::CImg<T> next();
};
#endif // !IMAGE_PARSER_H_ 
