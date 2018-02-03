#ifndef IMAGE_PARSER_H_ 
#define IMAGE_PARSER_H_

#include "CImg.h"
#include <string>
#include <vector>
 //T1 for the returned cimg format, T2 and T3 assume char*, char[] or string
template<typename T>
class ImageParser
{
	std::vector<std::string> filePaths; //
	std::string folderPath;
	std::vector<std::string>::const_iterator currentFilePath;
public:
	ImageParser() = delete;
	ImageParser(const std::string&, const std::string&);
	ImageParser(const ImageParser&) = delete;
	ImageParser(ImageParser&&) = delete;
	~ImageParser() = default;

	std::vector<std::string>::iterator begin() const;
	std::vector<std::string>::const_iterator begin() const;
	std::vector<std::string>::iterator end() const;
	std::vector<std::string>::const_iterator end() const;
	cimg_library::CImg<T> load(const std::string& imageName, const std::string& imagePath) const;
	cimg_library::CImg<T> next();
};
#endif // !IMAGE_PARSER_H_ 
