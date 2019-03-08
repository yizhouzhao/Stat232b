#ifndef VIDEO_TO_IMAGE_H
#define VIDEO_TO_IMAGE_H

#include <iostream>
#include <fstream>
#include <direct.h>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include "yolo_detection.h"
#include "haar_detection.h"

using namespace cv;

void Video2Image(std::string input_file, std::string output_folder, int skip = 10);

void VideoAddBox(std::string input_file, std::string box_folder, int skip = 10);


#endif // !VIDEO_TO_IMAGE_H
