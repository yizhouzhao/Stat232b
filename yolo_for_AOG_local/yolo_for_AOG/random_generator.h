#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <fstream>
#include <random>

using namespace cv;
using namespace std;

void generateRectanglesOnFrame(Mat frame, std::string filename="", int num = 20);

void writeRandomRectToFile(Mat frame, int num, std::string file_name);

#endif // !RANDOM_GENERATOR_H