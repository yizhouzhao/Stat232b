#ifndef CIRCLE_DETECTION_H
#define CIRCLE_DETECTION_H

#include <fstream>
//#include <sstream>
#include <iostream>

//#include "opencv2/objdetect/objdetect.hpp"
#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/Imgproc.hpp>
#include <opencv2/dnn.hpp>

using namespace cv;
using namespace dnn;
using namespace std;

void predictCircle(std::string filename, std::string writefile = "");



#endif // !CIRCLE_DETECTION_H
