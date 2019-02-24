#ifndef HAAR_DETECTION_H
#define HAAR_DETECTION_H

#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>

#include "yolo_detection.h"

using namespace std;
using namespace cv;

/** Function Headers */
//void detectAndDisplay(Mat frame);
//void predictImageHaar(const string& filename, bool eye_in_face = false, string writefile = "");

void predictImageHaarCascade(const string& filename, string class_name, string writefile = "", string channel = "alpha");

#endif // !HAAR_DETECTION_H