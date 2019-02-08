#ifndef HAAR_DETECTION_H
#define HAAR_DETECTION_H
#endif // !HAAR_DETECTION_H

#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>

using namespace std;
using namespace cv;

/** Function Headers */
void detectAndDisplay(Mat frame);

void predictImageHaar(const string& filename, bool face_instead_of_eye = true);
