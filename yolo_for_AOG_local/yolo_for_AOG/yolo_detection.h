#ifndef YOLO_DETECTION_H
#define YOLO_DETECTION_H

#include <fstream>
#include <sstream>
#include <iostream>

#include <opencv2/dnn.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

using namespace cv;
using namespace dnn;
using namespace std;

// Remove the bounding boxes with low confidence using non-maxima suppression
void postprocess(Mat& frame, const vector<Mat>& out, vector<string>& classes, string class_name, string writefile = "", string channel = "gamma");

// Draw the predicted bounding box
void drawPred(int classId, float conf, int left, int top, int right, int bottom, Mat& frame, vector<string>& classes);

// Get the names of the output layers
vector<String> getOutputsNames(const Net& net);

void predictImageYolov3(string filename, string class_name, string writefile = "", string channel = "gamma");

vector<string> getNamesOfClasses(const string& classesFile);

#endif //!YOLO_DETECTION_H