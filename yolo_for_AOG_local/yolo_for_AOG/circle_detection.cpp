#include "pch.h"
#include "circle_detection.h"

void predictCircle(std::string filename, std::string writefile) {
	// Open the image file
	ifstream ifile(filename);
	if (!ifile) throw("error");
	VideoCapture cap;
	cap.open(filename);
	string outputFile = filename.substr(0, filename.size() - 4) + "_circle_out_cpp.jpg";

	Mat frame;
	cap >> frame;
	Mat gray;
	cvtColor(frame, gray, COLOR_BGR2GRAY);

	GaussianBlur(gray, gray, Size(9, 9), 6, 6);
	vector<Vec3f> circles;
	//HoughCircles(gray, circles, HOUGH_GRADIENT,
	//	2, gray.rows / 4, 200, 100);
	HoughCircles(gray, circles, HOUGH_GRADIENT,
		2, 20, 450, 60, 0, 0);
	for (size_t i = 0; i < circles.size(); i++)
	{
		Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
		int radius = cvRound(circles[i][2]);
		// draw the circle center
		circle(frame, center, 3, Scalar(0, 255, 0), -1, 8, 0);
		// draw the circle outline
		circle(frame, center, radius, Scalar(0, 0, 255), 3, 8, 0);
	}

	namedWindow("circles", 1);
	imshow("circles", frame);
	waitKey(0);
}