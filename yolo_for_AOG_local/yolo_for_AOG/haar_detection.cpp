#include "pch.h"
#include "haar_detection.h"

/** Global variables */
//String face_cascade_name = "haarcascade_frontalface_alt.xml";
//String eyes_cascade_name = "haarcascade_eye_tree_eyeglasses.xml";

string window_name = "Capture - Face detection";
RNG rng(12345);

void predictImageHaar(const string& filename, bool face_instead_of_eye) {
	// PreDefined trained XML classifiers with facial features 
	CascadeClassifier cascade, nestedCascade;
	double scale = 1;

	//-- 1. Load the cascades
	if (face_instead_of_eye) {
		if (!cascade.load("data\\haarcascade_frontalface_alt.xml")) {
			printf("--(!)Error loading\n"); return;
		}
		if (!nestedCascade.load("data\\haarcascade_eye.xml")) {
				printf("--(!)Error loading\n"); return;
		}
	}

	// Open the image file
	ifstream ifile(filename);
	if (!ifile) throw("error");
	VideoCapture cap;
	cap.open(filename);
	string outputFile = filename.substr(0, filename.size() - 4) + "_haar_face_out_cpp.jpg";

	Mat frame;
	cap >> frame;

	//face detection part
	vector<Rect> faces, faces2;
	Mat gray, smallImg;
	cvtColor(frame, gray, COLOR_BGR2GRAY);
	double fx = 1 / scale;

	// Resize the Grayscale Image  
	resize(gray, smallImg, Size(), fx, fx, INTER_LINEAR);
	equalizeHist(smallImg, smallImg);

	// Detect faces of different sizes using cascade classifier  
	cascade.detectMultiScale(smallImg, faces, 1.1,
		2, 0 | CASCADE_SCALE_IMAGE, Size(30, 30));

	// Draw circles around the faces 
	for (size_t i = 0; i < faces.size(); i++)
	{
		Rect r = faces[i];
		Mat smallImgROI;
		vector<Rect> nestedObjects;
		Point center;
		Scalar color = Scalar(255, 0, 0); // Color for Drawing tool 
		int radius;

		double aspect_ratio = (double)r.width / r.height;
		if (0.75 < aspect_ratio && aspect_ratio < 1.3)
		{
			center.x = cvRound((r.x + r.width*0.5)*scale);
			center.y = cvRound((r.y + r.height*0.5)*scale);
			radius = cvRound((r.width + r.height)*0.25*scale);
			circle(frame, center, radius, color, 3, 8, 0);
		}
		else
			rectangle(frame, Point(cvRound(r.x*scale), cvRound(r.y*scale)),
				Point(cvRound((r.x + r.width - 1)*scale),
					cvRound((r.y + r.height - 1)*scale)), color, 3, 8, 0);
		if (nestedCascade.empty())
			continue;
		smallImgROI = smallImg(r);

		//
		std::vector<Rect> eyes;

		//-- In each face, detect eyes
		nestedCascade.detectMultiScale(gray, eyes, 1.1, 1, 0 | CASCADE_SCALE_IMAGE, Size(5, 5));

		imshow("small part", smallImgROI);
		waitKey(0);

		for (size_t j = 0; j < eyes.size(); j++)
		{
			Point center(faces[i].x + eyes[j].x + eyes[j].width*0.5, faces[i].y + eyes[j].y + eyes[j].height*0.5);
			int radius = cvRound((eyes[j].width + eyes[j].height)*0.25);
			circle(frame, center, radius, Scalar(255, 0, 0), 4, 8, 0);
		}

		//// Detection of eyes int the input image 
		//nestedCascade.detectMultiScale(smallImgROI, nestedObjects, 1.1, 2,
		//	0 | CASCADE_SCALE_IMAGE, Size(30, 30));

		//imshow("small part", smallImgROI);
		//waitKey(0);
		//
		//cout << "nestedObjects len: " << nestedObjects.size() << endl;
		//// Draw circles around eyes 
		//for (size_t j = 0; j < nestedObjects.size(); j++)
		//{
		//	Rect nr = nestedObjects[j];
		//	center.x = cvRound((r.x + nr.x + nr.width*0.5)*scale);
		//	center.y = cvRound((r.y + nr.y + nr.height*0.5)*scale);
		//	radius = cvRound((nr.width + nr.height)*0.25*scale);
		//	circle(frame, center, radius, color, 3, 8, 0);
		//}
	}

	// Write the frame with the detection boxes
	Mat detectedFrame;
	frame.convertTo(detectedFrame, CV_8U);
	imwrite(outputFile, detectedFrame);
}