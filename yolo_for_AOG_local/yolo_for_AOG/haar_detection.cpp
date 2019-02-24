#include "pch.h"
#include "haar_detection.h"
#include <algorithm>

/** Global variables */
//String face_cascade_name = "haarcascade_frontalface_alt.xml";
//String eyes_cascade_name = "haarcascade_eye_tree_eyeglasses.xml";

string window_name = "Capture - Face detection";
RNG rng(12345);

//void predictImageHaar(const string& filename, bool eye_in_face, string writefile) {
//	// PreDefined trained XML classifiers with facial features 
//	CascadeClassifier cascade, nestedCascade;
//	double scale = 1;
//
//	//-- 1. Load the cascades
//	if (!cascade.load("data\\haarcascade_frontalface_alt.xml")) {
//		printf("--(!)Error loading\n"); return;
//	}
//	if (!nestedCascade.load("data\\haarcascade_eye.xml")) {
//			printf("--(!)Error loading\n"); return;
//	}
//
//
//	// Open the image file
//	ifstream ifile(filename);
//	if (!ifile) throw("error");
//	VideoCapture cap;
//	cap.open(filename);
//	string outputFile = filename.substr(0, filename.size() - 4) + "_haar_face_out_cpp.jpg";
//
//	Mat frame;
//	cap >> frame;
//
//	//face detection part
//	vector<Rect> faces, faces2;
//	Mat gray, smallImg;
//	cvtColor(frame, gray, COLOR_BGR2GRAY);
//	double fx = 1 / scale;
//
//	// Resize the Grayscale Image  
//	resize(gray, smallImg, Size(), fx, fx, INTER_LINEAR);
//	equalizeHist(smallImg, smallImg);
//
//	// Detect faces of different sizes using cascade classifier  
//	cascade.detectMultiScale(smallImg, faces, 1.1,
//		2, 0 | CASCADE_SCALE_IMAGE, Size(30, 30));
//
//	// Draw circles around the faces 
//	for (size_t i = 0; i < faces.size(); i++)
//	{
//		Rect r = faces[i];
//		
//		ofstream file(writefile, std::ios_base::app);
//		if (file.is_open()) {
//			file << "alpha face " << r.x << " " << r.y << " " << r.x + r.width << " " << r.y + r.height << "\n";
//			file.close();
//		}
//
//		Mat smallImgROI;
//		vector<Rect> nestedObjects;
//		Point center;
//		Scalar color = Scalar(255, 0, 0); // Color for Drawing tool 
//		int radius;
//
//		double aspect_ratio = (double)r.width / r.height;
//		if (false) //0.75 < aspect_ratio && aspect_ratio < 1.3)
//		{
//			center.x = cvRound((r.x + r.width*0.5)*scale);
//			center.y = cvRound((r.y + r.height*0.5)*scale);
//			radius = cvRound((r.width + r.height)*0.25*scale);
//			circle(frame, center, radius, color, 3, 8, 0);
//		}
//		else
//			rectangle(frame, Point(cvRound(r.x*scale), cvRound(r.y*scale)),
//				Point(cvRound((r.x + r.width - 1)*scale),
//					cvRound((r.y + r.height - 1)*scale)), color, 3, 8, 0);
//		if (nestedCascade.empty())
//			continue;
//		
//		if (eye_in_face) {
//			smallImgROI = smallImg(r);
//
//			std::vector<Rect> eyes;
//
//			////-- In each gray image, detect eyes
//			nestedCascade.detectMultiScale(smallImgROI, eyes, 1.1, 1, 0 | CASCADE_SCALE_IMAGE, Size(8, 8));
//
//			//imshow("small part", smallImgROI);
//			//waitKey(0);
//
//			for (size_t j = 0; j < eyes.size(); j++)
//			{
//				//code for draw samllROI images with circles
//				//Point center(faces[i].x + eyes[j].x + eyes[j].width*0.5, faces[i].y + eyes[j].y + eyes[j].height*0.5);
//				//int radius = cvRound((eyes[j].width + eyes[j].height)*0.25);
//				//circle(frame, center, radius, Scalar(30, 40, 200), 4, 8, 0);
//
//				//code for draw samllROI images with rectangles
//				rectangle(frame, Point(cvRound((faces[i].x + eyes[j].x) * scale), cvRound((faces[i].y + eyes[j].y) * scale)),
//					Point(cvRound((faces[i].x + eyes[j].x + eyes[j].width - 1) * scale),
//						cvRound((faces[i].y + eyes[j].y + eyes[j].height - 1) * scale)), Scalar(30,40,220), 3, 8, 0);
//
//			//	//code for draw gray scale images with circles
//			//	//Point center(eyes[j].x + eyes[j].width*0.5, eyes[j].y + eyes[j].height*0.5);
//			//	//int radius = cvRound((eyes[j].width + eyes[j].height)*0.25);
//			//	//circle(frame, center, radius, Scalar(255, 0, 0), 4, 8, 0);
//
//				//rectangle(frame, Point(cvRound(eyes[j].x * scale), cvRound(eyes[j].y * scale)),
//				//	Point(cvRound((eyes[j].x + eyes[j].width - 1) * scale),
//				//		cvRound((eyes[j].y + eyes[j].height - 1) * scale)), Scalar(30,40,220), 3, 8, 0);
//			}
//
//			//// Detection of eyes int the input image 
//			//nestedCascade.detectMultiScale(smallImgROI, nestedObjects, 1.1, 2,
//			//	0 | CASCADE_SCALE_IMAGE, Size(30, 30));
//
//			//imshow("small part", smallImgROI);
//			//waitKey(0);
//			//
//			//cout << "nestedObjects len: " << nestedObjects.size() << endl;
//			//// Draw circles around eyes 
//			//for (size_t j = 0; j < nestedObjects.size(); j++)
//			//{
//			//	Rect nr = nestedObjects[j];
//			//	center.x = cvRound((r.x + nr.x + nr.width*0.5)*scale);
//			//	center.y = cvRound((r.y + nr.y + nr.height*0.5)*scale);
//			//	radius = cvRound((nr.width + nr.height)*0.25*scale);
//			//	circle(frame, center, radius, color, 3, 8, 0);
//			//}
//		}
//	}
//
//	if (!eye_in_face) {
//		std::vector<Rect> eyes;
//		////-- In each gray image, detect eyes
//		nestedCascade.detectMultiScale(gray, eyes, 1.1, 1, 0 | CASCADE_SCALE_IMAGE, Size(8, 8));
//
//		////imshow("small part", smallImgROI);
//		////waitKey(0);
//
//		// Initialize the parameters
//		float confThreshold = 0.5f; // Confidence threshold
//		float nmsThreshold = 0.15f;  // Non-maximum suppression threshold
//		std::vector<int> indices;
//		std::vector<float> confidences(eyes.size(), 1.0f);
//		NMSBoxes(eyes, confidences, confThreshold, nmsThreshold, indices);
//
//		for (size_t j = 0; j < indices.size(); j++)
//		{
//			int idx = indices[j];
//			//code for draw gray scale images with circles
//			//Point center(eyes[j].x + eyes[j].width*0.5, eyes[j].y + eyes[j].height*0.5);
//			//int radius = cvRound((eyes[j].width + eyes[j].height)*0.25);
//			//circle(frame, center, radius, Scalar(255, 0, 0), 4, 8, 0);
//			
//			ofstream file(writefile, std::ios_base::app);
//			if (file.is_open()) {
//				file << "beta eye " << eyes[idx].x << " " << eyes[idx].y << " " << eyes[idx].x + eyes[idx].width << " " << eyes[idx].y + eyes[idx].height << "\n";
//				file.close();
//			}
//
//				//draw rectangles
//			rectangle(frame, Point(cvRound(eyes[idx].x * scale), cvRound(eyes[idx].y * scale)),
//				Point(cvRound((eyes[idx].x + eyes[idx].width - 1) * scale),
//					cvRound((eyes[idx].y + eyes[idx].height - 1) * scale)), Scalar(30, 40, 230), 3, 8, 0);
//		}
//	}
//
//	// Write the frame with the detection boxes
//	Mat detectedFrame;
//	frame.convertTo(detectedFrame, CV_8U);
//	imwrite(outputFile, detectedFrame);
//}


void predictImageHaarCascade(const string& filename, string class_name, string writefile, string channel) {
	// PreDefined trained XML classifiers with facial features 
	CascadeClassifier cascade;
	double scale = 1;

	//-- 1. Load the cascades
	if (class_name == "face") {
		if (!cascade.load("data\\haarcascade_frontalface_alt.xml")) {
			printf("--(!)Error loading\n"); return;
		}
	}
	else if (class_name == "eye") {
		if (!cascade.load("data\\haarcascade_eye.xml")) {
			printf("--(!)Error loading\n"); return;
		}
	}
	else if (class_name == "plate") {
		if (!cascade.load("data\\haarcascade_russian_plate_number.xml")) {
			printf("--(!)Error loading\n"); return;
		}
	}

	else if (class_name == "cat face") {
		if (!cascade.load("data\\lbpcascade_frontalcatface.xml")) {
			printf("--(!)Error loading\n"); return;
		}
	}

	else if (class_name == "eye 2") {
		if (!cascade.load("data\\frontalEyes35x16.xml")) {
			printf("--(!)Error loading\n"); return;
		}
	}

	else if (class_name == "nose") {
		if (!cascade.load("data\\Nariz.xml")) {
			printf("--(!)Error loading\n"); return;
		}
	}

	else if (class_name == "mouth") {
		if (!cascade.load("data\\Mouth.xml")) {
			printf("--(!)Error loading\n"); return;
		}
	}

	else if (class_name == "lowerbody") {
		if (!cascade.load("data\\haarcascade_lowerbody.xml")) {
			printf("--(!)Error loading\n"); return;
		}
	}

	else if (class_name == "right eye") {
		if (!cascade.load("data\\right_eye.xml")) {
			printf("--(!)Error loading\n"); return;
		}
	}

	else if (class_name == "head and shouler") {
		if (!cascade.load("data\\HS.xml")) {
			printf("--(!)Error loading\n"); return;
		}
	}

	else {
		std::cout << "No current classifier for " + class_name << std::endl;
		return;
	}

	

	
	// Open the image file
	ifstream ifile(filename);
	if (!ifile) throw("error");
	VideoCapture cap;
	cap.open(filename);
	string outputFile = filename.substr(0, filename.size() - 4) + "_haar_cascade_" + class_name + "_out_cpp.jpg";

	Mat frame;
	cap >> frame;

	//face detection part
	vector<Rect> objects;
	Mat gray;
	cvtColor(frame, gray, COLOR_BGR2GRAY);
	double fx = 1 / scale;

	// Resize the Grayscale Image  
	//resize(gray, smallImg, Size(), fx, fx, INTER_LINEAR);
	//equalizeHist(smallImg, smallImg);

	// Detect faces of different sizes using cascade classifier  
	cascade.detectMultiScale(gray, objects, 1.1,
		2, 0 | CASCADE_SCALE_IMAGE, Size(30, 30));


	//std::vector<int> rejectLevels;
	//std::vector<double> levelWeights;
	//cascade.detectMultiScale(gray, objects, rejectLevels, levelWeights, 1.1,
	//	2, 0 | CASCADE_SCALE_IMAGE, Size(30, 30),Size(), true);

	//if nothing found
	if (objects.size() == 0) {
		std::cout << "No object: " << class_name << " found!" << std::endl;
		return;
	}

	//normalize weights into range [0,1]
	//double weights_min = *std::min_element(levelWeights.begin(), levelWeights.end());
	//double weights_range = *std::max_element(levelWeights.begin(), levelWeights.end()) - weights_min + 1e-6;
	//for (int i = 0; i < levelWeights.size(); ++i) {
	//	levelWeights[i] = (levelWeights[i] - weights_min) / weights_range;
	//	std::cout << "levelWeights: " << i << " " << levelWeights[i] << std::endl;
	//}

		// Initialize the parameters
	float confThreshold = 0.5f; // Confidence threshold
	float nmsThreshold = 0.15f;  // Non-maximum suppression threshold
	std::vector<int> indices;
	std::vector<float> confidences(objects.size(), 1.0f);
	NMSBoxes(objects, confidences, confThreshold, nmsThreshold, indices);

	class_name.erase(remove(class_name.begin(), class_name.end(), ' '), class_name.end());

	for (size_t i = 0; i < indices.size(); i++)
		{
			Rect r = objects[indices[i]];
				
			ofstream file(writefile, std::ios_base::app);
			if (file.is_open()) {
				file << channel << " " << class_name << " " << confidences[i] << " " << 
					r.x << " " << r.y << " " << r.x + r.width << " " << r.y + r.height << "\n";
				file.close();
			}
		
			Scalar color = Scalar(255, 0, 0); // Color for Drawing tool 
			rectangle(frame, Point(cvRound(r.x*scale), cvRound(r.y*scale)),
				Point(cvRound((r.x + r.width - 1)*scale),
					cvRound((r.y + r.height - 1)*scale)), color, 3, 8, 0);
	}

	//imshow("objects: ", frame);
	//waitKey(0);
	// Write the frame with the detection boxes
	Mat detectedFrame;
	frame.convertTo(detectedFrame, CV_8U);
	imwrite(outputFile, detectedFrame);
}