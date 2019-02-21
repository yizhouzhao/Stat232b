#include "pch.h"
#include "haar_detection.h"

/** Global variables */
//String face_cascade_name = "haarcascade_frontalface_alt.xml";
//String eyes_cascade_name = "haarcascade_eye_tree_eyeglasses.xml";

string window_name = "Capture - Face detection";
RNG rng(12345);

void predictImageHaar(const string& filename, bool eye_in_face, string writefile) {
	// PreDefined trained XML classifiers with facial features 
	CascadeClassifier cascade, nestedCascade;
	double scale = 1;

	//-- 1. Load the cascades
	if (!cascade.load("data\\haarcascade_frontalface_alt.xml")) {
		printf("--(!)Error loading\n"); return;
	}
	if (!nestedCascade.load("data\\haarcascade_eye.xml")) {
			printf("--(!)Error loading\n"); return;
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
		
		ofstream file(writefile, std::ios_base::app);
		if (file.is_open()) {
			file << "alpha face " << r.x << " " << r.y << " " << r.x + r.width << " " << r.y + r.height << "\n";
			file.close();
		}

		Mat smallImgROI;
		vector<Rect> nestedObjects;
		Point center;
		Scalar color = Scalar(255, 0, 0); // Color for Drawing tool 
		int radius;

		double aspect_ratio = (double)r.width / r.height;
		if (false) //0.75 < aspect_ratio && aspect_ratio < 1.3)
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
		
		if (eye_in_face) {
			smallImgROI = smallImg(r);

			std::vector<Rect> eyes;

			////-- In each gray image, detect eyes
			nestedCascade.detectMultiScale(smallImgROI, eyes, 1.1, 1, 0 | CASCADE_SCALE_IMAGE, Size(8, 8));

			//imshow("small part", smallImgROI);
			//waitKey(0);

			for (size_t j = 0; j < eyes.size(); j++)
			{
				//code for draw samllROI images with circles
				//Point center(faces[i].x + eyes[j].x + eyes[j].width*0.5, faces[i].y + eyes[j].y + eyes[j].height*0.5);
				//int radius = cvRound((eyes[j].width + eyes[j].height)*0.25);
				//circle(frame, center, radius, Scalar(30, 40, 200), 4, 8, 0);

				//code for draw samllROI images with rectangles
				rectangle(frame, Point(cvRound((faces[i].x + eyes[j].x) * scale), cvRound((faces[i].y + eyes[j].y) * scale)),
					Point(cvRound((faces[i].x + eyes[j].x + eyes[j].width - 1) * scale),
						cvRound((faces[i].y + eyes[j].y + eyes[j].height - 1) * scale)), Scalar(30,40,220), 3, 8, 0);

			//	//code for draw gray scale images with circles
			//	//Point center(eyes[j].x + eyes[j].width*0.5, eyes[j].y + eyes[j].height*0.5);
			//	//int radius = cvRound((eyes[j].width + eyes[j].height)*0.25);
			//	//circle(frame, center, radius, Scalar(255, 0, 0), 4, 8, 0);

				//rectangle(frame, Point(cvRound(eyes[j].x * scale), cvRound(eyes[j].y * scale)),
				//	Point(cvRound((eyes[j].x + eyes[j].width - 1) * scale),
				//		cvRound((eyes[j].y + eyes[j].height - 1) * scale)), Scalar(30,40,220), 3, 8, 0);
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
	}

	if (!eye_in_face) {
		std::vector<Rect> eyes;
		////-- In each gray image, detect eyes
		nestedCascade.detectMultiScale(gray, eyes, 1.1, 1, 0 | CASCADE_SCALE_IMAGE, Size(8, 8));

		////imshow("small part", smallImgROI);
		////waitKey(0);

		// Initialize the parameters
		float confThreshold = 0.5f; // Confidence threshold
		float nmsThreshold = 0.15f;  // Non-maximum suppression threshold
		std::vector<int> indices;
		std::vector<float> confidences(eyes.size(), 1.0f);
		NMSBoxes(eyes, confidences, confThreshold, nmsThreshold, indices);

		for (size_t j = 0; j < indices.size(); j++)
		{
			int idx = indices[j];
			//code for draw gray scale images with circles
			//Point center(eyes[j].x + eyes[j].width*0.5, eyes[j].y + eyes[j].height*0.5);
			//int radius = cvRound((eyes[j].width + eyes[j].height)*0.25);
			//circle(frame, center, radius, Scalar(255, 0, 0), 4, 8, 0);
			
			ofstream file(writefile, std::ios_base::app);
			if (file.is_open()) {
				file << "beta eye " << eyes[idx].x << " " << eyes[idx].y << " " << eyes[idx].x + eyes[idx].width << " " << eyes[idx].y + eyes[idx].height << "\n";
				file.close();
			}

				//draw rectangles
			rectangle(frame, Point(cvRound(eyes[idx].x * scale), cvRound(eyes[idx].y * scale)),
				Point(cvRound((eyes[idx].x + eyes[idx].width - 1) * scale),
					cvRound((eyes[idx].y + eyes[idx].height - 1) * scale)), Scalar(30, 40, 230), 3, 8, 0);
		}
	}

	// Write the frame with the detection boxes
	Mat detectedFrame;
	frame.convertTo(detectedFrame, CV_8U);
	imwrite(outputFile, detectedFrame);
}