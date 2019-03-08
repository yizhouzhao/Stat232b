// yolo_for_AOG.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <opencv2/opencv.hpp>
using namespace cv;
using namespace dnn;
#include "yolo_detection.h"
#include "haar_detection.h"
#include "random_generator.h"
#include "circle_detection.h"
#include "video_to_image.h"

int main()
{
    std::cout << "Hello World!\n"; 
	//vector<string> names_list = getNamesOfClasses("data\\coco.names");
	//std::cout << names_list[10] << std::endl;
	
	string filename = "C:\\Users\\Yizhou Zhao\\Desktop\\pic\\2.jpeg"; //"C:\\Users\\Yizhou Zhao\\Desktop\\pic\\Test 1.jpeg";
	//predictImageYolov3(filename,"person", "C:\\Users\\Yizhou Zhao\\Desktop\\pic\\classroom_33.txt"); //"C:\\Users\\Yizhou Zhao\\Desktop\\pic\\n_classroom_2.txt");

	//predictImageHaar(filename, false, "C:\\Users\\Yizhou Zhao\\Desktop\\pic\\n_classroom_2.txt");
	predictImageHaarCascade(filename, "head and shouler", "C:\\Users\\Yizhou Zhao\\Desktop\\pic\\classroom_N2.txt", "gamma");
	predictImageHaarCascade(filename, "face", "C:\\Users\\Yizhou Zhao\\Desktop\\pic\\classroom_N2.txt", "alpha");
	predictImageHaarCascade(filename, "nose", "C:\\Users\\Yizhou Zhao\\Desktop\\pic\\classroom_N2.txt", "beta");

	//Mat frame = Mat::zeros(600, 800, CV_8UC3);
	//generateRectanglesOnFrame(frame, "test.txt", 10);
	//imshow("pic", frame);
	//imwrite("C:\\Users\\Yizhou Zhao\\Desktop\\pic\\demo.jpg", frame);
	//waitKey(0);

	//predictCircle("C:\\Users\\Yizhou Zhao\\Desktop\\pic2\\cat1.jpg");
	//Video2Image("C:\\Users\\Yizhou Zhao\\Desktop\\pic2\\2.mp4", "C:\\Users\\Yizhou Zhao\\Desktop\\pic2\\video_box\\", 5);

	//VideoAddBox("C:\\Users\\Yizhou Zhao\\Desktop\\pic2\\2.mp4", "C:\\Users\\Yizhou Zhao\\Desktop\\pic2\\video_box\\", 5);

	return 0; 
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
