#include "pch.h"
#include "random_generator.h"

void generateRectanglesOnFrame(Mat frame, std::string filename, int num) {
	int image_height = frame.rows;
	int image_width = frame.cols;
	//Mat frame = Mat::zeros(image_height, image_width, CV_8UC3);
	//frame = cv::Scalar(0);
	std::random_device device;
	std::mt19937 generator(device());

	for (int n = 0; n < num; ++n) {
		//generate body part
		std::uniform_int_distribution<int> top_x_distribution(1, image_height);
		std::uniform_int_distribution<int> top_y_distribution(1, image_width);

		std::normal_distribution<double> rect_height(image_height / 3.0, 10);
		std::normal_distribution<double> rect_width(image_width / 3.0, 10);

		int top_x = top_x_distribution(generator);
		int top_y = top_y_distribution(generator);
		int bottom_x = std::min(top_x + std::max((int)rect_height(generator), 10), image_height - 1);
		int bottom_y = std::min(top_y + std::max((int)rect_width(generator), 10), image_width - 1);

		//std::cout << top_x << " " << top_y << " " << bottom_x << " " << bottom_y << " " << std::endl;
		Scalar color0 = Scalar(255, 75, 21);

		// Draw a rectangle ( 5th argument is not -ve)
		rectangle(frame, Point(top_y, top_x), Point(bottom_y, bottom_x), color0, 1, 8);
		//imshow("Image1", frame);
		//cv::waitKey(0);

		//generate face part
		std::uniform_int_distribution<int> top_x_1_distribution(top_x, bottom_x);
		std::uniform_int_distribution<int> top_y_1_distribution(top_y, bottom_y);

		std::normal_distribution<double> rect_height_1((bottom_x - top_x) / 2.0, 5);
		std::normal_distribution<double> rect_width_1((bottom_y - top_y) / 2.0, 5);

		int top_x_1 = top_x_1_distribution(generator);
		int top_y_1 = top_y_1_distribution(generator);
		int bottom_x_1 = std::min(top_x_1 + std::max((int)rect_height_1(generator), 8), bottom_x - 1);
		int bottom_y_1 = std::min(top_y_1 + std::max((int)rect_width_1(generator), 8), bottom_y - 1);

		Scalar color1 = Scalar(55, 75, 210);
		rectangle(frame, Point(top_y_1, top_x_1), Point(bottom_y_1, bottom_x_1), color1, 1, 8);
		//imshow("Image2", frame);
		//cv::waitKey(0);

		//generate eye parts
		std::uniform_int_distribution<int> top_x_21_distribution(top_x_1, bottom_x_1);
		std::uniform_int_distribution<int> top_y_21_distribution(top_y_1, bottom_y_1);

		std::normal_distribution<double> rect_height_21((bottom_x_1 - top_x_1) / 4.0, 3);
		std::normal_distribution<double> rect_width_21((bottom_y_1 - top_y_1) / 4.0, 3);

		int top_x_21 = top_x_21_distribution(generator);
		int top_y_21 = top_y_21_distribution(generator);
		int bottom_x_21 = std::min(top_x_21 + std::max((int)rect_height_21(generator), 3), bottom_x_1 - 1);
		int bottom_y_21 = std::min(top_y_21 + std::max((int)rect_width_21(generator), 3), bottom_y_1 - 1);

		Scalar color21 = Scalar(155, 20, 110);
		rectangle(frame, Point(top_y_21, top_x_21), Point(bottom_y_21, bottom_x_21), color21, 1, 8);
		//imshow("Image3", frame);
		//cv::waitKey(0);


		std::uniform_int_distribution<int> top_x_22_distribution(top_x_1, bottom_x_1);
		std::uniform_int_distribution<int> top_y_22_distribution(top_y_1, bottom_y_1);

		std::normal_distribution<double> rect_height_22((bottom_x_1 - top_x_1) / 4.0, 3);
		std::normal_distribution<double> rect_width_22((bottom_y_1 - top_y_1) / 4.0, 3);

		int top_x_22 = top_x_22_distribution(generator);
		int top_y_22 = top_y_22_distribution(generator);
		int bottom_x_22 = std::min(top_x_22 + std::max((int)rect_height_22(generator), 3), bottom_x_1 - 1);
		int bottom_y_22 = std::min(top_y_22 + std::max((int)rect_width_22(generator), 3), bottom_y_1 - 1);

		Scalar color22 = Scalar(155, 20, 110);
		rectangle(frame, Point(top_y_22, top_x_22), Point(bottom_y_22, bottom_x_22), color22, 1, 8);


		// Draw a filled rectangle ( 5th argument is -ve)
		//rectangle(frame, Point(115, 120), Point(170, 150), Scalar(100, 155, 25), -1, 8);
		//imshow("Image2", frame);

		if (n == 0) {
			ofstream file(filename);
			if (file.is_open()) file.close();
		}

		ofstream file(filename, std::ios_base::app);
		if (file.is_open()) {
			file << "gamma body " << top_x << " " << top_y << " " << bottom_x << " " << bottom_y << "\n";
			file << "alpha face " << top_x_1 << " " << top_y_1 << " " << bottom_x_1 << " " << bottom_y_1 << "\n";
			file << "beta eye " << top_x_21 << " " << top_y_21 << " " << bottom_x_21 << " " << bottom_y_21 << "\n";
			file << "beta eye " << top_x_22 << " " << top_y_22 << " " << bottom_x_22 << " " << bottom_y_22 << "\n";
			file.close();
		}
	}
}

void writeRandomRectToFile(Mat frame, int num, std::string file_name) {

}