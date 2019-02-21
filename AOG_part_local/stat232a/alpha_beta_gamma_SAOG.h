#ifndef ALPHA_BETA_GAMMA_SAOG_H
#define	ALPHA_BETA_GAMMA_SAOG_H
#define MATH_PI 3.1415926

#include "stdlib.h"
#include "math.h"
#include <fstream>
#include <sstream>

#include "AOG.h"
#include "plot_aog.h"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/Imgproc.hpp>

using namespace std;
using namespace AOG_LIB;

AOG<std::string, std::vector<double>> AlphaBetaGammaSAOG(std::string alpha_name, double alpha_weight,
	std::vector<std::string> beta_names, double beta_weight,
	std::string gamma_name, double gamma_weight){
	
	//initialize rules
	std::vector<Symbolic_Rule<std::string>> rules;

	Symbolic_State<std::string> gamma(gamma_name, false);

	Symbolic_State<std::string> alpha(alpha_name, false);
	std::vector<double> alpha_attr = { alpha_weight, gamma_weight, beta_weight};


	std::vector<Symbolic_State<std::string>> beta;
	for (unsigned int i = 0; i < beta_names.size(); ++i) {
		Symbolic_State<std::string> beta_i(beta_names[i], true);
		beta.push_back(beta_i);
	}


	std::vector<Symbolic_State<std::string>> top = { alpha };
	Symbolic_Rule<std::string> gamma2alpha(gamma, top);
	Symbolic_Rule<std::string> alpha2beta(alpha, beta);

	rules.push_back(gamma2alpha);
	rules.push_back(alpha2beta);

	AOG<std::string, std::vector<double>> aog(rules);
	aog.SetRoot(gamma);

	aog.SetVertexAttribute(1, alpha_attr);

	return aog;
}




//Plot alpha_beta_gamma_saog
Mat PlotAlphaBetaGammaSAOG(AOG<std::string, std::vector<double>>& aog) {
	Mat frame = PlotAOG(aog);
	//cv::imwrite("C:\\Users\\Yizhou Zhao\\Desktop\\example_aog.jpg", frame);
	Size s = frame.size();
	int arrow_x = 150;
	int arrow_y = s.width / 2 - 30;
	Point arrow_begin(20, arrow_x);
	Point arrow_end(arrow_y, arrow_x);
	arrowedLine(frame, arrow_begin, arrow_end, Scalar(0, 0, 0), 3, 8, 0, 0.1);

	Mat dst = frame;
	copyMakeBorder(frame, dst, 0, 50, 0, 0, 0, Scalar(255, 255, 255));
	std::vector<double> alpha_attr = aog.GetVertexContent(1)->GetAttribute();


	putText(dst, std::to_string(alpha_attr[0]).substr(0, 4), Point(40, 140), 
		FONT_HERSHEY_SIMPLEX, 0.45, Scalar(200, 50, 55), 2, cv::LINE_8);

	putText(dst, std::to_string(alpha_attr[1]).substr(0, 4), Point(140, 100),
		FONT_HERSHEY_SIMPLEX, 0.45, Scalar(200, 50, 55), 2, cv::LINE_8);

	putText(dst, std::to_string(alpha_attr[2]).substr(0, 4), Point(140, 210),
		FONT_HERSHEY_SIMPLEX, 0.45, Scalar(200, 50, 55), 2, cv::LINE_8);


	return dst;
}

//determine a rect is overlapped with another
bool RectOverlap(Rect rect1, Rect rect2) {
	Point l1(rect1.x, rect1.y);
	Point r1(rect1.x + rect1.width, rect1.y + rect1.height);

	Point l2(rect2.x, rect2.y);
	Point r2(rect2.x + rect2.width, rect2.y + rect2.height);
	
	//std::cout << "Rect 1 Top " << l1.x << " " << l1.y << std::endl;
	//std::cout << "Rect 1 Bottom " << r1.x << " " << r1.y << std::endl;
	//std::cout << "Rect 2 Top " << l2.x << " " << l2.y << std::endl;
	//std::cout << "Rect 2 Bottom " << r2.x << " " << r2.y << std::endl;

	// If one rectangle is on left side of other 
	if (l1.x > r2.x || l2.x > r1.x)
		return false;

	// If one rectangle is above other 
	if (l1.y > r2.y || l2.y > r1.y)
		return false;

	return true;
}

std::vector<std::vector<Rect>> ReadRectFromFile(std::string filename) {
	std::vector<std::vector<Rect>> rect_list;
	ifstream file(filename);
	if (file.is_open()) {
		std::vector<Rect> alpha_rects;
		std::vector<Rect> beta_rects;
		std::vector<Rect> gamma_rects;
		std::string line;
		while (std::getline(file, line)) {
			//random delete lines for debug
			//if (rand() % 100 < 60)
			//	continue;
			std::stringstream linestream(line);
			std::string channel;
			std::getline(linestream, channel, ' ');
			std::string name;
			std::getline(linestream, name, ' ');
			int top_x;
			int top_y;
			int bottom_x;
			int bottom_y;
			linestream >> top_x >> top_y >> bottom_x >> bottom_y;
			//std::cout <<name<< " " << channel << " " << top_x << " " << top_y << " " << bottom_x << " " << bottom_y << " " << std::endl;
			Rect rect(top_x, top_y, bottom_x - top_x, bottom_y - top_y);
			if (channel == "alpha")
				alpha_rects.push_back(rect);
			else if (channel == "beta")
				beta_rects.push_back(rect);
			else if (channel == "gamma")
				gamma_rects.push_back(rect);
		}
		rect_list.emplace_back(alpha_rects);
		rect_list.emplace_back(beta_rects);
		rect_list.emplace_back(gamma_rects);
	}
	return rect_list;
}

double NormalDensity1D(double x, double mean, double variance) {
	return 1.0 / sqrt(2 * MATH_PI * variance) * exp(-0.5 * (x - mean) * (x - mean) / variance);
}

std::vector<double> LearnMeanAndVariance(std::vector<Rect> channel1, std::vector<Rect> channel2) {
	double _mean_center_x = 0; //acturally the top left corner
	double _mean_center_y = 0;
	double _var_center_x = 0;
	double _var_center_y = 0;
	double _mean_scale_x = 0;
	double _mean_scale_y = 0;
	double _var_scale_x = 0;
	double _var_scale_y = 0;

	int _count = 0;
	int _count_non_duplicate = 0;
	for (size_t i = 0; i < channel1.size(); ++i) {
		int activate = 0;
		for (size_t j = 0; j < channel2.size(); ++j) {
			if (RectOverlap(channel1[i], channel2[j])) {
				activate = 1;
				_count++;
				double center_x = (channel2[j].x - channel1[i].x + .0) / channel1[i].width;
				double center_y = (channel2[j].y - channel1[i].y + .0) / channel1[i].height;
				double scale_x = (channel2[j].width + .0) / channel1[i].width;
				double scale_y = (channel2[j].height + .0) / channel1[i].height;

				if (_count > 1) {
					_var_center_x = (_count - 2.0) / (_count - 1) * _var_center_x +
						1.0 / _count * (center_x - _mean_center_x) * (center_x - _mean_center_x);
					_var_center_y = (_count - 2.0) / (_count - 1) * _var_center_y +
						1.0 / _count * (center_y - _mean_center_y) * (center_y - _mean_center_y);
					_var_scale_x = (_count - 2.0) / (_count - 1) * _var_scale_x +
						1.0 / _count * (scale_x - _mean_scale_x) * (scale_x - _mean_scale_x);
					_var_scale_y = (_count - 2.0) / (_count - 1) * _var_scale_y +
						1.0 / _count * (scale_y - _mean_scale_y) * (scale_y - _mean_scale_y);
				}

				_mean_center_x = (center_x + (_count - 1) * _mean_center_x) / _count;
				_mean_center_y = (center_y + (_count - 1) * _mean_center_y) / _count;
				_mean_scale_x = (scale_x + (_count - 1) * _mean_scale_x) / _count;
				_mean_scale_y = (scale_y + (_count - 1) * _mean_scale_y) / _count;
			}
		}
		_count_non_duplicate += activate;
	}
	
	return std::vector<double>({ double(_count), double(_count_non_duplicate), 
		_mean_center_x, _mean_center_y, _var_center_x, _var_center_y,
		_mean_scale_x, _mean_scale_y , _var_scale_x , _var_scale_y });
}


AOG<std::string, std::vector<double>> LearnAlphaBetaGammaSAOG(std::vector<Rect>& alpha_rects, std::vector<Rect>& beta_rects, std::vector<Rect>& gamma_rects) {
	std::vector<double> g2a_info = LearnMeanAndVariance(gamma_rects, alpha_rects);

	double gamma_to_alpha = g2a_info[1] / gamma_rects.size(); //????????????????????????????????????
	std::cout << "\n gamma_overlap_alpha_count " << gamma_to_alpha << std::endl;

	std::cout << "gamma:: position " << g2a_info[2] << " " << g2a_info[3] << " " << g2a_info[4] << " " << g2a_info[5] << std::endl;
	std::cout << "gamma:: scale " << g2a_info[6] << " " << g2a_info[7] << " " << g2a_info[8] << " " << g2a_info[9] << std::endl;

	std::vector<double> a2b_info = LearnMeanAndVariance(alpha_rects, beta_rects);

	std::cout << "\n beta_overlap_alpha_count " << a2b_info[0] << std::endl;
	double alpha_to_beta = a2b_info[1] / alpha_rects.size(); //????????????????????????????????????
	std::cout << "\n beta_overlap_alpha_count " << alpha_to_beta << std::endl;

	std::cout << "beta:: position " << a2b_info[2] << " " << a2b_info[3] << " " << a2b_info[4] << " " << a2b_info[5] << std::endl;
	std::cout << "beta:: scale " << a2b_info[6] << " " << a2b_info[7] << " " << a2b_info[8] << " " << a2b_info[9] << std::endl;

	std::vector<double> g2b_info = LearnMeanAndVariance(gamma_rects, beta_rects);

	
	double gamma_to_beta = g2b_info[1] / gamma_rects.size(); //????????????????????????????????????
	std::cout << "\n beta_overlap_gamma_count " << gamma_to_beta << std::endl;

	std::cout << "gamma->beta:: position " << g2b_info[2] << " " << g2b_info[3] << " " << g2b_info[4] << " " << g2b_info[5] << std::endl;
	std::cout << "gamma->beta:: scale " << g2b_info[6] << " " << g2b_info[7] << " " << g2b_info[8] << " " << g2b_info[9] << std::endl;


	std::vector<double> alpha_scores(alpha_rects.size(), 0);
	std::vector<double> beta_scores(beta_rects.size(), 0);
	std::vector<double> gamma_scores(gamma_rects.size(), 0);

	std::vector<std::vector<Rect>> group_objects;
	std::vector<double> group_scores;

	double g2a_scale = NormalDensity1D(0, 0, g2a_info[4]) * NormalDensity1D(0, 0, g2a_info[5]) *
		NormalDensity1D(0, 0, g2a_info[8]) * NormalDensity1D(0, 0, g2a_info[9]);
	double g2b_scale = NormalDensity1D(0, 0, g2b_info[4]) * NormalDensity1D(0, 0, g2b_info[5]) *
		NormalDensity1D(0, 0, g2b_info[8]) * NormalDensity1D(0, 0, g2b_info[9]);
	double a2b_scale = NormalDensity1D(0, 0, a2b_info[4]) * NormalDensity1D(0, 0, a2b_info[5]) *
		NormalDensity1D(0, 0, a2b_info[8]) * NormalDensity1D(0, 0, a2b_info[9]);

	std::cout << "scalerL: " << g2a_scale << " " << g2b_scale << " " << a2b_scale << " " << std::endl;

	for (size_t i = 0; i < gamma_rects.size(); ++i) {
		for (size_t j = 0; j < alpha_rects.size(); ++j) {
			if (RectOverlap(gamma_rects[i], alpha_rects[j])) {
				double center_x = (alpha_rects[j].x - gamma_rects[i].x + .0) / gamma_rects[i].width;
				double center_y = (alpha_rects[j].y - gamma_rects[i].y + .0) / gamma_rects[i].height;
				double scale_x = (alpha_rects[j].width + .0) / gamma_rects[i].width;
				double scale_y = (alpha_rects[j].height + .0) / gamma_rects[i].height;

				double add_score = NormalDensity1D(center_x, g2a_info[2], g2a_info[4]) *
					NormalDensity1D(center_y, g2a_info[3], g2a_info[5]) *
					NormalDensity1D(scale_x, g2a_info[6], g2a_info[8]) *
					NormalDensity1D(scale_y, g2a_info[7], g2a_info[9]);

				add_score /= g2a_scale;

				//std::cout << "add score: " << i << " " << j << " " << " " << add_score << std::endl;
				gamma_scores[i] += add_score * gamma_to_alpha;
				alpha_scores[j] += add_score * gamma_to_alpha;
			}
		}
	}

	for (size_t j = 0; j < alpha_rects.size(); ++j) {
		for (size_t k = 0; k < beta_rects.size(); ++k) {
			if (RectOverlap(alpha_rects[j], beta_rects[k])) {
				double center_x = (beta_rects[k].x - alpha_rects[j].x + .0) / alpha_rects[j].width;
				double center_y = (beta_rects[k].y - alpha_rects[j].y + .0) / alpha_rects[j].height;
				double scale_x = (beta_rects[k].width + .0) / alpha_rects[j].width;
				double scale_y = (beta_rects[k].height + .0) / alpha_rects[j].height;

				double add_score = NormalDensity1D(center_x, a2b_info[2], a2b_info[4]) *
					NormalDensity1D(center_y, a2b_info[3], a2b_info[5]) *
					NormalDensity1D(scale_x, a2b_info[6], a2b_info[8]) *
					NormalDensity1D(scale_y, a2b_info[7], a2b_info[9]);

				add_score /= a2b_scale;
				//std::cout << "add score: j k " << j  << " " << k << " " << " " << add_score << std::endl;
				beta_scores[k] += add_score  * alpha_to_beta;
				alpha_scores[j] += add_score  * alpha_to_beta;
			}
		}
	}
	
	for (size_t i = 0; i < gamma_rects.size(); ++i) {
		for (size_t k = 0; k < beta_rects.size(); ++k) {
			if (RectOverlap(gamma_rects[i], beta_rects[k])) {
				double center_x = (beta_rects[k].x - gamma_rects[i].x + .0) / gamma_rects[i].width;
				double center_y = (beta_rects[k].y - gamma_rects[i].y + .0) / gamma_rects[i].height;
				double scale_x = (beta_rects[k].width + .0) / gamma_rects[i].width;
				double scale_y = (beta_rects[k].height + .0) / gamma_rects[i].height;

				double add_score = NormalDensity1D(center_x, g2b_info[2], g2b_info[4]) *
					NormalDensity1D(center_y, g2b_info[3], g2b_info[5]) *
					NormalDensity1D(scale_x, g2b_info[6], g2b_info[8]) *
					NormalDensity1D(scale_y, g2b_info[7], g2b_info[9]);

				add_score /= g2b_scale;
				//std::cout << "add score: i, k " << i << " " << k << " " << " " << add_score << std::endl;
				beta_scores[k] += add_score * gamma_to_beta;
				gamma_scores[i] += add_score * gamma_to_beta;
			}
		}
	}

	Mat frame = Mat::zeros(3024, 4032, CV_8UC3);
	frame = cv::Scalar(255, 255, 255);

	for (int i = 0; i < gamma_rects.size(); i++) {
		std::cout << "gamma " << i << " score " << gamma_scores[i] << std::endl;
		Scalar color(180, 100, 20);
		rectangle(frame, gamma_rects[i], color, 4, 8, 0);
	}

	for (int j = 0; j < alpha_rects.size(); j++) {
		std::cout << "alpha " << j << " score " << alpha_scores[j] << std::endl;
		if (alpha_scores[j] < 0.5)
			continue;
		Scalar color(20, 180, 100);
		rectangle(frame, alpha_rects[j], color, 4, 8, 0);
	}

	for (int k = 0; k < beta_rects.size(); k++) {
		std::cout << "beta " << k << " score " << beta_scores[k] << std::endl;
		if (beta_scores[k] < 1) continue;
		Scalar color(100, 20, 180);
		rectangle(frame, beta_rects[k], color, 4, 8, 0);
	}

	cv::resize(frame, frame, Size(frame.cols / 4, frame.rows / 4));
	imshow("frame", frame);
	waitKey(0);

	//imwrite("C:\\Users\\Yizhou Zhao\\Desktop\\pic\\independent_no_background_filtered.jpg", frame);

	std::vector<Rect> reconstruced_gamma_rects;
	std::vector<Rect> reconstruced_alpha_rects;
	std::vector<Rect> reconstruced_beta_rects;
	for (int j = 0; j < alpha_rects.size(); j++) {
		if (alpha_scores[j] < 0.5)
			continue;

	}

	AOG<std::string, std::vector<double>> aog = AlphaBetaGammaSAOG("alpha", 1, { "beta" }, alpha_to_beta, "gamma", gamma_to_alpha);
	return aog;
}



#endif // !ALPHA_BETA_GAMMA_SAOG_H