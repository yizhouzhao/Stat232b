#ifndef ALPHA_BETA_GAMMA_SAOG_H
#define	ALPHA_BETA_GAMMA_SAOG_H
#define MATH_PI 3.1415926

#include "stdlib.h"
#include "math.h"
#include <fstream>
#include <sstream>
#include <algorithm> 


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

double RectOverlapArea(Rect rect1, Rect rect2) {
	Point l1(rect1.x, rect1.y);
	Point r1(rect1.x + rect1.width, rect1.y + rect1.height);

	Point l2(rect2.x, rect2.y);
	Point r2(rect2.x + rect2.width, rect2.y + rect2.height);

	int areaI = (min(r1.x, r2.x) -
		max(l1.x, l2.x)) *
		(min(r1.y, r2.y) -
			max(l1.y, l2.y));

	return (areaI + .0) / min(rect1.area(), rect2.area());
}

//determine a rect is inside another
bool RectInside(Rect rect1, Rect rect2) {
	if (rect1.x >= rect2.x &&
		rect1.y >= rect2.y &&
		rect2.x + rect2.width >= rect1.x + rect1.width &&
		rect2.y + rect2.height >= rect1.y + rect1.height
		)
		return true;

	else if (rect2.x >= rect1.x &&
		rect2.y >= rect1.y &&
		rect1.x + rect1.width >= rect2.x + rect2.width &&
		rect1.y + rect1.height >= rect2.y + rect2.height
		)
		return true;

	return false;
}

//std::vector<std::vector<Rect>> ReadRectFromFile(std::string filename) {
//	std::vector<std::vector<Rect>> rect_list;
//	std::ifstream file(filename);
//	if (file.is_open()) {
//		std::vector<Rect> alpha_rects;
//		std::vector<Rect> beta_rects;
//		std::vector<Rect> gamma_rects;
//		std::string line;
//		while (std::getline(file, line)) {
//			random delete lines for debug
//			if (rand() % 100 < 60)
//				continue;
//			std::stringstream linestream(line);
//			std::string channel;
//			std::getline(linestream, channel, ' ');
//			std::string name;
//			std::getline(linestream, name, ' ');
//			int top_x;
//			int top_y;
//			int bottom_x;
//			int bottom_y;
//			linestream >> top_x >> top_y >> bottom_x >> bottom_y;
//			std::cout <<name<< " " << channel << " " << top_x << " " << top_y << " " << bottom_x << " " << bottom_y << " " << std::endl;
//			Rect rect(top_x, top_y, bottom_x - top_x, bottom_y - top_y);
//			if (channel == "alpha")
//				alpha_rects.push_back(rect);
//			else if (channel == "beta")
//				beta_rects.push_back(rect);
//			else if (channel == "gamma")
//				gamma_rects.push_back(rect);
//		}
//		rect_list.emplace_back(alpha_rects);
//		rect_list.emplace_back(beta_rects);
//		rect_list.emplace_back(gamma_rects);
//	}
//	return rect_list;
//}

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

	std::unordered_map<std::string, double> information_map;

	int _count = 0;
	int _count_non_duplicate = 0;
	for (size_t i = 0; i < channel1.size(); ++i) {
		int activate = 0;
		for (size_t j = 0; j < channel2.size(); ++j) {
			if (RectInside(channel1[i], channel2[j])) {
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

	information_map["channel size"] = double(channel1.size());
	information_map["overlap count"] = double(_count_non_duplicate);
	information_map["top x mean"] = _mean_center_x;
	information_map["top y mean"] = _mean_center_y;
	information_map["scale x mean"] = _mean_scale_x;
	information_map["scale y mean"] = _mean_scale_y;
	information_map["top x variance"] = _var_center_x;
	information_map["top y variance"] = _var_center_y;
	information_map["scale x variance"] = _var_scale_x;
	information_map["scale y variance"] = _var_scale_y;

	return std::vector<double>({ double(channel1.size()), double(_count_non_duplicate),
		_mean_center_x, _mean_center_y, _var_center_x, _var_center_y,
		_mean_scale_x, _mean_scale_y , _var_scale_x , _var_scale_y });
}

AOG<std::string, std::vector<double>> LearnAndParseAlphaBetaGammaSAOG(const std::vector<Rect>& alpha_rects, const std::vector<Rect>& beta_rects,const std::vector<Rect>& gamma_rects,
	const std::vector<double>& alpha_confidence = std::vector<double>(),
	const std::vector<double>& beta_confidence = std::vector<double>(), 
	const std::vector<double>& gamma_confidence = std::vector<double>()) {
	//treshold for activate channel
	const double alpha_threshold = 1.2;
	const double beta_threshold = 1.2;
	const double gamma_threshold = 1.2;
	
	//learn the position and scale information from gamma to alpha channel
	std::vector<double> g2a_info = LearnMeanAndVariance(gamma_rects, alpha_rects);
	double gamma_to_alpha = g2a_info[1] / gamma_rects.size(); //????????????????????????????????????
	//std::cout << "\n gamma_overlap_alpha_count " << gamma_to_alpha << std::endl;

	//std::cout << "gamma:: position " << g2a_info[2] << " " << g2a_info[3] << " " << g2a_info[4] << " " << g2a_info[5] << std::endl;
	//std::cout << "gamma:: scale " << g2a_info[6] << " " << g2a_info[7] << " " << g2a_info[8] << " " << g2a_info[9] << std::endl;

	//learn the position and scale information from alpha to beta channel
	std::vector<double> a2b_info = LearnMeanAndVariance(alpha_rects, beta_rects);
	//std::cout << "\n beta_overlap_alpha_count " << a2b_info[0] << std::endl;
	double alpha_to_beta = a2b_info[1] / alpha_rects.size(); //????????????????????????????????????
	//std::cout << "\n beta_overlap_alpha_count " << alpha_to_beta << std::endl;
	//std::cout << "beta:: position " << a2b_info[2] << " " << a2b_info[3] << " " << a2b_info[4] << " " << a2b_info[5] << std::endl;
	//std::cout << "beta:: scale " << a2b_info[6] << " " << a2b_info[7] << " " << a2b_info[8] << " " << a2b_info[9] << std::endl;

	//learn the position and scale information from gamma to beta channel
	std::vector<double> g2b_info = LearnMeanAndVariance(gamma_rects, beta_rects);
	double gamma_to_beta = gamma_to_alpha * alpha_to_beta; //g2b_info[1] / gamma_rects.size(); //????????????????????????????????????
	//std::cout << "\n beta_overlap_gamma_count " << gamma_to_beta << std::endl;
	//std::cout << "gamma->beta:: position " << g2b_info[2] << " " << g2b_info[3] << " " << g2b_info[4] << " " << g2b_info[5] << std::endl;
	//std::cout << "gamma->beta:: scale " << g2b_info[6] << " " << g2b_info[7] << " " << g2b_info[8] << " " << g2b_info[9] << std::endl;


	std::vector<double> alpha_scores; 
	if (alpha_confidence.size() == 0) {
		alpha_scores = std::vector<double>(alpha_rects.size(), 1);
	}
	else {
		CV_Assert(alpha_confidence.size() == alpha_rects.size());
		alpha_scores = alpha_confidence;
	}

	std::vector<double> beta_scores;
	if (beta_confidence.size() == 0) {
		beta_scores = std::vector<double>(beta_rects.size(), 1);
	}
	else {
		CV_Assert(beta_confidence.size() == beta_rects.size());
		beta_scores = beta_confidence;
	}

	std::vector<double> gamma_scores;
	if (gamma_confidence.size() == 0) {
		gamma_scores = std::vector<double>(gamma_rects.size(), 1);
	}
	else {
		CV_Assert(gamma_confidence.size() == gamma_rects.size());
		gamma_scores = gamma_confidence;
	}


	std::map<std::pair<int, int>, double> gamma2alpha_score_record;
	std::map<std::pair<int, int>, double> gamma2beta_score_record;
	std::map<std::pair<int, int>, double> alpha2beta_score_record;

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
			if (RectOverlapArea(gamma_rects[i], alpha_rects[j]) > 0.9 && RectOverlap(gamma_rects[i], alpha_rects[j])) {
				double center_x = (alpha_rects[j].x - gamma_rects[i].x + .0) / gamma_rects[i].width;
				double center_y = (alpha_rects[j].y - gamma_rects[i].y + .0) / gamma_rects[i].height;
				double scale_x = (alpha_rects[j].width + .0) / gamma_rects[i].width;
				double scale_y = (alpha_rects[j].height + .0) / gamma_rects[i].height;

				double add_score = NormalDensity1D(center_x, g2a_info[2], g2a_info[4]) *
					NormalDensity1D(center_y, g2a_info[3], g2a_info[5]) *
					NormalDensity1D(scale_x, g2a_info[6], g2a_info[8]) *
					NormalDensity1D(scale_y, g2a_info[7], g2a_info[9]);

				add_score /= g2a_scale;
				gamma2alpha_score_record[{i, j}] = add_score;
				//std::cout << gamma2alpha_score_record[{i, j}]  << "!!!!!!!add score: " << i << " " << j << " " << " " << add_score << std::endl;
				gamma_scores[i] += add_score * gamma_to_alpha;
				alpha_scores[j] += add_score * gamma_to_alpha;
			}
		}
	}

	//for (auto i = 0; i < gamma_rects.size(); i++) {
	//	for (size_t j = 0; j < alpha_rects.size(); ++j) {
	//		std::cout << i << " " << j << " "<< gamma2alpha_score_record[{i, j}] << std::endl;
	//	}
	//}

	for (size_t j = 0; j < alpha_rects.size(); ++j) {
		for (size_t k = 0; k < beta_rects.size(); ++k) {
			if (RectInside(alpha_rects[j], beta_rects[k])) {
				double center_x = (beta_rects[k].x - alpha_rects[j].x + .0) / alpha_rects[j].width;
				double center_y = (beta_rects[k].y - alpha_rects[j].y + .0) / alpha_rects[j].height;
				double scale_x = (beta_rects[k].width + .0) / alpha_rects[j].width;
				double scale_y = (beta_rects[k].height + .0) / alpha_rects[j].height;

				double add_score = NormalDensity1D(center_x, a2b_info[2], a2b_info[4]) *
					NormalDensity1D(center_y, a2b_info[3], a2b_info[5]) *
					NormalDensity1D(scale_x, a2b_info[6], a2b_info[8]) *
					NormalDensity1D(scale_y, a2b_info[7], a2b_info[9]);

				add_score /= a2b_scale;
				alpha2beta_score_record[{j, k}] = add_score;
				//std::cout << "add score: j k " << j  << " " << k << " " << " " << add_score << std::endl;
				beta_scores[k] += add_score  * alpha_to_beta;
				alpha_scores[j] += add_score  * alpha_to_beta;
			}
		}
	}
	
	for (size_t i = 0; i < gamma_rects.size(); ++i) {
		for (size_t k = 0; k < beta_rects.size(); ++k) {
			if (RectInside(gamma_rects[i], beta_rects[k])) {
				double center_x = (beta_rects[k].x - gamma_rects[i].x + .0) / gamma_rects[i].width;
				double center_y = (beta_rects[k].y - gamma_rects[i].y + .0) / gamma_rects[i].height;
				double scale_x = (beta_rects[k].width + .0) / gamma_rects[i].width;
				double scale_y = (beta_rects[k].height + .0) / gamma_rects[i].height;

				double add_score = NormalDensity1D(center_x, g2b_info[2], g2b_info[4]) *
					NormalDensity1D(center_y, g2b_info[3], g2b_info[5]) *
					NormalDensity1D(scale_x, g2b_info[6], g2b_info[8]) *
					NormalDensity1D(scale_y, g2b_info[7], g2b_info[9]);

				add_score /= g2b_scale;
				gamma2beta_score_record[{i, k}] = add_score;
				//std::cout << "add score: i, k " << i << " " << k << " " << " " << add_score << std::endl;
				beta_scores[k] += add_score * gamma_to_beta;
				gamma_scores[i] += add_score * gamma_to_beta;
			}
		}
	}

	// visualize part
	Mat frame = Mat::zeros(3024, 4032, CV_8UC3);
	frame = cv::Scalar(255, 255, 255);
	for (int i = 0; i < gamma_rects.size(); i++) {
		std::cout << "gamma " << i << " score " << gamma_scores[i] << std::endl;
		Scalar color(180, 100, 20);
		rectangle(frame, gamma_rects[i], color, 4, 8, 0);
	}

	for (int j = 0; j < alpha_rects.size(); j++) {
		std::cout << "alpha " << j << " score " << alpha_scores[j] << std::endl;
		//if (alpha_scores[j] < alpha_threshold) continue;
		Scalar color(20, 180, 100);
		rectangle(frame, alpha_rects[j], color, 4, 8, 0);
	}

	for (int k = 0; k < beta_rects.size(); k++) {
		std::cout << "beta " << k << " score " << beta_scores[k] << std::endl;
		//if (beta_scores[k] < beta_threshold) continue;
		Scalar color(100, 20, 180);
		rectangle(frame, beta_rects[k], color, 4, 8, 0);
	}

	cv::resize(frame, frame, Size(frame.cols / 4, frame.rows / 4));
	imshow("frame", frame);
	//waitKey(0);

	//imwrite("C:\\Users\\Yizhou Zhao\\Desktop\\pic\\test_independent_no_background.jpg", frame);

	std::vector<Rect> reconstruced_gamma_rects;
	std::vector<Rect> reconstruced_alpha_rects;
	std::vector<Rect> reconstruced_beta_rects;

	std::vector<bool> visited_alpha_rects(alpha_rects.size(), false);
	std::vector<bool> visited_beta_rects(beta_rects.size(), false);
	std::vector<bool> visited_gamma_rects(gamma_rects.size(), false);

	//Reconstruct BETA & GAMMA from ALPHA
	for (size_t j = 0; j < alpha_rects.size(); j++) {
		if (alpha_scores[j] < alpha_threshold)
			continue;

		visited_alpha_rects[j] = true;

		int visited_gamma_index = -1;
		double visited_gamma_max_score = 0;

		for (size_t i = 0; i < gamma_rects.size(); ++i) {
			if ((RectOverlapArea(gamma_rects[i], alpha_rects[j]) > 0.9 && RectOverlap(gamma_rects[i], alpha_rects[j])) && (!visited_gamma_rects[i])) {
				if (visited_gamma_max_score < gamma2alpha_score_record[{i, j}]) {
					/*std::cout << "overlap or inside: " << i << " " << j << " " << "\n"
						<<
						RectOverlapArea(gamma_rects[i], alpha_rects[j]) << " " << RectInside(gamma_rects[i], alpha_rects[j]) << std::endl;*/
					//std::cout <<"scores: " <<gamma2alpha_score_record[{i, j}] << " vs " << gamma_scores[i] << " \n";
					visited_gamma_max_score = gamma2alpha_score_record[{i, j}]; //gamma_scores[i]; //gamma2alpha_score_record[{i, j}];
					visited_gamma_index = i;
				}
			}
		}
		//std::cout << "chosed: " << visited_gamma_index << " " << j << " " << visited_gamma_max_score << "\n";

		int visited_beta_index = -1;
		double visited_beta_max_score = 0;

		for (int k = 0; k < beta_rects.size(); ++k) {
			if (RectInside(beta_rects[k], alpha_rects[j]) && !visited_beta_rects[k]) {
				if (visited_beta_max_score < beta_scores[k]) {
					visited_beta_max_score = beta_scores[k];
					visited_beta_index = k;
				}
			}
		}
		
		//reconstruct alpha channel
		reconstruced_alpha_rects.push_back(alpha_rects[j]);

		//reconstruct gamma channel
		if (visited_gamma_index >= 0) {
			visited_gamma_rects[visited_gamma_index] = true;
			reconstruced_gamma_rects.push_back(gamma_rects[visited_gamma_index]);
		}
		else {
			//return std::vector<double>({ double(_count), double(_count_non_duplicate),
			//_mean_center_x, _mean_center_y, _var_center_x, _var_center_y,
			//_mean_scale_x, _mean_scale_y , _var_scale_x , _var_scale_y });

			double scale_x = alpha_rects[j].width / g2a_info[6];
			double scale_y = alpha_rects[j].height / g2a_info[7];
			double center_x = alpha_rects[j].x - g2a_info[2] * scale_x;
			double center_y = alpha_rects[j].y - g2a_info[3] * scale_y;
			
			reconstruced_gamma_rects.push_back(Rect(cvRound(center_x), cvRound(center_y),
				cvRound(scale_x), cvRound(scale_y)));
		}

		//reconstruct beta channel
		if (visited_beta_index >= 0) {
			visited_beta_rects[visited_beta_index] = true;
			reconstruced_beta_rects.push_back(beta_rects[visited_beta_index]);
		}
		else {
			//return std::vector<double>({ double(_count), double(_count_non_duplicate),
			//_mean_center_x, _mean_center_y, _var_center_x, _var_center_y,
			//_mean_scale_x, _mean_scale_y , _var_scale_x , _var_scale_y });

			double scale_x = alpha_rects[j].width * a2b_info[6];
			double scale_y = alpha_rects[j].height * a2b_info[7];
			double center_x = alpha_rects[j].x + a2b_info[2] * alpha_rects[j].width;
			double center_y = alpha_rects[j].y + a2b_info[3] * alpha_rects[j].height;

			reconstruced_beta_rects.push_back(Rect(cvRound(center_x), cvRound(center_y),
				cvRound(scale_x), cvRound(scale_y)));
		}
	}

	Mat frame2 = Mat::zeros(3024, 4032, CV_8UC3);
	//Mat frame2 = Mat::zeros(900, 1200, CV_8UC3);
	frame2 = cv::Scalar(255, 255, 255);
	//CV_Assert(reconstruced_gamma_rects.size() == reconstruced_beta_rects.size());
	//CV_Assert(reconstruced_alpha_rects.size() == reconstruced_beta_rects.size());

	//Reconstruct ALPHA from BETA_GAMMA
	for (size_t i = 0; i < gamma_rects.size(); i++) {
		if (gamma_scores[i] < gamma_threshold || visited_gamma_rects[i])
			continue;

		//std::cout << "gamma activated: " << i << " " << gamma_scores[i] << std::endl;

		Mat t_frame = Mat::zeros(3024, 4032, CV_8UC3);
		t_frame = cv::Scalar(255, 255, 255);

		rectangle(t_frame, gamma_rects[i], Scalar(0, 0, 0), 5, 8, 0);

		visited_gamma_rects[i] = true;
		double scale_x = gamma_rects[i].width * g2a_info[6];
		double scale_y = gamma_rects[i].height * g2a_info[7];
		double center_x = gamma_rects[i].x + g2a_info[2] * gamma_rects[i].width;
		double center_y = gamma_rects[i].y + g2a_info[3] * gamma_rects[i].height;

		Rect alpha_reconstructed(cvRound(center_x), cvRound(center_y),
			cvRound(scale_x), cvRound(scale_y));

		int visited_alpha_index = -1;
		double visited_alpha_max_score = 0;
		for (size_t j = 0; j < alpha_rects.size(); j++) {
			if (!visited_alpha_rects[j] && RectOverlap(alpha_reconstructed, alpha_rects[j])) {
				if (visited_alpha_max_score < alpha_scores[j]) {
					visited_alpha_index = j;
					visited_alpha_max_score = alpha_scores[j];
				}
			}
		}

		if (visited_alpha_index >= 0) {
			visited_alpha_rects[visited_alpha_index] = 1;
			alpha_reconstructed = alpha_rects[visited_alpha_index];
		}

		//rectangle(t_frame, alpha_reconstructed, Scalar(0, 0, 100), 5, 8, 0);
		//cv::resize(t_frame, t_frame, Size(t_frame.cols / 4, t_frame.rows / 4));
		//imshow("t_frame", t_frame);
		//waitKey(0);

		int visited_beta_index = -1;
		double visited_beta_max_score = 0;
		for (size_t k = 0; k < beta_rects.size(); ++k) {
			if (RectInside(alpha_reconstructed, beta_rects[k]) && !visited_beta_rects[k]) { // || 
				if (visited_beta_max_score < beta_scores[k]) {
					visited_beta_index = k;
					visited_beta_max_score = beta_scores[k];
				}
			}
		}

		if (visited_beta_index >= 0) {
			visited_beta_rects[visited_beta_index] = 1;
			reconstruced_beta_rects.push_back(beta_rects[visited_beta_index]);
			reconstruced_gamma_rects.push_back(gamma_rects[i]);
			reconstruced_alpha_rects.push_back(alpha_reconstructed);
		}
	}

	for (int i = 0; i < reconstruced_gamma_rects.size(); i++) {
		//std::cout << "gamma " << i << " score " << gamma_scores[i] << std::endl;
		Scalar color(180, 100, 20);
		rectangle(frame2, reconstruced_gamma_rects[i], color, 4, 8, 0);
	}

	for (int j = 0; j < reconstruced_alpha_rects.size(); j++) {
		//std::cout << "alpha " << j << " score " << alpha_scores[j] << std::endl;
		Scalar color(20, 180, 100);
		rectangle(frame2, reconstruced_alpha_rects[j], color, 4, 8, 0);
	}

	for (int k = 0; k < reconstruced_beta_rects.size(); k++) {
		//std::cout << "beta " << k << " score " << beta_scores[k] << std::endl;
		Scalar color(100, 20, 180);
		rectangle(frame2, reconstruced_beta_rects[k], color, 4, 8, 0);
	}

	cv::resize(frame2, frame2, Size(frame2.cols / 4, frame2.rows / 4));
	imshow("frame reconstructed", frame2);
	waitKey(0);
	//imwrite("C:\\Users\\Yizhou Zhao\\Desktop\\pic\\test1_independent_no_background_reconstrued.jpg", frame2);

	AOG<std::string, std::vector<double>> aog = AlphaBetaGammaSAOG("alpha", 1, { "beta" }, alpha_to_beta, "gamma", gamma_to_alpha);
	return aog;
}

std::unordered_map<std::string, double> LearnAttributes(std::vector<Rect> channel1, std::vector<Rect> channel2) {
	double _mean_center_x = 0; //acturally the top left corner
	double _mean_center_y = 0;
	double _var_center_x = 0;
	double _var_center_y = 0;
	double _mean_scale_x = 0;
	double _mean_scale_y = 0;
	double _var_scale_x = 0;
	double _var_scale_y = 0;

	std::unordered_map<std::string, double> information_map;

	int _count = 0;
	int _count_non_duplicate = 0;
	for (size_t i = 0; i < channel1.size(); ++i) {
		int activate = 0;
		for (size_t j = 0; j < channel2.size(); ++j) {
			if (RectInside(channel1[i], channel2[j])) {
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

	information_map["channel size"] = double(channel1.size());
	information_map["overlap count"] = double(_count_non_duplicate);
	information_map["top x mean"] = _mean_center_x;
	information_map["top y mean"] = _mean_center_y;
	information_map["scale x mean"] = _mean_scale_x;
	information_map["scale y mean"] = _mean_scale_y;
	information_map["top x variance"] = _var_center_x;
	information_map["top y variance"] = _var_center_y;
	information_map["scale x variance"] = _var_scale_x;
	information_map["scale y variance"] = _var_scale_y;

	return information_map;
}

AOG<std::string, std::unordered_map<std::string, double>> LearnAlphaBetaGammaSAOG(const std::vector<Rect>& alpha_rects, const std::vector<Rect>& beta_rects, const std::vector<Rect>& gamma_rects,
	const std::vector<double>& alpha_confidence = std::vector<double>(),
	const std::vector<double>& beta_confidence = std::vector<double>(),
	const std::vector<double>& gamma_confidence = std::vector<double>()) {

	//learn the position and scale information from gamma to alpha channel
	std::unordered_map<std::string, double> g2a_info = LearnAttributes(gamma_rects, alpha_rects);
	double gamma_to_alpha = g2a_info["overlap count"] / g2a_info["channel size"]; //????????????????????????????????????
	//std::cout << "\n gamma_overlap_alpha_count " << gamma_to_alpha << std::endl;

	//std::cout << "gamma:: position " << g2a_info[2] << " " << g2a_info[3] << " " << g2a_info[4] << " " << g2a_info[5] << std::endl;
	//std::cout << "gamma:: scale " << g2a_info[6] << " " << g2a_info[7] << " " << g2a_info[8] << " " << g2a_info[9] << std::endl;

	//learn the position and scale information from alpha to beta channel
	std::unordered_map<std::string, double> a2b_info = LearnAttributes(alpha_rects, beta_rects);
	//std::cout << "\n beta_overlap_alpha_count " << a2b_info[0] << std::endl;
	double alpha_to_beta = a2b_info["overlap count"] / a2b_info["channel size"]; //????????????????????????????????????
	//std::cout << "\n beta_overlap_alpha_count " << alpha_to_beta << std::endl;
	//std::cout << "beta:: position " << a2b_info[2] << " " << a2b_info[3] << " " << a2b_info[4] << " " << a2b_info[5] << std::endl;
	//std::cout << "beta:: scale " << a2b_info[6] << " " << a2b_info[7] << " " << a2b_info[8] << " " << a2b_info[9] << std::endl;

	//learn the position and scale information from gamma to beta channel
	std::unordered_map<std::string, double> g2b_info = LearnAttributes(gamma_rects, beta_rects);
	double gamma_to_beta = gamma_to_alpha * alpha_to_beta; //g2b_info[1] / gamma_rects.size(); //????????????????????????????????????
	//std::cout << "\n beta_overlap_gamma_count " << gamma_to_beta << std::endl;
	//std::cout << "gamma->beta:: position " << g2b_info[2] << " " << g2b_info[3] << " " << g2b_info[4] << " " << g2b_info[5] << std::endl;
	//std::cout << "gamma->beta:: scale " << g2b_info[6] << " " << g2b_info[7] << " " << g2b_info[8] << " " << g2b_info[9] << std::endl;

	std::vector<Symbolic_Rule<std::string>> rules;

	Symbolic_State<std::string> gamma("gamma", false);
	Symbolic_State<std::string> alpha("alpha", false);
	Symbolic_State<std::string> beta("beta", false);

	std::vector<Symbolic_State<std::string>> top = { alpha };
	std::vector<Symbolic_State<std::string>> bottom = { beta };

	Symbolic_Rule<std::string> gamma2alpha(gamma, top);
	Symbolic_Rule<std::string> alpha2beta(alpha, bottom);

	rules.push_back(gamma2alpha);
	rules.push_back(alpha2beta);

	AOG<std::string, std::unordered_map<std::string, double>> aog(rules);
	aog.SetRoot(gamma);

	VertexId alpha_id = aog.GetVertexIdByState(alpha);
	VertexId beta_id = aog.GetVertexIdByState(beta);
	VertexId gamma_id = aog.GetVertexIdByState(gamma);

	//std::cout << "alpha_id: " << alpha_id << " beta_id: " << beta_id << " gamma_id: " << gamma_id << std::endl;
	aog.SetVertexAttribute(alpha_id, a2b_info);
	aog.SetVertexAttribute(gamma_id, g2a_info);
	aog.SetVertexAttribute(beta_id, g2b_info);

	//unordered_map<VertexId, double> g2a_weights = aog.GetOutEdgeWeights(gamma_id, false);
	//g2a_weights[alpha_id] = gamma_to_alpha;

	//unordered_map<VertexId, double> a2b_weights = aog.GetOutEdgeWeights(alpha_id, false);
	//g2a_weights[alpha_id] = alpha_to_beta;

	//Mat frame = PlotAOG(aog);
	//cv::imshow("frame", frame);
	//cv::waitKey(0);

	return aog;
}

void ParseAlphaBetaGammaSAOG(const std::vector<Rect>& alpha_rects, const std::vector<Rect>& beta_rects, const std::vector<Rect>& gamma_rects, 
	const AOG<std::string, std::unordered_map<std::string, double>>& aog, std::string writefile,
	const std::vector<double>& alpha_confidence = std::vector<double>(),
	const std::vector<double>& beta_confidence = std::vector<double>(),
	const std::vector<double>& gamma_confidence = std::vector<double>()) {
	//treshold for activate channel
	const double alpha_threshold = 1.2;
	const double beta_threshold = 1.2;
	const double gamma_threshold = 1.2;

	std::vector<double> alpha_scores;
	if (alpha_confidence.size() == 0) {
		alpha_scores = std::vector<double>(alpha_rects.size(), 1);
	}
	else {
		CV_Assert(alpha_confidence.size() == alpha_rects.size());
		alpha_scores = alpha_confidence;
	}

	std::vector<double> beta_scores;
	if (beta_confidence.size() == 0) {
		beta_scores = std::vector<double>(beta_rects.size(), 1);
	}
	else {
		CV_Assert(beta_confidence.size() == beta_rects.size());
		beta_scores = beta_confidence;
	}

	std::vector<double> gamma_scores;
	if (gamma_confidence.size() == 0) {
		gamma_scores = std::vector<double>(gamma_rects.size(), 1);
	}
	else {
		CV_Assert(gamma_confidence.size() == gamma_rects.size());
		gamma_scores = gamma_confidence;
	}

	std::map<std::pair<int, int>, double> gamma2alpha_score_record;
	std::map<std::pair<int, int>, double> gamma2beta_score_record;
	std::map<std::pair<int, int>, double> alpha2beta_score_record;

	std::vector<std::vector<Rect>> group_objects;
	std::vector<double> group_scores;

	std::unordered_map<std::string, double> g2a_info = aog.GetVertexContent(0)->GetAttribute();
	std::unordered_map<std::string, double> a2b_info = aog.GetVertexContent(1)->GetAttribute();
	std::unordered_map<std::string, double> g2b_info = aog.GetVertexContent(2)->GetAttribute();

	double g2a_scale = NormalDensity1D(0, 0, g2a_info["top x variance"]) * NormalDensity1D(0, 0, g2a_info["top y variance"]) *
		NormalDensity1D(0, 0, g2a_info["scale x variance"]) * NormalDensity1D(0, 0, g2a_info["scale y variance"]);
	double g2b_scale = NormalDensity1D(0, 0, g2b_info["top x variance"]) * NormalDensity1D(0, 0, g2b_info["top y variance"]) *
		NormalDensity1D(0, 0, g2b_info["scale x variance"]) * NormalDensity1D(0, 0, g2b_info["scale y variance"]);
	double a2b_scale = NormalDensity1D(0, 0, a2b_info["top x variance"]) * NormalDensity1D(0, 0, a2b_info["top y variance"]) *
		NormalDensity1D(0, 0, a2b_info["scale x variance"]) * NormalDensity1D(0, 0, a2b_info["scale y variance"]);

	std::cout << "scalerL: " << g2a_scale << " " << g2b_scale << " " << a2b_scale << " " << std::endl;

	double gamma_to_alpha = g2a_info["overlap count"] / g2a_info["channel size"];
	for (size_t i = 0; i < gamma_rects.size(); ++i) {
		for (size_t j = 0; j < alpha_rects.size(); ++j) {
			if (RectOverlapArea(gamma_rects[i], alpha_rects[j]) > 0.9 && RectOverlap(gamma_rects[i], alpha_rects[j])) {
				double center_x = (alpha_rects[j].x - gamma_rects[i].x + .0) / gamma_rects[i].width;
				double center_y = (alpha_rects[j].y - gamma_rects[i].y + .0) / gamma_rects[i].height;
				double scale_x = (alpha_rects[j].width + .0) / gamma_rects[i].width;
				double scale_y = (alpha_rects[j].height + .0) / gamma_rects[i].height;

				double add_score = NormalDensity1D(center_x, g2a_info["top x mean"], g2a_info["top x variance"]) *
					NormalDensity1D(center_y, g2a_info["top y mean"], g2a_info["top y variance"]) *
					NormalDensity1D(scale_x, g2a_info["scale x mean"], g2a_info["scale x variance"]) *
					NormalDensity1D(scale_y, g2a_info["scale y mean"], g2a_info["scale y variance"]);

				add_score /= g2a_scale;
				gamma2alpha_score_record[{i, j}] = add_score;
				//std::cout << gamma2alpha_score_record[{i, j}]  << "!!!!!!!add score: " << i << " " << j << " " << " " << add_score << std::endl;
				gamma_scores[i] += add_score * gamma_to_alpha;
				alpha_scores[j] += add_score * gamma_to_alpha;
			}
		}
	}

	//for (auto i = 0; i < gamma_rects.size(); i++) {
	//	for (size_t j = 0; j < alpha_rects.size(); ++j) {
	//		std::cout << i << " " << j << " "<< gamma2alpha_score_record[{i, j}] << std::endl;
	//	}
	//}

	double alpha_to_beta = a2b_info["overlap count"] / a2b_info["channel size"];
	for (size_t j = 0; j < alpha_rects.size(); ++j) {
		for (size_t k = 0; k < beta_rects.size(); ++k) {
			if (RectInside(alpha_rects[j], beta_rects[k])) {
				double center_x = (beta_rects[k].x - alpha_rects[j].x + .0) / alpha_rects[j].width;
				double center_y = (beta_rects[k].y - alpha_rects[j].y + .0) / alpha_rects[j].height;
				double scale_x = (beta_rects[k].width + .0) / alpha_rects[j].width;
				double scale_y = (beta_rects[k].height + .0) / alpha_rects[j].height;

				double add_score = NormalDensity1D(center_x, a2b_info["top x mean"], a2b_info["top x variance"]) *
					NormalDensity1D(center_y, a2b_info["top y mean"], a2b_info["top y variance"]) *
					NormalDensity1D(scale_x, a2b_info["scale x mean"], a2b_info["scale x variance"]) *
					NormalDensity1D(scale_y, a2b_info["scale y mean"], a2b_info["scale y variance"]);

				add_score /= a2b_scale;
				alpha2beta_score_record[{j, k}] = add_score;
				//std::cout << "add score: j k " << j  << " " << k << " " << " " << add_score << std::endl;
				beta_scores[k] += add_score * alpha_to_beta;
				alpha_scores[j] += add_score * alpha_to_beta;
			}
		}
	}

	double gamma_to_beta = g2b_info["overlap count"] / g2b_info["channel size"];
	for (size_t i = 0; i < gamma_rects.size(); ++i) {
		for (size_t k = 0; k < beta_rects.size(); ++k) {
			if (RectInside(gamma_rects[i], beta_rects[k])) {
				double center_x = (beta_rects[k].x - gamma_rects[i].x + .0) / gamma_rects[i].width;
				double center_y = (beta_rects[k].y - gamma_rects[i].y + .0) / gamma_rects[i].height;
				double scale_x = (beta_rects[k].width + .0) / gamma_rects[i].width;
				double scale_y = (beta_rects[k].height + .0) / gamma_rects[i].height;

				double add_score = NormalDensity1D(center_x, g2b_info["top x mean"], g2b_info["top x variance"]) *
					NormalDensity1D(center_y, g2b_info["top y mean"], g2b_info["top y variance"]) *
					NormalDensity1D(scale_x, g2b_info["scale x mean"], g2b_info["scale x variance"]) *
					NormalDensity1D(scale_y, g2b_info["scale y mean"], g2b_info["scale y variance"]);

				add_score /= g2b_scale;
				gamma2beta_score_record[{i, k}] = add_score;
				//std::cout << "add score: i, k " << i << " " << k << " " << " " << add_score << std::endl;
				beta_scores[k] += add_score * gamma_to_beta;
				gamma_scores[i] += add_score * gamma_to_beta;
			}
		}
	}

	// visualize part
	Mat frame = Mat::zeros(3024, 4032, CV_8UC3);
	//Mat frame = Mat::zeros(900, 1200, CV_8UC3);
	frame = cv::Scalar(255, 255, 255);
	for (int i = 0; i < gamma_rects.size(); i++) {
		std::cout << "gamma " << i << " score " << gamma_scores[i] << std::endl;
		Scalar color(180, 100, 20);
		rectangle(frame, gamma_rects[i], color, 4, 8, 0);
	}

	for (int j = 0; j < alpha_rects.size(); j++) {
		std::cout << "alpha " << j << " score " << alpha_scores[j] << std::endl;
		//if (alpha_scores[j] < alpha_threshold) continue;
		Scalar color(20, 180, 100);
		rectangle(frame, alpha_rects[j], color, 4, 8, 0);
	}

	for (int k = 0; k < beta_rects.size(); k++) {
		std::cout << "beta " << k << " score " << beta_scores[k] << std::endl;
		//if (beta_scores[k] < beta_threshold) continue;
		Scalar color(100, 20, 180);
		rectangle(frame, beta_rects[k], color, 4, 8, 0);
	}

	cv::resize(frame, frame, Size(frame.cols / 4, frame.rows / 4));
	imshow("frame", frame);
	waitKey(0);

	imwrite("C:\\Users\\Yizhou Zhao\\Desktop\\pic\\test_3_independent_no_background.jpg", frame);

	std::vector<Rect> reconstruced_gamma_rects;
	std::vector<Rect> reconstruced_alpha_rects;
	std::vector<Rect> reconstruced_beta_rects;

	std::vector<bool> visited_alpha_rects(alpha_rects.size(), false);
	std::vector<bool> visited_beta_rects(beta_rects.size(), false);
	std::vector<bool> visited_gamma_rects(gamma_rects.size(), false);

	//Reconstruct BETA & GAMMA from ALPHA
	for (size_t j = 0; j < alpha_rects.size(); j++) {
		if (alpha_scores[j] < alpha_threshold)
			continue;

		visited_alpha_rects[j] = true;

		int visited_gamma_index = -1;
		double visited_gamma_max_score = 0;

		for (size_t i = 0; i < gamma_rects.size(); ++i) {
			if ((RectOverlapArea(gamma_rects[i], alpha_rects[j]) > 0.9 && RectOverlap(gamma_rects[i], alpha_rects[j])) && (!visited_gamma_rects[i])) {
				if (visited_gamma_max_score < gamma2alpha_score_record[{i, j}]) {
					/*std::cout << "overlap or inside: " << i << " " << j << " " << "\n"
						<<
						RectOverlapArea(gamma_rects[i], alpha_rects[j]) << " " << RectInside(gamma_rects[i], alpha_rects[j]) << std::endl;*/
						//std::cout <<"scores: " <<gamma2alpha_score_record[{i, j}] << " vs " << gamma_scores[i] << " \n";
					visited_gamma_max_score = gamma2alpha_score_record[{i, j}]; //gamma_scores[i]; //gamma2alpha_score_record[{i, j}];
					visited_gamma_index = i;
				}
			}
		}
		//std::cout << "chosed: " << visited_gamma_index << " " << j << " " << visited_gamma_max_score << "\n";

		int visited_beta_index = -1;
		double visited_beta_max_score = 0;

		for (int k = 0; k < beta_rects.size(); ++k) {
			if (RectInside(beta_rects[k], alpha_rects[j]) && !visited_beta_rects[k]) {
				if (visited_beta_max_score < beta_scores[k]) {
					visited_beta_max_score = beta_scores[k];
					visited_beta_index = k;
				}
			}
		}

		//reconstruct alpha channel
		reconstruced_alpha_rects.push_back(alpha_rects[j]);

		//reconstruct gamma channel
		if (visited_gamma_index >= 0) {
			visited_gamma_rects[visited_gamma_index] = true;
			reconstruced_gamma_rects.push_back(gamma_rects[visited_gamma_index]);
		}
		else {
			//return std::vector<double>({ double(_count), double(_count_non_duplicate),
			//_mean_center_x, _mean_center_y, _var_center_x, _var_center_y,
			//_mean_scale_x, _mean_scale_y , _var_scale_x , _var_scale_y });

			double scale_x = alpha_rects[j].width / g2a_info["scale x mean"];
			double scale_y = alpha_rects[j].height / g2a_info["scale y mean"];
			double center_x = alpha_rects[j].x - g2a_info["top x mean"] * scale_x;
			double center_y = alpha_rects[j].y - g2a_info["top y mean"] * scale_y;

			reconstruced_gamma_rects.push_back(Rect(cvRound(center_x), cvRound(center_y),
				cvRound(scale_x), cvRound(scale_y)));
		}

		//reconstruct beta channel
		if (visited_beta_index >= 0) {
			visited_beta_rects[visited_beta_index] = true;
			reconstruced_beta_rects.push_back(beta_rects[visited_beta_index]);
		}
		else {
			//return std::vector<double>({ double(_count), double(_count_non_duplicate),
			//_mean_center_x, _mean_center_y, _var_center_x, _var_center_y,
			//_mean_scale_x, _mean_scale_y , _var_scale_x , _var_scale_y });

			double scale_x = alpha_rects[j].width * a2b_info["scale x mean"];
			double scale_y = alpha_rects[j].height * a2b_info["scale y mean"];
			double center_x = alpha_rects[j].x + a2b_info["top x mean"] * alpha_rects[j].width;
			double center_y = alpha_rects[j].y + a2b_info["top y mean"] * alpha_rects[j].height;

			reconstruced_beta_rects.push_back(Rect(cvRound(center_x), cvRound(center_y),
				cvRound(scale_x), cvRound(scale_y)));
		}
	}

	Mat frame2 = Mat::zeros(3024, 4032, CV_8UC3);
	//Mat frame2 = Mat::zeros(900, 1200, CV_8UC3);
	frame2 = cv::Scalar(255, 255, 255);
	//CV_Assert(reconstruced_gamma_rects.size() == reconstruced_beta_rects.size());
	//CV_Assert(reconstruced_alpha_rects.size() == reconstruced_beta_rects.size());

	//Reconstruct ALPHA from BETA_GAMMA
	for (size_t i = 0; i < gamma_rects.size(); i++) {
		if (gamma_scores[i] < gamma_threshold || visited_gamma_rects[i])
			continue;

		//std::cout << "gamma activated: " << i << " " << gamma_scores[i] << std::endl;

		Mat t_frame = Mat::zeros(3024, 4032, CV_8UC3);
		t_frame = cv::Scalar(255, 255, 255);

		rectangle(t_frame, gamma_rects[i], Scalar(0, 0, 0), 5, 8, 0);

		visited_gamma_rects[i] = true;
		double scale_x = gamma_rects[i].width * g2a_info["scale x mean"];
		double scale_y = gamma_rects[i].height * g2a_info["scale y mean"];
		double center_x = gamma_rects[i].x + g2a_info["top x mean"] * gamma_rects[i].width;
		double center_y = gamma_rects[i].y + g2a_info["top y mean"] * gamma_rects[i].height;

		Rect alpha_reconstructed(cvRound(center_x), cvRound(center_y),
			cvRound(scale_x), cvRound(scale_y));

		int visited_alpha_index = -1;
		double visited_alpha_max_score = 0;
		for (size_t j = 0; j < alpha_rects.size(); j++) {
			if (!visited_alpha_rects[j] && RectOverlap(alpha_reconstructed, alpha_rects[j])) {
				if (visited_alpha_max_score < alpha_scores[j]) {
					visited_alpha_index = j;
					visited_alpha_max_score = alpha_scores[j];
				}
			}
		}

		if (visited_alpha_index >= 0) {
			visited_alpha_rects[visited_alpha_index] = 1;
			alpha_reconstructed = alpha_rects[visited_alpha_index];
		}

		//rectangle(t_frame, alpha_reconstructed, Scalar(0, 0, 100), 5, 8, 0);
		//cv::resize(t_frame, t_frame, Size(t_frame.cols / 4, t_frame.rows / 4));
		//imshow("t_frame", t_frame);
		//waitKey(0);

		int visited_beta_index = -1;
		double visited_beta_max_score = 0;
		for (size_t k = 0; k < beta_rects.size(); ++k) {
			if (RectInside(alpha_reconstructed, beta_rects[k]) && !visited_beta_rects[k]) { // || 
				if (visited_beta_max_score < beta_scores[k]) {
					visited_beta_index = k;
					visited_beta_max_score = beta_scores[k];
				}
			}
		}

		if (visited_beta_index >= 0) {
			visited_beta_rects[visited_beta_index] = 1;
			reconstruced_beta_rects.push_back(beta_rects[visited_beta_index]);
			reconstruced_gamma_rects.push_back(gamma_rects[i]);
			reconstruced_alpha_rects.push_back(alpha_reconstructed);
		}
	}

	for (int i = 0; i < reconstruced_gamma_rects.size(); i++) {
		//std::cout << "gamma " << i << " score " << gamma_scores[i] << std::endl;
		Scalar color(180, 100, 20);
		rectangle(frame2, reconstruced_gamma_rects[i], color, 4, 8, 0);
	}

	for (int j = 0; j < reconstruced_alpha_rects.size(); j++) {
		//std::cout << "alpha " << j << " score " << alpha_scores[j] << std::endl;
		Scalar color(20, 180, 100);
		rectangle(frame2, reconstruced_alpha_rects[j], color, 4, 8, 0);
	}

	for (int k = 0; k < reconstruced_beta_rects.size(); k++) {
		//std::cout << "beta " << k << " score " << beta_scores[k] << std::endl;
		Scalar color(100, 20, 180);
		rectangle(frame2, reconstruced_beta_rects[k], color, 4, 8, 0);
	}

	cv::resize(frame2, frame2, Size(frame2.cols / 4, frame2.rows / 4));
	imshow("frame2", frame2);
	waitKey(0);

	//imwrite("C:\\Users\\Yizhou Zhao\\Desktop\\pic\\test_3_independent_no_background_reconstrued.jpg", frame2);

	std::ofstream file(writefile);
	if (file.is_open()) {
		for (size_t i = 0; i < reconstruced_gamma_rects.size(); ++i) {
			stringstream ss;
			file << "gamma gamma 1 " << reconstruced_gamma_rects[i].x << " " <<
				reconstruced_gamma_rects[i].y << " "
				<< reconstruced_gamma_rects[i].x + reconstruced_gamma_rects[i].width << " "
				<< reconstruced_gamma_rects[i].y + reconstruced_gamma_rects[i].height << "\n";
			file << "alpha alpha 1 " << reconstruced_alpha_rects[i].x << " " <<
				reconstruced_alpha_rects[i].y << " "
				<< reconstruced_alpha_rects[i].x + reconstruced_alpha_rects[i].width << " "
				<< reconstruced_alpha_rects[i].y + reconstruced_alpha_rects[i].height << "\n";
			file << "beta beta 1 " << reconstruced_beta_rects[i].x << " " <<
				reconstruced_beta_rects[i].y << " "
				<< reconstruced_beta_rects[i].x + reconstruced_beta_rects[i].width << " "
				<< reconstruced_beta_rects[i].y + reconstruced_beta_rects[i].height << "\n";
		}
	}
}

//AOG<std::string, std::vector<double>> LearnAndParseAlphaBetaGammaSAOG(const std::vector<std::vector<Rect>>& multi_channels, const std::vector<std::vector<double>>& multi_confidences) {
//	return LearnAndParseAlphaBetaGammaSAOG(multi_channels[0], multi_channels[1], multi_channels[2],
//		multi_confidences[0], multi_confidences[1], multi_confidences[2]);
//}

void LearnAndParseAlphaBetaGammaSAOG2(const std::vector<std::vector<Rect>>& multi_channels, const std::vector<std::vector<double>>& multi_confidences, std::string writefile) {
	AOG<std::string, std::unordered_map<std::string, double>> aog = LearnAlphaBetaGammaSAOG(multi_channels[0], multi_channels[1], multi_channels[2],
		multi_confidences[0], multi_confidences[1], multi_confidences[2]);
	ParseAlphaBetaGammaSAOG(multi_channels[0], multi_channels[1], multi_channels[2], aog, writefile,
		multi_confidences[0], multi_confidences[1], multi_confidences[2]);
}

//std::vector<double>
//std::unordered_map<std::string, double>
AOG<std::string, std::unordered_map<std::string, double>> LearnAlphaBetaGammaAOGFromFile(std::string filename) {
	std::vector<std::vector<Rect>> rect_list;
	std::vector<std::vector<double>> confidence_list;
	std::ifstream file(filename);
	if (file.is_open()) {
		std::vector<Rect> alpha_rects;
		std::vector<Rect> beta_rects;
		std::vector<Rect> gamma_rects;

		std::vector<double> alpha_confidence;
		std::vector<double> beta_confidence;
		std::vector<double> gamma_confidence;

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
			double confidence;
			int top_x;
			int top_y;
			int bottom_x;
			int bottom_y;
			linestream >> confidence >> top_x >> top_y >> bottom_x >> bottom_y;
			//std::cout <<name<< " " << channel << " " << top_x << " " << top_y << " " << bottom_x << " " << bottom_y << " " << std::endl;
			Rect rect(top_x, top_y, bottom_x - top_x, bottom_y - top_y);
			if (channel == "alpha") {
				alpha_rects.push_back(rect);
				alpha_confidence.push_back(confidence);
			}
				
			else if (channel == "beta") {
				beta_rects.push_back(rect);
				beta_confidence.push_back(confidence);
			}
				
			else if (channel == "gamma") {
				gamma_rects.push_back(rect);
				gamma_confidence.push_back(confidence);
			}
				
		}
		rect_list.emplace_back(alpha_rects);
		rect_list.emplace_back(beta_rects);
		rect_list.emplace_back(gamma_rects);

		confidence_list.emplace_back(alpha_confidence);
		confidence_list.emplace_back(beta_confidence);
		confidence_list.emplace_back(gamma_confidence);

	}
	//AOG<std::string, std::vector<double>> aog_1 = LearnAndParseAlphaBetaGammaSAOG(rect_list, confidence_list);
	//LearnAndParseAlphaBetaGammaSAOG2(rect_list, confidence_list, writefile);
	AOG<std::string, std::unordered_map<std::string, double>> aog = LearnAlphaBetaGammaSAOG(rect_list[0], rect_list[1], rect_list[2],
		confidence_list[0], confidence_list[1], confidence_list[2]);

	return aog;
}

void LearnAndParseAlphaBetaGammaAOGFromFile(std::string input_file, std::string target_file, std::string write_file) {
	AOG<std::string, std::unordered_map<std::string, double>> aog = LearnAlphaBetaGammaAOGFromFile(input_file);

	std::vector<std::vector<Rect>> out_rect_list;
	std::vector<std::vector<double>> out_confidence_list;
	std::ifstream outfile(target_file);
	if (outfile.is_open()) {
		std::vector<Rect> alpha_rects;
		std::vector<Rect> beta_rects;
		std::vector<Rect> gamma_rects;

		std::vector<double> alpha_confidence;
		std::vector<double> beta_confidence;
		std::vector<double> gamma_confidence;

		std::string line;
		while (std::getline(outfile, line)) {
			//random delete lines for debug
			//if (rand() % 100 < 60)
			//	continue;
			std::stringstream linestream(line);
			std::string channel;
			std::getline(linestream, channel, ' ');
			std::string name;
			std::getline(linestream, name, ' ');
			double confidence;
			int top_x;
			int top_y;
			int bottom_x;
			int bottom_y;
			linestream >> confidence >> top_x >> top_y >> bottom_x >> bottom_y;
			//std::cout <<name<< " " << channel << " " << top_x << " " << top_y << " " << bottom_x << " " << bottom_y << " " << std::endl;
			Rect rect(top_x, top_y, bottom_x - top_x, bottom_y - top_y);
			if (channel == "alpha") {
				alpha_rects.push_back(rect);
				alpha_confidence.push_back(confidence);
			}

			else if (channel == "beta") {
				beta_rects.push_back(rect);
				beta_confidence.push_back(confidence);
			}

			else if (channel == "gamma") {
				gamma_rects.push_back(rect);
				gamma_confidence.push_back(confidence);
			}

		}
		out_rect_list.emplace_back(alpha_rects);
		out_rect_list.emplace_back(beta_rects);
		out_rect_list.emplace_back(gamma_rects);

		out_confidence_list.emplace_back(alpha_confidence);
		out_confidence_list.emplace_back(beta_confidence);
		out_confidence_list.emplace_back(gamma_confidence);

	}
	
	ParseAlphaBetaGammaSAOG(out_rect_list[0], out_rect_list[1], out_rect_list[2], aog, write_file,
		out_confidence_list[0], out_confidence_list[1], out_confidence_list[2]);
}

#include <boost/filesystem.hpp>
void LearnAndParseVideoImagesFromFolder(std::string input_file, std::string target_folder, std::string write_folder) {
	AOG<std::string, std::unordered_map<std::string, double>> aog = LearnAlphaBetaGammaAOGFromFile(input_file);

	std::vector<std::string> parseList;	
	std::vector<std::string> writeList;
	using namespace boost::filesystem;
	for (directory_iterator itr(target_folder); itr != directory_iterator(); ++itr)
	{
		stringstream ss;
		ss << itr->path().filename();
		cout << itr->path().filename() << ' '; // display filename only
		if (is_regular_file(itr->status())) cout << " [" << file_size(itr->path()) << ']';
		cout << '\n';

		std::string filename = ss.str();
		filename = filename.substr(1, filename.size() - 2);
		parseList.push_back(target_folder + "\\" + filename);
		writeList.push_back(write_folder + "\\" + filename);
	}

	for (size_t i = 0; i < parseList.size(); ++i) {
		std::vector<std::vector<Rect>> out_rect_list;
		std::vector<std::vector<double>> out_confidence_list;
		std::ifstream outfile(parseList[i]);
		if (outfile.is_open()) {
			std::vector<Rect> alpha_rects;
			std::vector<Rect> beta_rects;
			std::vector<Rect> gamma_rects;

			std::vector<double> alpha_confidence;
			std::vector<double> beta_confidence;
			std::vector<double> gamma_confidence;

			std::string line;
			while (std::getline(outfile, line)) {
				//random delete lines for debug
				//if (rand() % 100 < 60)
				//	continue;
				std::stringstream linestream(line);
				std::string channel;
				std::getline(linestream, channel, ' ');
				std::string name;
				std::getline(linestream, name, ' ');
				double confidence;
				int top_x;
				int top_y;
				int bottom_x;
				int bottom_y;
				linestream >> confidence >> top_x >> top_y >> bottom_x >> bottom_y;
				//std::cout <<name<< " " << channel << " " << top_x << " " << top_y << " " << bottom_x << " " << bottom_y << " " << std::endl;
				Rect rect(top_x, top_y, bottom_x - top_x, bottom_y - top_y);
				if (channel == "alpha") {
					alpha_rects.push_back(rect);
					alpha_confidence.push_back(confidence);
				}

				else if (channel == "beta") {
					beta_rects.push_back(rect);
					beta_confidence.push_back(confidence);
				}

				else if (channel == "gamma") {
					gamma_rects.push_back(rect);
					gamma_confidence.push_back(confidence);
				}

			}
			out_rect_list.emplace_back(alpha_rects);
			out_rect_list.emplace_back(beta_rects);
			out_rect_list.emplace_back(gamma_rects);

			out_confidence_list.emplace_back(alpha_confidence);
			out_confidence_list.emplace_back(beta_confidence);
			out_confidence_list.emplace_back(gamma_confidence);

		}
		
		std::cout << "writefile: " << writeList[i] << std::endl;
		std::cout << "out_rect_list len " << out_rect_list.size() << " out_confidence_list " << out_confidence_list.size() << std::endl;
		ParseAlphaBetaGammaSAOG(out_rect_list[0], out_rect_list[1], out_rect_list[2], aog, writeList[i],
			out_confidence_list[0], out_confidence_list[1], out_confidence_list[2]);
	}
}
#endif // !ALPHA_BETA_GAMMA_SAOG_H