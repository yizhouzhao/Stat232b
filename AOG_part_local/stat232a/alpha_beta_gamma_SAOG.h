#ifndef ALPHA_BETA_GAMMA_SAOG_H
#define	ALPHA_BETA_GAMMA_SAOG_H

#include "stdlib.h"
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

AOG<std::string, std::vector<double>> LearnAlphaBetaGammaSAOG(std::vector<Rect> alpha_rects, std::vector<Rect> beta_rects, std::vector<Rect> gamma_rects) {
	int gamma_overlap_alpha_count = 0;
	for (size_t i = 0; i < gamma_rects.size(); ++i) {
		for (size_t j = 0; j < alpha_rects.size(); ++j) {
			if (RectOverlap(gamma_rects[i], alpha_rects[j])) {
				gamma_overlap_alpha_count++;
				break;
			}
		}
	}

	//std::cout << "\n gamma_overlap_alpha_count " << gamma_overlap_alpha_count << std::endl;
	double gamma_to_alpha = gamma_overlap_alpha_count / (gamma_rects.size() + .0);

	//std::cout << "gamma_overlap_ratio " << gamma_to_alpha << std::endl;

	int beta_overlap_alpha_count = 0;
	for (size_t i = 0; i < beta_rects.size(); ++i) {
		for (size_t j = 0; j < alpha_rects.size(); ++j) {
			if (RectOverlap(beta_rects[i], alpha_rects[j])) {
				beta_overlap_alpha_count++;
				break;
			}
		}
	}

	double beta_to_alpha = beta_overlap_alpha_count / (beta_rects.size() + .0);

	AOG<std::string, std::vector<double>> aog = AlphaBetaGammaSAOG("alpha", 1, { "beta" }, beta_to_alpha, "gamma", gamma_to_alpha);
	return aog;
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



#endif // !ALPHA_BETA_GAMMA_SAOG_H