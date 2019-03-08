#ifndef multichannel_alpha_beta_gamma_saog_h
#define multichannel_alpha_beta_gamma_saog_h

#include "alpha_beta_gamma_saog.h"

void ParseAlphaBetaGammaSAOG_MultiBeta(const std::vector<Rect>& alpha_rects, const std::vector<Rect>& beta_rects, const std::vector<Rect>& gamma_rects,
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

	cv::resize(frame, frame, Size(frame.cols, frame.rows));
	imshow("frame", frame);
	waitKey(0);

	//imwrite("C:\\Users\\Yizhou Zhao\\Desktop\\pic\\test_3_independent_no_background.jpg", frame);

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

		for (int k = 0; k < beta_rects.size(); ++k) {
			if (RectInside(beta_rects[k], alpha_rects[j]) && !visited_beta_rects[k] && beta_scores[k] > beta_threshold) {
				reconstruced_beta_rects.push_back(beta_rects[k]);
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
	}

	Mat frame2 = Mat::zeros(3024, 4032, CV_8UC3);
	//Mat frame2 = Mat::zeros(900, 1200, CV_8UC3);
	frame2 = cv::Scalar(255, 255, 255);
	//CV_Assert(reconstruced_gamma_rects.size() == reconstruced_beta_rects.size());
	//CV_Assert(reconstruced_alpha_rects.size() == reconstruced_beta_rects.size());

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

	cv::resize(frame2, frame2, Size(frame2.cols, frame2.rows));
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
		file.close();
	}
}


void LearnAndParseAlphaBetaGammaAOGFromFile_MultiBeta(std::string input_file, std::string target_file, std::string write_file) {
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

	ParseAlphaBetaGammaSAOG_MultiBeta(out_rect_list[0], out_rect_list[1], out_rect_list[2], aog, write_file,
		out_confidence_list[0], out_confidence_list[1], out_confidence_list[2]);
}



#endif // !multichannel_alpha_beta_gamma_saog
