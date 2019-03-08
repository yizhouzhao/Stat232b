#include "pch.h"
#include "video_to_image.h"

#include <boost/filesystem.hpp>
void Video2Image(std::string input_file, std::string output_folder, int skip) {
	//_mkdir(output_folder.c_str());
	std::ifstream ifile(input_file);
	if (!ifile) throw("error: cannot open video");
	ifile.close();

	VideoCapture cap;
	VideoWriter video;
	Mat frame;

	cap.open(input_file);
	std::string output_file = input_file;
	//output_file.replace(output_file.end() - 4, output_file.end(), "_0.avi");
	//video.open(output_file, VideoWriter::fourcc('M', 'J', 'P', 'G'), 28, Size(cap.get(CAP_PROP_FRAME_WIDTH), cap.get(CAP_PROP_FRAME_HEIGHT)));

	int frame_count = 0;
	int time_start = getTickCount();
	while (waitKey(1) < 0) {
		cap >> frame;

		if (frame.empty()) {
			std::cout << "Done processing !!! In total: " << frame_count << " time elapse: "
				<< (getTickCount() - time_start) / getTickFrequency() << std::endl;
			std::cout << "Output file is stored as " << output_file << std::endl;
			waitKey(3000);
			break;
		}
		
		if (frame_count % skip == 0) {
			std::cout << "frame: " << frame_count;
			std::string temp_image_path = "temp.jpg";
			imwrite(temp_image_path, frame);
			std::string write_path = output_folder + "\\" + std::to_string(frame_count / skip) + ".txt";
			std::cout << "write_path: " << write_path << std::endl;
			std::ofstream ofile(write_path);
			ofile.close();
			predictImageYolov3(temp_image_path, "person", write_path, "gamma");
			//predictImageHaarCascade(temp_image_path, "head and shouler", write_path, "gamma");
			predictImageHaarCascade(temp_image_path, "face", write_path, "alpha");
			predictImageHaarCascade(temp_image_path, "nose", write_path, "beta");
		}

		frame_count++;
	}

	cap.release();
}

void VideoAddBox(std::string input_file, std::string box_folder, int skip) {
	//_mkdir(output_folder.c_str());
	std::ifstream ifile(input_file);
	if (!ifile) throw("error: cannot open video");

	VideoCapture cap;
	VideoWriter video;
	Mat frame;

	cap.open(input_file);
	std::string output_file = input_file;
	output_file.replace(output_file.end() - 4, output_file.end(), "_parsing_output.avi");
	video.open(output_file, VideoWriter::fourcc('M', 'J', 'P', 'G'), 5, Size(cap.get(CAP_PROP_FRAME_WIDTH), cap.get(CAP_PROP_FRAME_HEIGHT)));

	std::vector<std::string> boxList;
	using namespace boost::filesystem;
	for (directory_iterator itr(box_folder); itr != directory_iterator(); ++itr)
	{
		stringstream ss;
		ss << itr->path().filename();
		cout << itr->path().filename() << ' '; // display filename only
		if (is_regular_file(itr->status())) cout << " [" << file_size(itr->path()) << ']';
		cout << '\n';

		std::string filename = ss.str();
		filename = filename.substr(1, filename.size() - 2);
		boxList.push_back(box_folder + "\\" + filename);
	}


	int frame_count = 0;
	int time_start = getTickCount();
	while (waitKey(1) < 0) {
		cap >> frame;

		if (frame.empty()) {
			std::cout << "Done processing !!! In total: " << frame_count << " time elapse: "
				<< (getTickCount() - time_start) / getTickFrequency() << std::endl;
			std::cout << "Output file is stored as " << output_file << std::endl;
			waitKey(3000);
			break;
		}

		if (frame_count % skip == 0) {
			int idx = frame_count / skip;
			std::cout << "file here " << idx << " " << boxList.size() << std::endl;
			std::ifstream box_file(boxList[idx]);
			std::string line;
			while (std::getline(box_file, line)) {
				//random delete lines for debug
				//if (rand() % 100 < 60)
				//	continue;
				if (line.size() == 0) break;
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
				std::cout <<name<< " " << channel << " " << top_x << " " << top_y << " " << bottom_x << " " << bottom_y << " " << std::endl;
				Rect rect(top_x, top_y, bottom_x - top_x, bottom_y - top_y);
				if (channel == "alpha") {
					Scalar color(20, 180, 100);
					rectangle(frame, rect, color, 4, 8, 0);
				}
				else if (channel == "beta") {
					Scalar color(100, 20, 180);
					rectangle(frame, rect, color, 4, 8, 0);
				}
				else if (channel == "gamma") {
					Scalar color(180, 100, 20);
					rectangle(frame, rect, color, 4, 8, 0);
				}
			}

			Mat detectedFrame;
			frame.convertTo(detectedFrame, CV_8U);

			video.write(detectedFrame);
		}
		frame_count++;
	}

	video.release();
	std::cout << "Done processing F !!!\n";
}