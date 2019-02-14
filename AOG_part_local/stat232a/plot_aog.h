#ifndef PLOT_AOG_H
#define PLOT_AOG_H
#endif // !PLOT_AOG_H

#include <iomanip>
#include <sstream> 
#include <math.h> 

#include "AOG.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/Imgproc.hpp>

using namespace AOG_LIB;
using namespace cv;

/* draws dashed circle */
void dash_circle(Mat& img, Point center, int radius, const void* color, int fill)
{
	Size size = img.size();
	size_t step = img.step;
	int pix_size = (int)img.elemSize();
	uchar* ptr = img.ptr();
	int err = 0, dx = radius, dy = 0, plus = 1, minus = (radius << 1) - 1;
	int inside = center.x >= radius && center.x < size.width - radius &&
		center.y >= radius && center.y < size.height - radius;

#define ICV_PUT_POINT( ptr, x )     \
        memcpy( ptr + (x)*pix_size, color, pix_size );

	while (dx >= dy)
	{
		int mask;
		int y11 = center.y - dy, y12 = center.y + dy, y21 = center.y - dx, y22 = center.y + dx;
		int x11 = center.x - dx, x12 = center.x + dx, x21 = center.x - dy, x22 = center.x + dy;

		if (inside)
		{
			uchar *tptr0 = ptr + y11 * step;
			uchar *tptr1 = ptr + y12 * step;

			if (!fill && dy % 2 == 0)
			{
				ICV_PUT_POINT(tptr0, x11);
				ICV_PUT_POINT(tptr1, x11);
				ICV_PUT_POINT(tptr0, x12);
				ICV_PUT_POINT(tptr1, x12);
			}
			else
			{
				//ICV_HLINE(tptr0, x11, x12, color, pix_size);
				//ICV_HLINE(tptr1, x11, x12, color, pix_size);
			}

			tptr0 = ptr + y21 * step;
			tptr1 = ptr + y22 * step;

			if (!fill && dy % 2 == 0)
			{
				ICV_PUT_POINT(tptr0, x21);
				ICV_PUT_POINT(tptr1, x21);
				ICV_PUT_POINT(tptr0, x22);
				ICV_PUT_POINT(tptr1, x22);
			}
			else
			{
				//ICV_HLINE(tptr0, x21, x22, color, pix_size);
				//ICV_HLINE(tptr1, x21, x22, color, pix_size);
			}
		}
		else if (x11 < size.width && x12 >= 0 && y21 < size.height && y22 >= 0)
		{
			if (fill)
			{
				x11 = std::max(x11, 0);
				x12 = MIN(x12, size.width - 1);
			}

			if ((unsigned)y11 < (unsigned)size.height)
			{
				uchar *tptr = ptr + y11 * step;

				if (!fill && dy % 2 == 0)
				{
					if (x11 >= 0)
						ICV_PUT_POINT(tptr, x11);
					if (x12 < size.width)
						ICV_PUT_POINT(tptr, x12);
				}
				//else
					//ICV_HLINE(tptr, x11, x12, color, pix_size);
			}

			if ((unsigned)y12 < (unsigned)size.height)
			{
				uchar *tptr = ptr + y12 * step;

				if (!fill && dy % 2 == 0)
				{
					if (x11 >= 0)
						ICV_PUT_POINT(tptr, x11);
					if (x12 < size.width)
						ICV_PUT_POINT(tptr, x12);
				}
				//else
					//ICV_HLINE(tptr, x11, x12, color, pix_size);
			}

			if (x21 < size.width && x22 >= 0)
			{
				if (fill)
				{
					x21 = std::max(x21, 0);
					x22 = MIN(x22, size.width - 1);
				}

				if ((unsigned)y21 < (unsigned)size.height)
				{
					uchar *tptr = ptr + y21 * step;

					if (!fill && dy % 2 == 0)
					{
						if (x21 >= 0)
							ICV_PUT_POINT(tptr, x21);
						if (x22 < size.width)
							ICV_PUT_POINT(tptr, x22);
					}
					//else
						//ICV_HLINE(tptr, x21, x22, color, pix_size);
				}

				if ((unsigned)y22 < (unsigned)size.height)
				{
					uchar *tptr = ptr + y22 * step;

					if (!fill && dy % 2 == 0)
					{
						if (x21 >= 0)
							ICV_PUT_POINT(tptr, x21);
						if (x22 < size.width)
							ICV_PUT_POINT(tptr, x22);
					}
					//else
						//ICV_HLINE(tptr, x21, x22, color, pix_size);
				}
			}
		}
		dy++;
		err += plus;
		plus += 2;

		mask = (err <= 0) - 1;

		err -= minus & mask;
		dx += mask;
		minus -= mask & 2;
	}

#undef  ICV_PUT_POINT
}

void dash_line(Mat& img, Point pt1, Point pt2, Scalar color, int gap = 5) {
	double dist = sqrt(pow(pt1.x - pt2.x, 2) + pow(pt1.y - pt2.y, 2));
	std::vector<Point> pts;

	for (int i = 0; i < dist; i += gap) {
		double r = i / dist;
		int x(pt1.x * (1 - r) + pt2.x * r + 0.5);
		int y(pt1.y * (1 - r) + pt2.y * r + 0.5);
		Point p(x, y);
		pts.push_back(p);
	}

	for (auto& p : pts) {
		cv::circle(img, p, 1, color, -1);
	}
}


template<class StateType, class AttributeType>
void PlotAOG(const AOG<StateType, AttributeType>& aog, VertexId root, HersheyFonts font = cv::FONT_HERSHEY_SIMPLEX) {
	//auto root_parent = aog.ParentsVertices(root);
	//if (!root_parent.empty() && !aog.GetVertexContent(root_parent[0])->IsAnd())
	//	root = root_parent[0];

	std::vector<std::vector<VertexId>> layer_list;
	std::vector<VertexId> layer0(1, root);
	layer_list.push_back(layer0);

	std::vector<VertexId> leaf_nodes;
	std::map<VertexId, int> num_of_child_leaf_nodes;
	int height(0);
	bool has_new = true;

	std::vector<VertexId> children = aog.ChildrenVertices(root);
	for (int j = 0; j < children.size(); ++j) {
		auto content = aog.GetStateByVertexId(children[j]).GetContent();
		std::cout << "root child content: " << content << " ";
	}


	while (has_new)
	{
		has_new = false;
		std::vector<VertexId> this_layer = layer_list.back();
		int size = this_layer.size();
		std::vector<VertexId> next_layer;
		for (int i = 0; i < size; i++) {
			VertexId parent_id = this_layer[i];
			num_of_child_leaf_nodes[parent_id] = 0;

			std::vector<VertexId> children = aog.ChildrenVertices(parent_id);

			if (children.size() == 0) {
				leaf_nodes.push_back(parent_id);
				num_of_child_leaf_nodes[parent_id] = 1;
			}

			else {
				has_new = true;
				for (int j = 0; j < children.size(); ++j) {
					auto content = aog.GetStateByVertexId(children[j]).GetContent();
					next_layer.push_back(children[j]);
					std::cout << "child content: " << content << " ";
				}
					
			}
		}
		layer_list.emplace_back(next_layer);
		height++;
		std::cout << height <<std::endl;
	}

	for (int i = layer_list.size() - 2; i >= 0; --i) {
		std::vector<VertexId> this_layer = layer_list[i];
		int size = this_layer.size();
		for (int j = 0; j < size; ++j) {
			VertexId parent_id = this_layer[j];
			std::vector<VertexId> children = aog.ChildrenVertices(parent_id);
			int children_size = children.size();
			//std::cout << "children size: " << children_size << std::endl;
			for (int k = 0; k < children_size; ++k) {
				num_of_child_leaf_nodes[parent_id] += num_of_child_leaf_nodes[children[k]];
			}
		}
	}

	for (int i = 0; i < layer_list.size(); i++) {
		for (int j = 0; j < layer_list[i].size(); j++) {
			std::cout << num_of_child_leaf_nodes[layer_list[i][j]] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "leaf total: " << leaf_nodes.size() << std::endl;
	std::cout << "height: " << height << std::endl;

	int height_per_layer = 100;
	int width_per_node = 60;

	Mat frame = Mat::zeros(height * height_per_layer, (leaf_nodes.size() + 1) * width_per_node, CV_8UC3);
	frame = cv::Scalar(255, 255, 255);

	Scalar color(0, 0, 0);
	for (int i = 0; i < layer_list.size(); i++) {
		int cum_nodes = 0;
		int cum_children = 0;
		for (int j = 0; j < layer_list[i].size(); j++) {
			VertexId cur_vertex_id = layer_list[i][j];
			std::vector<VertexId> parents_id = aog.ParentsVertices(cur_vertex_id);
			//if(parents_id.size() > 0 && !aog.GetVertexContent(parents_id[0])->IsAnd()) continue;

			int nodes  = num_of_child_leaf_nodes[cur_vertex_id];
			Point center((nodes / 2.0 + cum_nodes + 1) * width_per_node, i * height_per_layer + 50);
			auto content = aog.GetStateByVertexId(cur_vertex_id).GetContent();
			auto and_or = aog.GetVertexContent(cur_vertex_id)->IsAnd();
			std::cout << content << " id: " << cur_vertex_id << " and or: " << and_or << " " << center.x << " " << center.y << " ";
			int radius = width_per_node / 3;
			//rectangle(frame, Point(100, 100), Point(200, 200), color, 5, 8);
			//dash_circle(frame, center, radius, color, 2, cv::LINE_8);
			
			Point arrow_start(center.x, center.y + width_per_node / 2.5);
			if (and_or == 0) {
				dash_circle(frame, center, radius, &color, 0);
				std::unordered_map<VertexId, double> edge_weights = aog.GetOutEdgeWeights(cur_vertex_id, true);
				for (auto w : edge_weights)
				{
					VertexId cur_child = w.first;
					int child_nodes = num_of_child_leaf_nodes[cur_child];
					Point arrow_end((child_nodes / 2.0 + cum_children + 1) * width_per_node, 
						(i + 1) * height_per_layer + 50 - width_per_node / 2.5);
					//cv::arrowedLine(frame, arrow_start, arrow_end, color, 1, 8, 0, 0.04);
					dash_line(frame, arrow_start, arrow_end, color);

					Point arrow_middle((arrow_start.x + 2 * arrow_end.x) / 3, (arrow_start.y + 2 * arrow_end.y) / 3);
					double weight = std::floor((w.second * 100) + .5) / 100;
					std::stringstream ss;
					ss << std::fixed << std::setprecision(2) << weight;
					std::string mystring = ss.str();

					std::cout << "weight: " << std::to_string(weight) << std::endl;
					putText(frame, mystring, arrow_middle, font, 0.4, (255, 255, 255), 1, cv::LINE_8);

					cum_children += child_nodes;
				}

			}
				
			else{
				circle(frame, center, radius, color, 2, 8, 0);
				std::vector<VertexId> children = aog.ChildrenVertices(cur_vertex_id);

				Point arrow_start(center.x, center.y + width_per_node / 2.5);
				int children_size = children.size();
				for (int k = 0; k < children_size; ++k) {
					VertexId cur_child = children[k];
					int child_nodes = num_of_child_leaf_nodes[cur_child];
					Point arrow_end((child_nodes / 2.0 + cum_children + 1) * width_per_node,
						(i + 1) * height_per_layer + 50 - width_per_node / 2.5);
					cv::line(frame, arrow_start, arrow_end, color, 1, 8, 0);

					cum_children += child_nodes;
				}
			}
				
			Point text_center(center.x - width_per_node / 8, center.y);
			putText(frame, content, text_center, font, 0.45, Scalar(80, 50, 255), 1, cv::LINE_8);

			cum_nodes += nodes;
			//std::cout << num_of_child_leaf_nodes[layer_list[i][j]] << " ";
		}
		std::cout << std::endl;
	}


	Point p1(0, 0);                            // start & end points 
	Point p2(300, 300);
	LineIterator it(frame, p1, p2, 8);            // get a line iterator
	for (int i = 0; i < it.count; i++, it++)
		if (i % 5 != 0) { (*it)[0] = 255; }         // every 5'th pixel gets dropped, blue stipple line


	cv::imshow("image-", frame);
	cv::imwrite("C:\\Users\\Yizhou Zhao\\Desktop\\demo.jpg", frame);
	cv::waitKey(0);


}



//template<class StateType, class AttributeType>
//int GetAOGHeight(const AOG<StateType, AttributeType>& aog, VertexId root) {
//	std::list<VectexId> queue;
//	queue.push_back(root);
//	
//	int height = 0;
//	while (!queue.empty()) {
//		int size = queue.size();
//		while (size--) {
//			auto front = queue.front();
//			queue.pop_front();
//
//		}
//	}
//}

