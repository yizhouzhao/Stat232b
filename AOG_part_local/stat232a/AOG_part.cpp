// stat232a.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

/*
#include "pch.h"
#include <iostream>
#include <memory>
#include <random>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include <stack>
#include <utility>
#include <algorithm>
//#include<boost/functional/hash.hpp>

#include "AOG_Edge.h"
#include "Symbolic_State.h"
#include "Symbolic_Rule.h"
#include "Core/Graph.hpp"
#include "AOG_Vertex.h"
#include "AOG.h"
*/

#include "pch.h"
#include <string>
#include <iostream>
#include "AOG.h"
#include "plot_aog.h"
#include "alpha_beta_gamma_saog.h"

using namespace std;
using namespace AOG_LIB;

template <class T>
using SequenceType = std::vector<AOG_LIB::Symbolic_State<T> >;

double SimplePotentialFunc(AOG_Vertex<string, string>& vertex_self,
	const vector<shared_ptr<AOG_Vertex<string, string> > >& neighbors)
{
	double potential = 0;
	for (auto neighbor : neighbors)
		if (neighbor->GetAttribute() == vertex_self.GetAttribute())
			potential += neighbor->GetPotential();
	return potential;
}

//output the sample result to a file
void SampleToFile(string filename, const AOG<string, string>& t_aog, VertexId sample_root, int num_samples = 100)
{
	//open file
	ofstream outFile;
	outFile.open(filename);
	if (outFile.is_open())
	{
		while (num_samples)
		{
			vector<VertexId > res;
			double sample_prob = 1;
			t_aog.Sample(sample_root, res, sample_prob);
			num_samples--;
			for (auto vid : res)
			{
				outFile << t_aog.GetStateByVertexId(vid).GetContent() << " ";

			}
			outFile << "\n";
		}
	}
	else
	{
		std::cerr << "Error opening file!\n";
		throw exception();
	}
	outFile.close();
}


int main()
{ 
	
	// an example
	/*
	//Define an AOG
	vector<Symbolic_Rule<string> > rules;
	Symbolic_State<string> root("O", false);
	Symbolic_State<string> A("A", false);
	Symbolic_State<string> B("B", false);
	Symbolic_State<string> C("C", false);
	vector<Symbolic_State<string> > level1A = { A };
	vector<Symbolic_State<string> > level1B = { B };
	vector<Symbolic_State<string> > level1C = { C };
	Symbolic_Rule<string> level1Arule(root, level1A);
	Symbolic_Rule<string> level1Brule(root, level1B);
	Symbolic_Rule<string> level1Crule(root, level1C);
	rules.push_back(level1Arule);
	rules.push_back(level1Brule);
	rules.push_back(level1Crule);

	// A
	Symbolic_State<string> a1("a1", false);
	Symbolic_State<string> a2("a2", false);
	Symbolic_State<string> a3("a3", false);
	vector<Symbolic_State<string> > a1a2a3 = { a1, a2, a3 };
	Symbolic_Rule<string> A_a1a2a3(A, a1a2a3);
	rules.push_back(A_a1a2a3);

	Symbolic_State<string> a1A1("a1A1", true);
	Symbolic_State<string> a1A2("a1A2", true);
	Symbolic_State<string> a1B1("a1B1", true);
	Symbolic_State<string> a1B2("a1B2", true);
	vector<Symbolic_State<string> > a1A1a1A2 = { a1A1, a1A2 };
	vector<Symbolic_State<string> > a1B1a1B2 = { a1B1, a1B2 };
	Symbolic_Rule<string> a1_a1A1a1A2(a1, a1A1a1A2);
	Symbolic_Rule<string> a1_a1B1a1B2(a1, a1B1a1B2);
	rules.push_back(a1_a1A1a1A2);
	rules.push_back(a1_a1B1a1B2);

	Symbolic_State<string> a21("a21", true);
	Symbolic_State<string> a22("a22", true);
	vector<Symbolic_State<string> > va21 = { a21 };
	vector<Symbolic_State<string> > va22 = { a22 };
	Symbolic_Rule<string> a2_a21(a2, va21);
	Symbolic_Rule<string> a2_a22(a2, va22);

	rules.push_back(a2_a21);
	rules.push_back(a2_a22);


	Symbolic_State<string> a31("a31", true);
	Symbolic_State<string> a32("a32", true);
	vector<Symbolic_State<string> > va31 = { a31 };
	vector<Symbolic_State<string> > va32 = { a32 };
	Symbolic_Rule<string> a3_a31(a3, va31);
	Symbolic_Rule<string> a3_a32(a3, va32);
	rules.push_back(a3_a31);
	rules.push_back(a3_a32);


	// B
	Symbolic_State<string> b1("b1", false);
	Symbolic_State<string> b2("b2", false);
	vector<Symbolic_State<string> > b1b2 = { b1, b2 };
	Symbolic_Rule<string> B_b1b2(B, b1b2);
	rules.push_back(B_b1b2);

	Symbolic_State<string> B1A1("b1A1", true);
	Symbolic_State<string> B1A2("b1A2", true);
	Symbolic_State<string> B1A3("b1A3", true);
	Symbolic_State<string> B1B1("b1B1", true);
	Symbolic_State<string> B1B2("b1B2", true);

	vector<Symbolic_State<string> > B1A1B1A2B1A3 = { B1A1, B1A2, B1A3 };
	vector<Symbolic_State<string> > B1B1B1B2 = { B1B1, B1B2 };
	Symbolic_Rule<string> b1_B1A1B1A2B1A3(b1, B1A1B1A2B1A3);
	Symbolic_Rule<string> b1_B1B1B1B2(b1, B1B1B1B2);
	rules.push_back(b1_B1A1B1A2B1A3);
	rules.push_back(b1_B1B1B1B2);

	Symbolic_State<string> B2A1("b2A1", true);
	Symbolic_State<string> B2A2("b2A2", true);
	Symbolic_State<string> B2B1("b2B1", true);
	Symbolic_State<string> B2B2("b2B2", true);
	vector<Symbolic_State<string> > B2A1B2A2 = { B2A1, B2A2 };
	vector<Symbolic_State<string> > B2B1B2B2 = { B2B1, B2B2 };
	Symbolic_Rule<string> b2_B2A1B2A2(b2, B2A1B2A2);
	Symbolic_Rule<string> b2_B2B1B2B2(b2, B2B1B2B2);
	rules.push_back(b2_B2A1B2A2);
	rules.push_back(b2_B2B1B2B2);


	// C
	Symbolic_State<string> c1("c1", false);
	Symbolic_State<string> c2("c2", false);
	vector<Symbolic_State<string> > c1c2 = { c1, c2 };
	Symbolic_Rule<string> C_c1c2(C, c1c2);
	rules.push_back(C_c1c2);

	Symbolic_State<string> C1A1("c1A1", true);
	Symbolic_State<string> C1A2("c1A2", true);
	Symbolic_State<string> C1A3("c1A3", true);
	Symbolic_State<string> C1B1("c1B1", true);
	Symbolic_State<string> C1B2("c1B2", true);

	vector<Symbolic_State<string> > C1A1C1A2C1A3 = { C1A1, C1A2, C1A3 };
	vector<Symbolic_State<string> > C1B1C1B2 = { C1B1, C1B2 };
	Symbolic_Rule<string> c1_C1A1C1A2C1A3(c1, C1A1C1A2C1A3);
	Symbolic_Rule<string> c1_C1B1C1B2(c1, C1B1C1B2);
	rules.push_back(c1_C1A1C1A2C1A3);
	rules.push_back(c1_C1B1C1B2);

	Symbolic_State<string> C2A1("c2A1", true);
	Symbolic_State<string> C2A2("c2A2", true);
	Symbolic_State<string> C2B1("c2B1", true);
	Symbolic_State<string> C2B2("c2B2", true);
	vector<Symbolic_State<string> > C2A1C2A2 = { C2A1, C2A2 };
	vector<Symbolic_State<string> > C2B1C2B2 = { C2B1, C2B2 };
	Symbolic_Rule<string> c2_C2A1C2A2(c2, C2A1C2A2);
	Symbolic_Rule<string> c2_C2B1C2B2(c2, C2B1C2B2);
	rules.push_back(c2_C2A1C2A2);
	rules.push_back(c2_C2B1C2B2);

	AOG<string, string> t_aog(rules);
	t_aog.SetRoot(root);
	VertexId root_id = t_aog.GetVertexIdByState(root);
	t_aog.SetVertexAttribute(root_id, "good");
	t_aog.SetVertexPotentialFunc(root_id, SimplePotentialFunc);
	t_aog.AddVertexNeighbor(root_id, t_aog.GetVertexIdByState(A));
	t_aog.AddVertexNeighbor(root_id, t_aog.GetVertexIdByState(B));
	t_aog.AddVertexNeighbor(root_id, t_aog.GetVertexIdByState(C));

	t_aog.SetVertexPotential(t_aog.GetVertexIdByState(A), 1);
	t_aog.SetVertexAttribute(t_aog.GetVertexIdByState(A), "good");
	t_aog.SetVertexPotential(t_aog.GetVertexIdByState(B), 2);
	t_aog.SetVertexAttribute(t_aog.GetVertexIdByState(B), "bad");
	t_aog.SetVertexPotential(t_aog.GetVertexIdByState(C), 3);
	t_aog.SetVertexAttribute(t_aog.GetVertexIdByState(C), "good");
	t_aog.UpdateVertexPotential(root_id);
	cout << "root node potential is: " << t_aog.GetVertexContent(root_id)->GetPotential() << "\n";

	// A's weight
	unordered_map<VertexId, double> weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(root), true);
	for (auto it = weights.begin(); it != weights.end(); it++)
	{
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == A)
			it->second = 0.5;
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == B)
			it->second = 0.3;
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == C)
			it->second = 0.2;
	}
	t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(root), weights);


	weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(a1), true);
	for (auto it = weights.begin(); it != weights.end(); it++)
	{
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == a1A1)
			it->second = 0.7;
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == a1B1)
			it->second = 0.3;
	}
	t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(a1), weights);


	weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(a2), true);
	for (auto it = weights.begin(); it != weights.end(); it++)
	{
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == a21)
			it->second = 0.46;
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == a22)
			it->second = 0.54;
	}
	t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(a2), weights);


	weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(a3), true);
	for (auto it = weights.begin(); it != weights.end(); it++)
	{
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == a31)
			it->second = 0.76;
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == a32)
			it->second = 0.24;
	}
	t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(a3), weights);

	// B's weight
	weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(b1), true);
	for (auto it = weights.begin(); it != weights.end(); it++)
	{
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == B1A1)
			it->second = 0.21;
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == B1B1)
			it->second = 0.79;
	}
	t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(b1), weights);

	weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(b2), true);
	for (auto it = weights.begin(); it != weights.end(); it++)
	{
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == B2A1)
			it->second = 0.64;
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == B2B1)
			it->second = 0.36;
	}
	t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(b2), weights);


	// C's weight
	weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(c1), true);
	for (auto it = weights.begin(); it != weights.end(); it++)
	{
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == C1A1)
			it->second = 0.71;
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == C1B1)
			it->second = 0.29;
	}
	t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(c1), weights);

	weights = t_aog.GetOutEdgeWeights(t_aog.GetVertexIdByState(c2), true);
	for (auto it = weights.begin(); it != weights.end(); it++)
	{
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == C2A1)
			it->second = 0.34;
		if (t_aog.GetStateByVertexId(t_aog.ChildrenVertices(it->first)[0]) == C2B1)
			it->second = 0.66;
	}
	t_aog.SetOutEdgeWeights(t_aog.GetVertexIdByState(c2), weights);

	//t_aog.SaveGraph("./", "grammar_example.txt");
	//t_aog.Visualize("./", "./grammar_example_vis.txt");

	//SampleToFile("./grammar_example_samples.txt", t_aog, t_aog.GetVertexIdByState(root));
	*/
	

	////another simple example
	//vector<Symbolic_Rule<string> > rules;
	//Symbolic_State<string> root("O", false);
	//Symbolic_State<string> A("A", true);
	//Symbolic_State<string> B("B", true);
	//Symbolic_State<string> C("C", true);
	//vector<Symbolic_State<string> > level1A = { A };
	//vector<Symbolic_State<string> > level1B = { B };
	//vector<Symbolic_State<string> > level1C = { C };
	//Symbolic_Rule<string> level1Arule(root, level1A);
	//Symbolic_Rule<string> level1Brule(root, level1B);
	//Symbolic_Rule<string> level1Crule(root, level1C);
	//rules.push_back(level1Arule);
	//rules.push_back(level1Brule);
	//rules.push_back(level1Crule);

	//AOG<string, string> t_aog(rules);
	//t_aog.SetRoot(root);
	//
	//Mat frame = PlotAOG(t_aog);
	//imshow("image", frame);
	//waitKey(0);

	////AlphaBetaGammaSAOG<string, string> example_aog(rules);

	AOG<std::string, std::vector<double>> example_aog = AlphaBetaGammaSAOG("face", 0.45,
		{ "right eye", "left eye" }, 0.8, "body", 0.10);

	////example_aog.SaveGraph("./", "abc_grammar_example.txt");

	//Mat frame2 = PlotAlphaBetaGammaSAOG(example_aog);
	//imshow("image2", frame2);
	//imwrite("C:\\Users\\Yizhou Zhao\\Desktop\\pic\\saog.jpg", frame2);
	//waitKey(0);

	std::vector<std::vector<Rect>> rect_list = ReadRectFromFile("classroom_31.txt");
	//std::cout << rect_list[0][0].x;

	AOG<std::string, std::vector<double>> learn_aog = 
		LearnAlphaBetaGammaSAOG(rect_list[0], rect_list[1], rect_list[2]);

	Mat frame3 = PlotAlphaBetaGammaSAOG(learn_aog);
	//imshow("image3", frame3);
	//waitKey(0);

	std::cout << NormalDensity1D(0, 0, 1) << std::endl;

	std::cout << "Hello World!\n";

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
