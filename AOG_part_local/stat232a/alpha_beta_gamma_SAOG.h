#ifndef ALPHA_BETA_GAMMA_SAOG_H
#define	ALPHA_BETA_GAMMA_SAOG_H
#endif // !ALPHA_BETA_GAMMA_SAOG_H

#include "AOG.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/Imgproc.hpp>

using namespace AOG_LIB;

//template <class StateType, class AttributeType>
//class AlphaBetaGammaSAOG : public AOG<StateType, AttributeType> {
//public:
//	//default constructor
//	AlphaBetaGammaSAOG() {}
//
//	//construct from a set of rules. Currently it assigns the weight of 1.0 to all edges
//	AlphaBetaGammaSAOG(const std::vector<Symbolic_Rule<StateType>>& rules) {
//		this->has_root_ = false;
//		for (auto rule : rules)
//			this->AddRule(rule);
//	}
//	
//	//
//};

AOG<std::string, std::vector<double>> AlphaBetaGammaSAOG(std::string alpha_name, double alpha_weight,
	std::vector<std::string> beta_names, std::vector<double> beta_weights,
	std::string gamma_name, double gamma_weights){
	
	//initialize rules
	std::vector<Symbolic_Rule<std::string>> rules;

	Symbolic_State<std::string> gamma(gamma_name, false);
	Symbolic_State<std::string> alpha(alpha_name, false);
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
	return aog;
}