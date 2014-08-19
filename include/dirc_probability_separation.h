#include "dirc_probability_spread.h"
#include "dirc_point.h"
#include <vector>

#ifndef DIRC_PROBABILITY_SPERATION
#define DIRC_PROBABILITY_SPERATION

class DircProbabilitySeparation
{
protected:
	DircProbabilitySpread* neg_dens;
	DircProbabilitySpread* pos_dens;
public:
	DircProbabilitySeparation(DircProbabilitySpread* ineg_dens, DircProbabilitySpread* ipos_dens);
	
	//make virtual for polymorphism later
	virtual double spread_function(double neg, double pos);
	double get_log_likelihood_spread_diff(std::vector<dirc_point> inpoints);
};

#endif