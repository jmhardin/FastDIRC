#include "dirc_point.h"
#include "dirc_spread_radius.h"
#include <vector>
#include <TRandom3.h>
#ifndef DIRC_SPREAD_GAUSSIAN
#define DIRC_SPREAD_GAUSSIAN
// class DircSpreadGaussian : public DircSpreadRadius
class DircSpreadGaussian : public DircProbabilitySpread
{
private:
	double lin_slope, r_trans, sigma2, max_val;
	TRandom3 *rand_gen;
public:
	//by reference - later for speed
	DircSpreadGaussian(\
		double isigma, \
		std::vector<dirc_point> isupport,\
		bool itest_time_dir = true);
	void support_spread(double spread_sig);
	void support_x_weight();
	double support_spread_function(dirc_point support, dirc_point test);
	double radius_spread_function(double r);
};
#endif
