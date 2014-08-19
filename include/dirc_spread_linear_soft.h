#include "dirc_point.h"
#include "dirc_spread_radius.h"
#include <vector>
#ifndef DIRC_SPREAD_LINEAR_SOFT
#define DIRC_SPREAD_LINEAR_SOFT
class DircSpreadLinearSoft : public DircSpreadRadius
{
private:
	double lin_slope, r_trans, sigma2, max_val;
public:
	//by reference - later for speed
	DircSpreadLinearSoft(\
		double ilin_slope, \
		double ir_trans, \
		double isigma, \
		std::vector<dirc_point> isupport,\
		bool itest_time_dir = true);
	
	double radius_spread_function(double r);
};
#endif
