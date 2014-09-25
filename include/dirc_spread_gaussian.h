#include "dirc_point.h"
#include "dirc_spread_radius.h"
#include <vector>
#include <TRandom3.h>
#ifndef DIRC_SPREAD_GAUSSIAN
#define DIRC_SPREAD_GAUSSIAN
// class DircSpreadGaussian : public DircSpreadRadius
// class DircSpreadGaussian : public DircProbabilitySpread
class DircSpreadGaussian
{
private:
	double lin_slope, r_trans, sigma2, max_val;
	TRandom3 *rand_gen;
	std::vector<dirc_point> support_points;
	bool test_time_dir;
	double get_weight(dirc_point inpoint);
public:
	//by reference - later for speed
	DircSpreadGaussian(\
		double isigma, \
		std::vector<dirc_point> isupport,\
		bool itest_time_dir = true);
	void support_spread(double spread_sig);
	void support_x_weight();

	
	//I am so sorry.  I did this disgusting thing for speed
	inline double radius_spread_function(double r2) __attribute__((always_inline))
	{
		if (r2 < 5*sigma2)
		{
			return exp(-r2/sigma2);
// 	  		return std::max(0.0,(1-r2/(sigma2))); //Epanechnikov
	  // 		return std::max(0.0,(1-r*r/sigma2)*(1-r*r/sigma2)); //quartic
			  // 		return std::min(sqrt(sigma2),sqrt(sigma2)/r);
		}
		else
		{
			 return 0;
		}
	};
	inline double support_spread_function(dirc_point support, dirc_point test)__attribute__((always_inline))
	{
		double dx2,dy2;
		dx2 = support.x - test.x;
		dx2 *= dx2;
		dy2 = support.y - test.y;
		dy2 *= dy2;
		return radius_spread_function(dx2+dy2);
	};

	double get_single_log_likelihood(dirc_point inpoint);
	double get_log_likelihood(std::vector<dirc_point> inpoints);
};
#endif
