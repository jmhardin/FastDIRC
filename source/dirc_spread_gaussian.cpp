#include "../include/dirc_point.h"
#include "../include/dirc_spread_gaussian.h"
#include <vector>
#include <math.h>

DircSpreadGaussian::DircSpreadGaussian(\
	double isigma, \
	std::vector<dirc_point> isupport,\
	bool itest_time_dir /*=true*/)\
	: DircSpreadRadius(isupport)
{
	sigma2 = isigma*isigma;
	rand_gen = new TRandom3(123);
}
void DircSpreadGaussian::support_spread(double spread_sig)
{
	unsigned int start_size = support_points.size();
	int mult_add = 1;
	for (unsigned int i = 0; i < start_size; i++)
	{
		for (int j = 0; j < mult_add; j++)
		{
			dirc_point new_point;
			new_point.x = support_points[i].x + rand_gen->Gaus(0,spread_sig);
			new_point.y = support_points[i].y + rand_gen->Gaus(0,spread_sig);
			new_point.t = support_points[i].t;
			support_points.push_back(new_point);
		}
	}
}
void DircSpreadGaussian::support_x_weight()
{
	double x;
	for (unsigned int i = 0; i < support_points.size(); i++)
	{
		x = support_points[i].x;
		support_points[i].weight = 1 + sqrt(fabs(x));
	}
}
double DircSpreadGaussian::radius_spread_function(double r)
{
	if (r < 20*sqrt(sigma2))
	{
		return exp(-r*r/sigma2);
// 		return std::max(0.0,(1-r*r/sigma2)); //Epanechnikov
// 		return std::max(0.0,(1-r*r/sigma2)*(1-r*r/sigma2)); //quartic
		 // 		return std::min(sqrt(sigma2),sqrt(sigma2)/r);
	}
	else
	{
		return 0;
	}
}
