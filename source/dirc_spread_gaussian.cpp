#include "../include/dirc_point.h"
#include "../include/dirc_spread_gaussian.h"
#include <vector>
#include <math.h>

// DircSpreadGaussian::DircSpreadGaussian(\
// 	double isigma, \
// 	std::vector<dirc_point> isupport,\
// 	bool itest_time_dir /*=true*/)\
// 	: DircProbabilitySpread(isupport)
// {
// 	sigma2 = isigma*isigma;
// 	rand_gen = new TRandom3(123);
// }
// void DircSpreadGaussian::support_spread(double spread_sig)
// {
// 	unsigned int start_size = support_points.size();
// 	int mult_add = 1;
// 	for (unsigned int i = 0; i < start_size; i++)
// 	{
// 		for (int j = 0; j < mult_add; j++)
// 		{
// 			dirc_point new_point;
// 			new_point.x = support_points[i].x + rand_gen->Gaus(0,spread_sig);
// 			new_point.y = support_points[i].y + rand_gen->Gaus(0,spread_sig);
// 			new_point.t = support_points[i].t;
// 			support_points.push_back(new_point);
// 		}
// 	}
// }
// void DircSpreadGaussian::support_x_weight()
// {
// 	double x;
// 	for (unsigned int i = 0; i < support_points.size(); i++)
// 	{
// 		x = support_points[i].x;
// 		support_points[i].weight = 1 + sqrt(fabs(x));
// 	}
// }
// double DircSpreadGaussian::support_spread_function(dirc_point support, dirc_point test)
// {
// 	double dx2,dy2;
// 	dx2 = support.x - test.x;
// 	dx2 *= dx2;
// 	dy2 = support.y - test.y;
// 	dy2 *= dy2;
// 	return radius_spread_function(dx2+dy2);
// 	
// }
// double DircSpreadGaussian::radius_spread_function(double r2)
// {
// 	if (r2 < 5*sigma2)
// 	{
// 		return exp(-r2/sigma2);
// // 		return std::max(0.0,(1-r*r/sigma2)); //Epanechnikov
// // 		return std::max(0.0,(1-r*r/sigma2)*(1-r*r/sigma2)); //quartic
// 		 // 		return std::min(sqrt(sigma2),sqrt(sigma2)/r);
// 	}
// 	else
// 	{
// 		return 0;
// 	}
// }
DircSpreadGaussian::DircSpreadGaussian(\
	double isigma, \
	std::vector<dirc_point> isupport,\
	bool itest_time_dir /*=true*/)
{
	sigma2 = isigma*isigma;
	rand_gen = new TRandom3(123);
	support_points = isupport;
	test_time_dir = itest_time_dir;
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
// double DircSpreadGaussian::support_spread_function(dirc_point support, dirc_point test)
// {
// 	double dx2,dy2;
// 	dx2 = support.x - test.x;
// 	dx2 *= dx2;
// 	dy2 = support.y - test.y;
// 	dy2 *= dy2;
// 	return radius_spread_function(dx2+dy2);
// 	
// }
// double DircSpreadGaussian::radius_spread_function(double r2)
// {
// 	if (r2 < 5*sigma2)
// 	{
// 		return exp(-r2/sigma2);
// // 		return std::max(0.0,(1-r*r/sigma2)); //Epanechnikov
// // 		return std::max(0.0,(1-r*r/sigma2)*(1-r*r/sigma2)); //quartic
// 		 // 		return std::min(sqrt(sigma2),sqrt(sigma2)/r);
// 	}
// 	else
// 	{
// 		return 0;
// 	}
// }
double DircSpreadGaussian::get_single_log_likelihood(dirc_point inpoint)
{
	double tprob = 0;
	double log_mult = 4;
	double weight = 1;
// 	weight = get_weight(inpoint);
	for (unsigned int j = 0; j < support_points.size(); j++)
	{
		if (test_time_dir && support_points[j].t*inpoint.t < -.00001)
		{
			//Only look at if they went up or down for now
			continue;
		}
		tprob += support_spread_function(support_points[j],inpoint);
	}
	tprob /= support_points.size();
	tprob = weight*log_mult*log(tprob + 10e-3);
	
	return tprob;
}
double DircSpreadGaussian::get_log_likelihood(std::vector<dirc_point> inpoints)
{
	double tprob = 0;
	double rval = 0;
	int eval_count = 0;
	double log_mult = 4;
	double weight = 1;
	for (unsigned int i = 0; i < inpoints.size(); i++)
	{
		tprob = 0;
// 		weight = get_weight(inpoints[i]);
		for (unsigned int j = 0; j < support_points.size(); j++)
		{
			if (test_time_dir && support_points[j].t*inpoints[i].t < -.00001)
			{
				//Only look at if they went up or down for now
				continue;
			}
			tprob += support_spread_function(support_points[j],inpoints[i]);
			eval_count++;
		}
		tprob /= support_points.size();
		
		rval += weight*log_mult*log(tprob);
	}
// 	printf("eval_count: %d\n",eval_count);
	rval -= log(inpoints.size());
	
	return rval;
}
double DircSpreadGaussian::get_weight(dirc_point inpoint)
{
	return 1;
// 	return 1 + sqrt(sqrt(fabs(inpoint.x)));
}