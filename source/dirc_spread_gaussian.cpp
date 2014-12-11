#include "../include/dirc_point.h"
#include "../include/dirc_spread_gaussian.h"
#include <vector>
#include <math.h>

DircSpreadGaussian::DircSpreadGaussian(\
	double isigma, \
	std::vector<dirc_point> isupport,\
	double x_unc,\
	double y_unc,\
	double t_unc)
{
	sigma2 = isigma*isigma;
	sigma2inv = 1/sigma2;
	
	
	x_sig2inv = 1/(x_unc*x_unc);
	y_sig2inv = 1/(y_unc*y_unc);
	t_sig2inv = 1/(t_unc*t_unc);
	
	
	rand_gen = new TRandom3(123);
	support_points = isupport;
	
	printf("presps: %u\n",support_points.size());
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

double DircSpreadGaussian::get_single_log_likelihood(dirc_point inpoint)
{
	double tprob = 0;
	double log_mult = 4;
	double weight = 1;
// 	weight = get_weight(inpoint);
	for (unsigned int j = 0; j < support_points.size(); j++)
	{
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
		printf("LL loop %d\n",i);
		tprob = 0;
		for (unsigned int j = 0; j < support_points.size(); j++)
		{
			printf("LL inner loop %d\n",j);
			tprob += support_spread_function(support_points[j],inpoints[i]);
			eval_count++;
		}
		tprob /= support_points.size();
		
		rval += weight*log_mult*log(tprob+10e-3);
	}
	rval -= log(inpoints.size());
	
	return rval;
}
double DircSpreadGaussian::get_log_likelihood_new_support(std::vector<dirc_point> &inpoints, std::vector<dirc_point> &t_support)
{
  //Ideally, I find a way tp do this without replicating code....
  //maybe change above func to call this
	double tprob = 0;
	double rval = 0;
	int eval_count = 0;
	double log_mult = 4;
	double weight = 1;
	for (unsigned int i = 0; i < inpoints.size(); i++)
	{
		tprob = 0;
// 		weight = get_weight(inpoints[i]);
		for (unsigned int j = 0; j < t_support.size(); j++)
		{
			tprob += support_spread_function(t_support[j],inpoints[i]);
			eval_count++;
		}
		tprob /= t_support.size();
		
		rval += weight*log_mult*log(tprob+10e-3);
	}
	rval -= log(inpoints.size());
	
	return rval;
}
double DircSpreadGaussian::get_weight(dirc_point inpoint)
{
	return 1;
}
