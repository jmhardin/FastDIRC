#include "../include/dirc_progressive_separation.h"
#include "../include/dirc_spread_gaussian.h"
int max_sim_phots;
int step_sim_phots;
double x,y,phi,theta;
double BAR, tracking_unc, ckov_unc;


DircProgressiveSeparation::DircProgressiveSeparation(\
	DircOpticalSim* imodel,\
	int imax_phots,\
	int istep_phots,\
	double isigma, \
	double x_unc,\
	double y_unc,\
	double t_unc,\
	double im1 /*= .4937*/,\
	double im2 /*= .1396*/,\
	double ithresh /*= 20*/,\
)
{
	max_sim_phots = imax_phots;
	step_sim_phots = istep_phots;
  
	mass_1 = im1;
	mass_2 = im2;
	ll_threshold = ithresh;
	
	//be good - try not to abuse side effects here...
	dirc_model = imodel;
	
	std::vector<dirc_points> no_hits;
	
	spread_func = new DircSpreadGaussian(\
		isigma,\
		no_hits,\
		x_unc,\
		y_unc,\
		t_unc);
}
	
void DircProgressiveSeparation::set_max_step_phots(int im, int is)
{
	max_sim_phots = im;
	step_sim_phots = is;
}

void DircProgressiveSeparation::set_masses(double im1,double im2)
{
	mass_1 = im1;
	mass_2 = im2;
}
void set_threshold(double ithresh)
{
	ll_threshold = ithresh;
}

double DircProgressiveSeparation::ll_diff(\
	std::vector<dirc_point> &hit_points, \
	int num_support,\
	double beta_1,\
	double beta_2)
{
	double ll_1,ll_2;
	//NOTE: arg 2 is ignored for beta >= 0
	std::vector<dirc_point> support_1 = dirc_model->sim_rand_n_photons(\
		n_sim_phots,\
		0,\
		BAR,\
		x,\
		y,\
		theta,\
		phi,\
		tracking_unc,\
		ckov_unc,\
		beta_1);
	
	std::vector<dirc_point> support_2 = dirc_model->sim_rand_n_photons(\
		n_sim_phots,\
		0,\
		BAR,\
		x,\
		y,\
		theta,\
		phi,\
		tracking_unc,\
		ckov_unc,\
		beta_2);
	
	ll_1 = spread_func->get_log_likelihood_new_support(hit_points, support_1);
	ll_2 = spread_func->get_log_likelihood_new_support(hit_points, support_2);
	
	return ll_2 - ll_1;
}

double get_ll_progressive(\
	std::vector<dirc_point> &hit_points,\
	int iBAR,\
	double iE,\
	double ix,\
	double iy,\
	double itheta,\
	double iphi,\
	double itracking_unc,\
	double ickov_unc,\
)
{
	double beta_1 = dirc_model->get_beta(iE,mass_1);
	double beta_2 = dirc_model->get_beta(iE,mass_2);
	
	x = ix;
	y = iy;
	theta = itheta;
	phi = iphi;
	tracking_unc = itracking_unc;
	ckov_unc = ickov_unc;
	
	double cur_ll = 0;
	
	for (int i = 0; i < max_sim_phots; i++)
	{
		cur_ll += ll_diff(\
			std::vector<dirc_point> &hit_points, \
			int step_sim_phots,\
			double beta_1,\
			double beta_2);
		
		if (fabs(cur_ll) > ll_threshold)
		{
			//Assume we've IDed it and return;
			return cur_ll;
		}
	}
}
	
double get_ll_max(\
	std::vector<dirc_point> &hit_points,\
	int iBAR,\
	double iE,\
	double ix,\
	double iy,\
	double itheta,\
	double iphi,\
	double itracking_unc,\
	double ickov_unc,\
)
{
	double beta_1 = dirc_model->get_beta(iE,mass_1);
	double beta_2 = dirc_model->get_beta(iE,mass_2);
	
	x = ix;
	y = iy;
	theta = itheta;
	phi = iphi;
	tracking_unc = itracking_unc;
	ckov_unc = ickov_unc;
	
	return ll_diff(\
		std::vector<dirc_point> &hit_points, \
		int max_sim_phots,\
		double beta_1,\
		double beta_2);
		
}
