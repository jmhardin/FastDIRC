#include "dirc_optical_sim.h"
#include "dirc_point.h"
#include "dirc_spread_gaussian.h"

#ifndef DIRC_PROGRESSIVE_SEPARATION
#define DIRC_PROGRESSIVE_SEPARATION

class DircProgressiveSeparation
{
protected:
	int max_sim_phots;
	int step_sim_phots;
	double ll_threshold;
	double E,x,y,phi,theta;
	double mass_1,mass_2;
	double BAR, tracking_unc, ckov_unc;
	
	DircOpticalSim* dirc_model;
	DircSpreadGaussian* spread_func;
	
	double ll_diff(\
		std::vector<dirc_point> &hit_points, \
		int num_support,\
		double beta_1,\
		double beta_2);
	double ll_pts(\
		std::vector<dirc_point> &hit_points, \
		int num_support,\
		double beta);
	
	void get_inf_and_ll(\
		int &num_neg_inf,\
		double &ll_sum,\
		std::vector<double> prob_vals);
public:
	DircProgressiveSeparation(\
		DircOpticalSim* imodel,\
		int imax_phots,\
		int istep_phots,\
		double isigma, \
		double x_unc,\
		double y_unc,\
		double t_unc,\
		double im1 /*= .4937*/,\
		double im2 /*= .1396*/,\
		double ithresh /*= 20*/);
	
	void set_masses(double im1,double im2);
	void set_threshold(double ithresh);
	void set_max_step_phots(int im, int is);
	
	double get_ll_progressive(\
		std::vector<dirc_point> &hit_points,\
		int iBAR,\
		double iE,\
		double ix,\
		double iy,\
		double itheta,\
		double iphi,\
		double itracking_unc,\
		double ickov_unc);
	double get_ll_max(\
		std::vector<dirc_point> &hit_points,\
		int iBAR,\
		double iE,\
		double ix,\
		double iy,\
		double itheta,\
		double iphi,\
		double itracking_unc,\
		double ickov_unc);
};
#endif