#include "dirc_point.h"
#include <utility>
#include <vector> 

#ifndef DIRC_LUT
#define DIRC_LUT
struct lut_entry
{
	double phi;
	double theta;
	double time;
};

class DircLUT
{
	private:

		DircLUTEnum* pt_to_ind;

		std::vector<std::vector<lut_entry> > lu_table;

		//Needed?  Private?
		std::vector<lut_entry> get_base_phi_theta(dirc_point pt);
	public:
		DircLUT(DircLUTEnum* ipt_to_ind);

		void add_table_pt(dirc_point pt, double phi, double theta);

		void get_ckov_theta_all(std::vector<double> &rval, std::vector<double> &ret_dt, std::vector<dirc_point> pt,double inc_phi, double inc_theta, double inc_y);
		//Needed?  Private?
		void get_base_phi_theta_all(std::vector<lut_entry> &rval, std::vector<dirc_point> pts);
	/*
		void get_ckov_theta_single_oval_cut(std::vector<double> &rval, \
				std::vector<double> &ret_dt, \
				std::vector<dirc_point> pts, \
				double inc_phi, \
				double inc_theta, \
				double inc_y, \
				double center_ang, \
				double center_ang_spread,\
				double center_time_spread);
	*/
		void get_chromcorr_m_b_single_oval_cut(
				std::vector<dirc_point> pts, \
				double inc_phi, \
				double inc_theta, \
				double inc_y, \
				double center_ang, \
				double center_ang_spread_sq,\
				double time_spread_sq,\
				double expected_angle,\
				double &m_direct,\
				double &b_direct,\
				double &m_indirect,\
				double &b_indirect);
		void get_ckov_theta_single_oval_cut(
				std::vector<double> &rval, \
				std::vector<double> &ret_dt_dl, \
				std::vector<dirc_point> pts, \
				double inc_phi, \
				double inc_theta, \
				double inc_y, \
				double center_ang, \
				double center_ang_spread_sq,\
				double time_spread_sq,\
				double m_direct=0,\
				double b_direct=0,\
				double m_indirect=0,\
				double b_indirect=0);



};
#endif
