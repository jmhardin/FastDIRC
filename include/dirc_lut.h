#include "dirc_point.h"
#include <utility>
#include <vector> 

#ifndef DIRC_LUT
#define DIRC_LUT
class DircLUT
{
private:

	DircLUTEnum* pt_to_ind;

	std::vector<std::vector<std::pair<double,double> > > lu_table;

	//Needed?  Private?
	std::vector<std::pair<double,double> > get_base_phi_theta(dirc_point pt);
public:
	DircLUT(DircLUTEnum* ipt_to_ind);
	
	void add_table_pt(dirc_point pt, double phi, double theta);

	void get_ckov_theta_all(std::vector<double> &rval, std::vector<dirc_point> pt,double inc_phi, double inc_theta);
	//Needed?  Private?
	void get_base_phi_theta_all(std::vector<std::pair<double,double> > &rval, std::vector<dirc_point> pts);

};
#endif
