#include "../include/dirc_lut_enum.h"
#include "../include/dirc_lut.h"
#include <math.h>
#include <stdio.h>

DircLUT::DircLUT(DircLUTEnum* ipt_to_ind)
{
	pt_to_ind = ipt_to_ind;

	lu_table.resize(pt_to_ind->get_max_enum());
	for (unsigned int i = 0; i < lu_table.size(); i++)
	{
		lu_table[i].clear();
	}
}
void DircLUT::add_table_pt(dirc_point pt, double phi, double theta)
{
	std::pair<double,double> add_pair;
	add_pair.first = phi;
	add_pair.second = theta;
	int ind = pt_to_ind->return_enum(pt);

	if (ind >= 0)
	{
		lu_table[ind].push_back(add_pair);
	}
}
std::vector<std::pair<double,double> > DircLUT::get_base_phi_theta(dirc_point pt)
{
	int ind = pt_to_ind->return_enum(pt);

	if (ind >= 0)
	{
		return lu_table[ind];
	}
	else
	{
		// :(  Something is rotten in the state of Denmark
		std::vector<std::pair<double,double> > empty;
		empty.clear();
		return empty;
	}
}
void DircLUT::get_base_phi_theta_all(std::vector<std::pair<double,double> > &rval, std::vector<dirc_point> pts)
{
//	std::vector<std::pair<double,double> > rval;
	rval.clear();
	
	for (unsigned int i = 0; i < pts.size(); i++)
	{
		int ind = pt_to_ind->return_enum(pts[i]);
		if (ind >= 0)
		{
			for (unsigned int j = 0; j < lu_table[ind].size(); j++)
			{
				rval.push_back(lu_table[ind][j]);
			}
		}
	}
//	return rval;
}

void DircLUT::get_ckov_theta_all(std::vector<double> &rval, std::vector<dirc_point> pts,double inc_phi, double inc_theta)
{
	std::vector<std::pair<double,double> > table_res;
	get_base_phi_theta_all(table_res, pts);
	
//	std::vector<double> rval;
	rval.clear();
	
	double px,py,pz;

	px = sin(inc_phi/57.3)*sin(inc_theta/57.3);
	py = cos(inc_phi/57.3)*sin(inc_theta/57.3);
	pz = cos(inc_theta/57.3);

	//possible orientations;
	double pxs[8];
	double pys[8];
	double pzs[8];
	int ind = 0;
	//I'm pretty sure there are only 8 geometries to test here.
	for (int i = -1; i <=1; i+=2)
	{
		for (int j = -1; j <=1; j+=2)
		{
			for (int k = -1; k <=1; k+=2)
			{
				pxs[ind] = i*px;
				pys[ind] = j*py;
				pzs[ind] = k*pz;
				//printf("p: %d %12.04f %12.04f %12.04f v: %12.04f %12.04f, %12.04f\n",ind,pxs[j],pys[j],pzs[j]);
				ind++;
			}
		}
	}


	double fill_val = 0;
	double vx,vy,vz;
	double vphi,vtheta;
	for (unsigned int i = 0; i < table_res.size(); i++)
	{
		vphi = table_res[i].first;
		vtheta = table_res[i].second;
		vx = sin(vphi)*sin(vtheta);
		vy = cos(vtheta);
		vz = cos(vphi)*sin(vtheta);
		for (int j = 0; j < 8; j++)
		{
			//printf("p: %12.04f %12.04f %12.04f v: %12.04f %12.04f, %12.04f\n",pxs[j],pys[j],pzs[j],vx,vy,vz);
			fill_val = acos(vx*pxs[j] + vy*pys[j] + vz*pzs[j]);
			rval.push_back(fill_val);
		}	
	}
//	return rval;
}
