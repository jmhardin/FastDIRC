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
	lut_entry add_entry;
	add_entry.phi = phi;
	add_entry.theta = theta;
	add_entry.time = pt.t;

	int ind = pt_to_ind->return_enum(pt);

	if (ind >= 0)
	{
		lu_table[ind].push_back(add_entry);
	}
}
std::vector<lut_entry> DircLUT::get_base_phi_theta(dirc_point pt)
{
	int ind = pt_to_ind->return_enum(pt);

	if (ind >= 0)
	{
		return lu_table[ind];
	}
	else
	{
		// :(  Something is rotten in the state of Denmark
		std::vector<lut_entry> empty;
		empty.clear();
		return empty;
	}
}
void DircLUT::get_base_phi_theta_all(std::vector<lut_entry> &rval, std::vector<dirc_point> pts)
{
	//returns lut_entries with the time "missing" from each
	rval.clear();
	
	for (unsigned int i = 0; i < pts.size(); i++)
	{
		int ind = pt_to_ind->return_enum(pts[i]);
		if (ind >= 0)
		{
			for (unsigned int j = 0; j < lu_table[ind].size(); j++)
			{
				lut_entry push_entry = lu_table[ind][j];
				push_entry.time -= pts[i].t;
				rval.push_back(push_entry);
			}
		}
	}
//	return rval;
}
/*
void DircLUT::get_ckov_theta_all(std::vector<double> &rval, std::vector<dirc_point> pts,double inc_phi, double inc_theta,double inc_y)
{
	//report inc y as distance from top of bar
	double bar_length = 4900; //hardcoded for now, but improving
	double quartz_index = 1.47; //hardcoded for now, but improving
	double c_mm_ns = 300; //hardcoded for now, but improving
	std::vector<lut_entry> table_res;
	get_base_phi_theta_all(table_res, pts);
	
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
	double vt,vt_direct,vt_indirect;
	double y_direct = .5*bar_length - inc_y;
	double y_indirect = 1.5*bar_length + inc_y;

	double internal_refl_limit = sqrt(quartz_index*quartz_index-1)/quartz_index;

	double time_cut = 3;//very loose 3ns cut;
	double fabs_off = 0;//fudge factor;
	for (unsigned int i = 0; i < table_res.size(); i++)
	{
		vphi = table_res[i].phi;
		vtheta = table_res[i].theta;
		vt = table_res[i].time;
		
		vx = sin(vphi)*sin(vtheta);
		vy = cos(vtheta);
		vz = cos(vphi)*sin(vtheta);

		vt_direct = vt + quartz_index*y_direct/(c_mm_ns*vy);		
		vt_indirect = vt + quartz_index*y_indirect/(c_mm_ns*vy);	
		
		//not totally internal reflected - do this cut before hand
		if (vz > internal_refl_limit || vx > internal_refl_limit) continue;

		if (fabs(fabs(vt_direct)-fabs_off) > time_cut && fabs(fabs(vt_indirect)-fabs_off) > time_cut)
		{
			//printf("%12.04f %12.04f %12.04f\n",vt,fabs(vt_direct),fabs(vt_indirect));	
			continue;//Does not pass time cut
		}
		//printf("%12.04f %12.04f %12.04f\n",vt,vt_direct,vt_indirect);	

		for (int j = 0; j < 8; j++)
		{
			//printf("p: %12.04f %12.04f %12.04f v: %12.04f %12.04f, %12.04f\n",pxs[j],pys[j],pzs[j],vx,vy,vz);
			fill_val = acos(vx*pxs[j] + vy*pys[j] + vz*pzs[j]);
			
			rval.push_back(fill_val);
		}	
	}
//	return rval;
}
*/
void DircLUT::get_ckov_theta_all(std::vector<double> &rval, std::vector<double> &ret_dt, std::vector<dirc_point> pts,double inc_phi, double inc_theta,double inc_y)
{
	ret_dt.clear();//syncronized return of time deviations
	rval.clear();//DO NOT FORGET!
	
	//report inc y as distance from top of bar
	double bar_length = 4900; //hardcoded for now, but improving
	double quartz_index = 1.47; //hardcoded for now, but improving
	double c_mm_ns = 300; //hardcoded for now, but improving
	
	double px,py,pz;

	px = sin(inc_phi/57.3)*sin(inc_theta/57.3);
	py = cos(inc_phi/57.3)*sin(inc_theta/57.3);
	pz = cos(inc_theta/57.3);

	//possible orientations;
	double pxs[4];
	double pys[4];
	double pzs[4];
	int ind = 0;
	//I'm pretty sure there are only 8 geometries to test here.
	for (int i = -1; i <=1; i+=2)
	{
		for (int j = -1; j <=1; j+=2)
		{
			pxs[ind] = j*px;
			pys[ind] = py;//account for this based on timing and direction
			pzs[ind] = j*pz;
			//printf("p: %d %12.04f %12.04f %12.04f v: %12.04f %12.04f, %12.04f\n",ind,pxs[j],pys[j],pzs[j]);
			ind++;
		}
	}


	double fill_val = 0;
	double vx,vy,vz;
	double vphi,vtheta;
	double vt,vt_direct,vt_indirect;
	double y_direct = .5*bar_length - inc_y;
	double y_indirect = 1.5*bar_length + inc_y;

	double internal_refl_limit = sqrt(quartz_index*quartz_index-1)/quartz_index;

	double time_cut = 10;//very loose 3ns cut;
	double fabs_off = 0;//fudge factor;
	double angle_mid = 47/57.3;
	double angle_loose = 10/57.3;
	int passed_ind = 0;
	int passed_refl = 0;
	int passed_time = 0;
	int rval_filled = 0;
	for (unsigned int i = 0; i < pts.size(); i++)
	{
		int ind = pt_to_ind->return_enum(pts[i]);
		if (ind < 0) 
		{
			continue;
		}
		double avg_vx=0;
		double avg_vy=0;
		double avg_vz=0;
		double num_in_avg=0;
		
		int time_direction = 0;	
		passed_ind++;
		for (unsigned int j = 0; j < lu_table[ind].size(); j++)
		{
			vphi = lu_table[ind][j].phi;
			vtheta = lu_table[ind][j].theta;
			vt = lu_table[ind][j].time - pts[i].t;

			vx = sin(vphi)*sin(vtheta);
			vy = cos(vtheta);
			vz = cos(vphi)*sin(vtheta);

			vt_direct = vt + quartz_index*y_direct/(c_mm_ns*vy);		
			vt_indirect = vt + quartz_index*y_indirect/(c_mm_ns*vy);	

			//not totally internal reflected - do this cut before hand
			if (vz > internal_refl_limit || vx > internal_refl_limit) 
			{
				continue;
			}

			passed_refl++;

			if (fabs(fabs(vt_direct)-fabs_off) < time_cut)
			{
				time_direction = 1;
			}
			else if (fabs(fabs(vt_indirect)-fabs_off) < time_cut)
			{
				time_direction = -1;
			}
			else
			{
				//printf("%12.04f %12.04f %12.04f\n",vt,fabs(vt_direct),fabs(vt_indirect));	
				continue;//Does not pass time cut
			}
			passed_time++;
			//printf("%12.04f %12.04f %12.04f\n",vt,vt_direct,vt_indirect);	

			for (int k = 0; k < 4; k++)
			{
				//printf("p: %12.04f %12.04f %12.04f v: %12.04f %12.04f, %12.04f\n",pxs[k],pys[k],pzs[k],vx,vy,vz);
				fill_val = acos(vx*pxs[k] + vy*pys[k]*time_direction + vz*pzs[k]);
				if (fabs(fill_val - angle_mid) < angle_loose)
				{
					if (fill_val != fill_val)
					{
						printf("WTH Over.\n");
						continue;
					}
					//sane angle, add it to the average
				//	rval.push_back(fill_val);
					avg_vx += vx;
					avg_vy += vy;
					avg_vz += vz;
					num_in_avg += 1;
					rval.push_back(fill_val);
					rval_filled++;
				
					if (time_direction == 1)
					{
						ret_dt.push_back(vt_direct);
					}
					else if (time_direction == -1)
					{
						ret_dt.push_back(vt_indirect);
					}
	
					//rval.push_back(fill_val);
				}
			}
		}
	}
}
