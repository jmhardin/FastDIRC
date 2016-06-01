#include "../include/dirc_lut_enum.h"
#include "../include/dirc_lut.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
void DircLUT::get_ckov_theta_all(std::vector<double> &rval, std::vector<double> &ret_dt_dl, std::vector<dirc_point> pts,double inc_phi, double inc_theta,double inc_y)
{
	ret_dt_dl.clear();//syncronized return of time deviations divided by propagation distance in ps/mm
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
	for (int i = -1; i <=1; i+=2)
	{
		for (int j = -1; j <=1; j+=2)
		{
			pxs[ind] = i*px;
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

	double time_cut = 10;//very loose 10ns cut;
	double fabs_off = 0;//fudge factor;
	double angle_mid = 47.1/57.3;
	double angle_loose = 3/57.3;
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
		bool added_angle = false;
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
			//vt_direct = vt + quartz_index*y_direct/(c_mm_ns*vy);		
			//t_indirect = vt + quartz_index*y_indirect/(c_mm_ns*vy);	

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
			double avg_angle_pt = 0;
			double avg_angle_count = 0;
			for (int k = 0; k < 4; k++)
			{
				//printf("p: %12.04f %12.04f %12.04f v: %12.04f %12.04f, %12.04f\n",pxs[k],pys[k],pzs[k],vx,vy,vz);
				fill_val = acos(vx*pxs[k] + vy*pys[k]*time_direction + vz*pzs[k]);

				//printf("%12.04f\n",fill_val);
				if (fabs(fill_val - angle_mid) < angle_loose)
				{
					if (fill_val != fill_val)
					{
						printf("WTH Over.\n");
						continue;
					}
					avg_vx += vx;
					avg_vy += vy;
					avg_vz += vz;
					num_in_avg += 1;

					avg_angle_pt += fill_val;
					avg_angle_count += 1;

					rval.push_back(fill_val);
					//printf("%10d %12.04f\n",ind,fill_val*57.3);
					//rval_filled++;

					if (time_direction == 1)
					{
						//ret_dt_dl.push_back(vt_direct);
						ret_dt_dl.push_back(vt_direct*vy/y_direct*1000);
					}
					else if (time_direction == -1)
					{
						//ret_dt_dl.push_back(vt_indirect);
						ret_dt_dl.push_back(vt_indirect*vy/y_indirect*1000);
					}
					added_angle = true; 
					//break;
				}
			}
		}
	}
	//printf("%d / %d\n",rval.size(),pts.size());
}
void DircLUT::get_chromcorr_m_b_single_oval_cut(
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
		double &b_indirect)
{
	std::vector<double> dtdl_direct;
	std::vector<double> dang_direct;
	std::vector<double> dtdl_indirect;
	std::vector<double> dang_indirect;

	//center_ang_spread_sq /= 12;
	//time_spread_sq /= 12;

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
			pxs[ind] = i*px;
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
	int passed_ind = 0;
	int passed_refl = 0;
	int passed_time = 0;
	int rval_filled = 0;
	double used_dt = 0;
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
		bool added_angle = false;
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

			if (fabs(fabs(vt_direct)) < time_cut)
			{
				time_direction = 1;
				used_dt = vt_direct;
			}
			else if (fabs(fabs(vt_indirect)) < time_cut)
			{
				time_direction = -1;
				used_dt = vt_indirect;
			}
			else
			{
				//printf("%12.04f %12.04f %12.04f\n",vt,fabs(vt_direct),fabs(vt_indirect));	
				continue;//Does not pass time cut
			}
			passed_time++;
			//printf("%12.04f %12.04f %12.04f\n",vt,vt_direct,vt_indirect);	
			double avg_angle_pt = 0;
			double avg_angle_count = 0;
			for (int k = 0; k < 4; k++)
			{
				//printf("p: %12.04f %12.04f %12.04f v: %12.04f %12.04f, %12.04f\n",pxs[k],pys[k],pzs[k],vx,vy,vz);
				fill_val = acos(vx*pxs[k] + vy*pys[k]*time_direction + vz*pzs[k]);
				if ((fill_val-center_ang)*(fill_val-center_ang)/center_ang_spread_sq + used_dt*used_dt/time_spread_sq > 1)
				{
					continue;
				}

				if (fill_val != fill_val)
				{
					printf("WTH Over.\n");
					continue;
				}
				avg_vx += vx;
				avg_vy += vy;
				avg_vz += vz;
				num_in_avg += 1;

				avg_angle_pt += fill_val;
				avg_angle_count += 1;

				added_angle = true; 
				//break;
			}
			if (avg_angle_count > 0)
			{
				if (time_direction == 1)
				{
					dtdl_direct.push_back(vt_direct*vy/y_direct);//note that this is not scaled for display benefit
					dang_direct.push_back(avg_angle_pt/avg_angle_count-expected_angle);//need to "know" angle to do correction
				}
				else if (time_direction == -1)
				{
					dtdl_indirect.push_back(vt_indirect*vy/y_indirect);//note that this is not scaled for display benefit
					dang_indirect.push_back(avg_angle_pt/avg_angle_count-expected_angle);//need to "know" angle to do correction
				}
			}
		}
	}
	//now do the linear regression
	double mean_dtdl_direct=0;
	double mean_dang_direct=0;
	double mean_dtdl_indirect=0;
	double mean_dang_indirect=0;
	for (unsigned int i = 0; i < dtdl_direct.size(); i++)
	{
		mean_dtdl_direct += dtdl_direct[i];
	}
	mean_dtdl_direct /= dtdl_direct.size();
	for (unsigned int i = 0; i < dang_direct.size(); i++)
	{
		mean_dang_direct += dang_direct[i];
	}
	mean_dang_direct /= dang_direct.size();
	for (unsigned int i = 0; i < dtdl_indirect.size(); i++)
	{
		mean_dtdl_indirect += dtdl_indirect[i];
	}
	mean_dtdl_indirect /= dtdl_indirect.size();
	for (unsigned int i = 0; i < dang_indirect.size(); i++)
	{
		mean_dang_indirect += dang_indirect[i];
	}
	mean_dang_indirect /= dang_indirect.size();


	double slope_direct = 0;
	double intercept_direct = 0;
	double slope_indirect = 0;
	double intercept_indirect = 0;

	double dtdl_dev2_direct = 0; //Sum[(x-meanx)^2]
	double dtdldang_dev_direct = 0;	//Sum[(x-meanx)(y-meany)]
	double dtdl_dev2_indirect = 0; //Sum[(x-meanx)^2]
	double dtdldang_dev_indirect = 0;	//Sum[(x-meanx)(y-meany)]

	for (unsigned int i = 0; i < dtdl_direct.size(); i++)
	{
		dtdl_dev2_direct += (dtdl_direct[i] - mean_dtdl_direct)*(dtdl_direct[i] - mean_dtdl_direct);
		dtdldang_dev_direct += (dang_direct[i] - mean_dang_direct)*(dtdl_direct[i] - mean_dtdl_direct);
	}
	for (unsigned int i = 0; i < dtdl_indirect.size(); i++)
	{
		dtdl_dev2_indirect += (dtdl_indirect[i] - mean_dtdl_indirect)*(dtdl_indirect[i] - mean_dtdl_indirect);
		dtdldang_dev_indirect += (dang_indirect[i] - mean_dang_indirect)*(dtdl_indirect[i] - mean_dtdl_indirect);
	}
	m_direct = dtdldang_dev_direct/dtdl_dev2_direct;
	b_direct = mean_dang_direct - m_direct*mean_dtdl_direct;
	m_indirect = dtdldang_dev_indirect/dtdl_dev2_indirect;
	b_indirect = mean_dang_indirect - m_indirect*mean_dtdl_indirect;
}
void DircLUT::get_ckov_theta_single_oval_cut(
		std::vector<double> &rval, \
		std::vector<double> &ret_dt_dl, \
		std::vector<dirc_point> pts, \
		double inc_phi, \
		double inc_theta, \
		double inc_y, \
		double center_ang, \
		double center_ang_spread_sq,\
		double time_spread_sq,\
		double m_direct,/*=0*/\
		double b_direct,/*=0*/\
		double m_indirect,/*=0*/\
		double b_indirect/*=0*/)
{
	ret_dt_dl.clear();//syncronized return of time deviations
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
			pxs[ind] = i*px;
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
	int passed_ind = 0;
	int passed_refl = 0;
	int passed_time = 0;
	int rval_filled = 0;
	double used_dt = 0;
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
		bool added_angle = false;
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

			if (fabs(fabs(vt_direct)) < time_cut)
			{
				time_direction = 1;
				used_dt = vt_direct;
			}
			else if (fabs(fabs(vt_indirect)) < time_cut)
			{
				time_direction = -1;
				used_dt = vt_indirect;
			}
			else
			{
				//printf("%12.04f %12.04f %12.04f\n",vt,fabs(vt_direct),fabs(vt_indirect));	
				continue;//Does not pass time cut
			}
			passed_time++;
			//printf("%12.04f %12.04f %12.04f\n",vt,vt_direct,vt_indirect);	
			double avg_angle_pt = 0;
			double avg_angle_count = 0;
			double chrom_corr = 0;
			for (int k = 0; k < 4; k++)
			{
				//printf("p: %12.04f %12.04f %12.04f v: %12.04f %12.04f, %12.04f\n",pxs[k],pys[k],pzs[k],vx,vy,vz);
				fill_val = acos(vx*pxs[k] + vy*pys[k]*time_direction + vz*pzs[k]);
				if ((fill_val-center_ang)*(fill_val-center_ang)/center_ang_spread_sq + used_dt*used_dt/time_spread_sq > 1)
				{
					continue;
				}

				if (fill_val != fill_val)
				{
					printf("WTH Over.\n");
					continue;
				}
				avg_vx += vx;
				avg_vy += vy;
				avg_vz += vz;
				num_in_avg += 1;

				avg_angle_pt += fill_val;
				avg_angle_count += 1;

				added_angle = true; 
				//break;
			}
			if (avg_angle_count > 0)
			{
				//if (added_angle == true) break;
				if (time_direction == 1)
				{
					//ret_dt.push_back(vt_direct);
					ret_dt_dl.push_back(vt_direct*vy/y_direct*1000);
					chrom_corr = m_direct*vt_direct*vy/y_direct + b_direct;
				}
				else if (time_direction == -1)
				{
					//ret_dt.push_back(vt_indirect);
					ret_dt_dl.push_back(vt_indirect*vy/y_indirect*1000);
					chrom_corr = m_indirect*vt_indirect*vy/y_indirect + b_indirect;
				}
				rval.push_back(avg_angle_pt/avg_angle_count - chrom_corr);
			}
		}
	}
}
