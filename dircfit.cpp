#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

#include <stdlib.h>
#include <math.h>

#include "include/dirc_goptical_sim.h"
#include "include/dirc_point.h"
#include "include/dirc_probability_spread.h"
#include "include/dirc_probability_separation.h"
#include "include/dirc_spread_radius.h"
#include "include/dirc_spread_relative.h"
#include "include/dirc_spread_linear_soft.h"
#include "include/dirc_spread_gaussian.h"
#include "include/dirc_digitizer.h"
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TRandom3.h>

#include <Goptical/Math/Vector>
#include <Goptical/Math/VectorPair>
#include <Goptical/Math/Transform>

using namespace Goptical;

std::vector<dirc_point> fold_x(std::vector<dirc_point> inpoints)
{
	std::vector<dirc_point> outvec;
	for (unsigned int i = 0; i < inpoints.size(); i++)
	{
		dirc_point tpoint = inpoints[i];
		tpoint.x = fabs(tpoint.x);
		outvec.push_back(tpoint);
	}
	return outvec;
}

int main(int nargs, char* argv[])
{
	double in_num = 0;
	

	if (nargs > 1)
	{
		in_num = atof(argv[1]);
		printf("in number: %12.04f\n",in_num);
	}
	
	double resx = 6;
	double resy = 6;
	double minx = -8000;
	double maxx = -minx;
	double miny = -800;
	double maxy = -miny;
	
	double rad_to_deg = 57.2958;
	
	double res_enhance = 1;

	double tracking_unc = .0000*57.3; //mrad
// 	double ckov_unc = .0077*57.3; //chromatic + optical aberation = 7.7mrad
	double ckov_unc = .003*57.3; //transport = 3mrad

	
	double energy = 5.0;
	double kmass = .4937;
	double pimass = .1396;
	
	double particle_x = 0;
	double particle_y = 0;
	double particle_theta = 4;
	double particle_phi = 40;
	
	int num_runs = 4000;
	
	int n_sim_phots = 40;
	
	int refraction_sim_n = 0;
	
	int n_phi_phots = 120000;
// 	int n_phi_phots =2;
	int n_z_phots = 4;	
	double sfunc_m = 25;
	double sfunc_r = 35;
	double sfunc_sig = 6;
	
	bool out_layout = false;
	if (out_layout == true)
	{
	  n_phi_phots = 100;
	  n_z_phots = 4;
	}
	
	//20k,Fold_x, 6, and 1.5 result in .8935 - 3.3 mrad or better
	
	//Good perpendicular version needs up_down_sep=false
	
	double pdf_unc_red_fac = 1;
	double wedge_uncertainty = 0/57.3;
	double mirror_angle_change = 0;
	double mirror_angle_change_unc = 0;
	double box_rot = 0;
	double box_rot_unc = 0;
	double bar_box_box_angle = 0/57.3;
	double mirror_r_difference = 0;
	double wedge_non_uniformity = 0;
	double pmt_offset = 0;
	
	double sm_xl = -10000000;
	double sm_xr = -sm_xl;
	
	double overlap_x = -1;
	
	bool three_seg_mirror = true;
	
	
	double liquid_absorbtion = 0*-log(.7)/1000;
	double liquid_index = in_num;
	
	bool coverage_plot = false;
	int num_cov = 100000;
	
	bool up_down_sep = false;
	
	TRandom3 spread_ang(23);
	
	
	DircGopticalSim *dirc_model = new DircGopticalSim(\
		4357,\
		-1200 + mirror_r_difference,\
		300.38,\
		74.11 + mirror_angle_change,\
		600,\
		47.87 + box_rot);
	
	double pion_beta = dirc_model->get_beta(energy,pimass);
	double kaon_beta = dirc_model->get_beta(energy,kmass);
	
	double pion_angle = rad_to_deg*acos(1/(1.47*pion_beta));
	double kaon_angle = rad_to_deg*acos(1/(1.47*kaon_beta));
	
	printf("Pion emission angle nominal: %12.02f\n",pion_angle);
	printf("Kaon emission angle nominal: %12.02f\n",kaon_angle);
	
	
	char* rootfilename = new char[256];
	sprintf(rootfilename,"fitdirc.root");	

	TFile* tfile = new TFile(rootfilename,"RECREATE");
	
	TH1F *ll_diff_pion = new TH1F("ll_diff_pion","difference of log likelihood real = pion",8000,-1000,1000);
	TH1F *ll_diff_kaon = new TH1F("ll_diff_kaon","difference of log likelihood real = kaon",8000,-1000,1000);
	
	TH1F *phot_found_pion = new TH1F("phot_found_pion","number of photons found on pion angle", 1001,-.5,1000.5);
	TH1F *phot_found_kaon = new TH1F("phot_found_kaon","number of photons found on kaon angle", 1001,-.5,1000.5);
	
	TH1F *ref_pion_before = new TH1F("ref_pion_before","Angle of Pion Photons going into interface", 9000,0,90);
	TH1F *ref_kaon_before = new TH1F("ref_kaon_before","Angle of Kaon Photons going into interface", 9000,0,90);
	TH1F *ref_pion_after = new TH1F("ref_pion_after","Angle of Pion Photons going out of interface", 9000,0,90);
	TH1F *ref_kaon_after = new TH1F("ref_kaon_after","Angle of Kaon Photons going out of interface", 9000,0,90);

	TH2F *ref_theta_cphi_pion = new TH2F("ref_theta_cphi_pion","Emitted angle of Pion Photons versus angle into interface", 7200,0,360,1800,0,90);
	TH2F *ref_theta_cphi_kaon = new TH2F("ref_theta_cphi_kaon","Emitted angle of Kaon Photons versus angle into interface", 7200,0,360,1800,0,90);
	
	TH2F *pion_dist_xy = new TH2F("pion_dist_xy","xy val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH2F *kaon_dist_xy = new TH2F("kaon_dist_xy","xy val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

	maxy *= 5;
	
	TH2F *pion_coverage_xy = new TH2F("pion_coverage_xy","xy val of generated points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH2F *kaon_coverage_xy = new TH2F("kaon_coverage_xy","xy val of generated points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

	
	TH1F *pion_cerenkov = new TH1F("pion_cerenkov","pion cerenkov angle distribution",300,45,48);
	TH1F *kaon_cerenkov = new TH1F("kaon_cerenkov","kaon cerenkov angle distribution",300,45,48);

	TH1F *liquid_dist = new TH1F("liquid_dist","distance travelled in liquid (mm)",1500,0,1500);
	dirc_model ->set_store_traveled(true);
	
	for (int i = 0; i < 10; i++)
	{
		double tmp;
		pion_cerenkov->Fill(dirc_model->get_cerenkov_angle_rand(pion_beta,ckov_unc, tmp));
		kaon_cerenkov->Fill(dirc_model->get_cerenkov_angle_rand(kaon_beta,ckov_unc, tmp));
	}

	
	
	
	DircDigitizer digitizer(\
		minx,\
		maxx,\
		resx,\
		miny,\
		maxy,\
		resy);
	std::vector<dirc_point> hit_points_pion;
	std::vector<dirc_point> hit_points_kaon;
	
// 	dirc_model->set_wedge_mirror_rand(wedge_non_uniformity);

// 	dirc_model->set_upper_wedge_angle_diff(0);
	
	dirc_model->set_liquid_absorbtion(liquid_absorbtion);
	dirc_model->set_liquid_index(liquid_index);
	dirc_model->set_three_seg_mirror(three_seg_mirror);
	dirc_model->set_sidemirror(sm_xr,sm_xl);
	
	
	dirc_model->set_pmt_offset(pmt_offset);
	
	for (int i = 0; i < 1; i++)
	{
// 		dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
// 		dirc_model->set_bar_box_angle(bar_box_box_angle);
		hit_points_pion = dirc_model->sim_reg_n_photons(\
			n_phi_phots,\
			n_z_phots,\
			false,\
			out_layout,\
			pion_angle,\
			particle_x,\
			particle_y,\
			particle_theta,\
			particle_phi,\
			0,\
			ckov_unc/pdf_unc_red_fac,\
			pion_beta,\
			up_down_sep);

		hit_points_kaon = dirc_model->sim_reg_n_photons(\
			n_phi_phots,\
			n_z_phots,\
			false,\
			out_layout,\
			kaon_angle,
			particle_x,\
			particle_y,\
			particle_theta,\
			particle_phi,\
			0,\
			ckov_unc/pdf_unc_red_fac,\
			kaon_beta,\
			up_down_sep);
	}

// 	DircProbabilitySpread* pdf_pion = new DircSpreadLinearSoft(sfunc_m, sfunc_r, sfunc_sig, hit_points_pion);
// 	DircProbabilitySpread* pdf_kaon = new DircSpreadLinearSoft(sfunc_m, sfunc_r, sfunc_sig, hit_points_kaon);

	DircProbabilitySpread* pdf_pion = new DircSpreadGaussian(sfunc_sig, hit_points_pion);
	DircProbabilitySpread* pdf_kaon = new DircSpreadGaussian(sfunc_sig, hit_points_kaon);
// 	
	//Digitizing and linear is eqivalent to "counting," but slower
// 	digitizer.digitize_points(hit_points_pion);
// 	digitizer.digitize_points(hit_points_kaon);
// 	DircProbabilitySpread* pdf_pion = new DircSpreadLinearSoft(sfunc_m, sfunc_r, sfunc_sig, hit_points_pion);
//  DircProbabilitySpread* pdf_kaon = new DircSpreadLinearSoft(sfunc_m, sfunc_r, sfunc_sig, hit_points_kaon);

	DircProbabilitySeparation* sep_pdfs = new DircProbabilitySeparation(pdf_kaon,pdf_pion);
	
	printf("Begining Run\n");
	double llc, llf;
	std::vector<dirc_point> sim_points;
	std::vector<dirc_point> confound_points;

	
	for (int i = 0; i < num_runs; i++)
	{
// 		dirc_model->set_pmt_angle(spread_ang.Gaus(47.87,box_rot_unc));
		dirc_model->set_focus_mirror_angle(spread_ang.Gaus(74.11,mirror_angle_change_unc));
// 		dirc_model->set_upper_wedge_angle_diff(spread_ang.Gaus(0,wedge_uncertainty));
// 		dirc_model->set_bar_box_angle(spread_ang.Gaus(0,bar_box_box_angle));
		
		printf("\r                                                    ");
		printf("\rrunning iter %d/%d  ",i+1,num_runs);
		
		if (overlap_x > 0)
		{
			confound_points = dirc_model->sim_rand_n_photons(\
			      n_sim_phots,\
			      false,\
			      false,\
			      pion_angle,\
			      particle_x,\
			      particle_y,\
			      particle_theta,\
			      particle_phi,\
			      tracking_unc,\
			      ckov_unc,\
			      1,\
			      up_down_sep);
			//assume electrom has v=c
			
			for (unsigned int i = 0; i < confound_points.size(); i++)
			{
			      confound_points[i].x += overlap_x;
			}
		}
		
		
		fflush(stdout);
		sim_points = dirc_model->sim_rand_n_photons(\
			n_sim_phots,\
			false,\
			false,\
			pion_angle ,
			particle_x,\
			particle_y,\
			particle_theta,\
			particle_phi,\
			tracking_unc,\
			ckov_unc,\
			pion_beta,\
			up_down_sep);
		for (unsigned int i = 0; i < confound_points.size(); i++)
		{
			sim_points.push_back(confound_points[i]);
		}
		digitizer.digitize_points(sim_points);
		
		llc = pdf_pion->get_log_likelihood(sim_points);
		llf = pdf_kaon->get_log_likelihood(sim_points);
		

		ll_diff_pion->Fill(llc-llf);
				
		phot_found_pion->Fill(sim_points.size());
		
		sim_points = dirc_model->sim_rand_n_photons(\
			n_sim_phots,\
			false,\
			false,\
			kaon_angle,\
			particle_x,\
			particle_y,\
			particle_theta,\
			particle_phi,\
			tracking_unc,\
			ckov_unc,\
			kaon_beta,\
			up_down_sep);
		
		for (unsigned int i = 0; i < confound_points.size(); i++)
		{
			sim_points.push_back(confound_points[i]);
		}
		
		digitizer.digitize_points(sim_points);
				
		llc = pdf_pion->get_log_likelihood(sim_points);
		llf = pdf_kaon->get_log_likelihood(sim_points);

		ll_diff_kaon->Fill(llc-llf);
		
// 		ll_diff_kaon->Fill(sep_pdfs->get_log_likelihood_spread_diff(sim_points));
		
		phot_found_kaon->Fill(sim_points.size());
	}
	
	printf("\nRun Completed\n");
	if (coverage_plot == true)
	{
		printf("starting coverage ouput\n");

		for (int i = 0; i < num_cov; i++)
		{
			printf("\r                                                                                                                           ");
			printf("\rcoverage iter %d/%d",i+1,num_cov);
			
			fflush(stdout);
			
			particle_theta = spread_ang.Uniform(0,14);
			particle_phi = spread_ang.Uniform(0,360);
			energy = spread_ang.Uniform(2,4.5);
			pion_beta = dirc_model->get_beta(energy,pimass);
			kaon_beta = dirc_model->get_beta(energy,kmass);
			
			sim_points = dirc_model->sim_rand_n_photons(\
				n_sim_phots,\
				false,\
				false,\
				pion_angle ,
				particle_x,\
				particle_y,\
				particle_theta,\
				particle_phi,\
				tracking_unc,\
				ckov_unc,\
				pion_beta,\
				up_down_sep);
			
			for (unsigned int j = 0; j < sim_points.size(); j++)
			{
				pion_coverage_xy->Fill(sim_points[j].x,sim_points[j].y);
			}
			
			sim_points = dirc_model->sim_rand_n_photons(\
				n_sim_phots,\
				false,\
				false,\
				pion_angle ,
				particle_x,\
				particle_y,\
				particle_theta,\
				particle_phi,\
				tracking_unc,\
				ckov_unc,\
				kaon_beta,\
				up_down_sep);
			
			for (unsigned int j = 0; j < sim_points.size(); j++)
			{
				kaon_coverage_xy->Fill(sim_points[j].x,sim_points[j].y);
			}
			
		}
		printf("\nfinished output of coverage plots\n");
	}
	else
	{
		pion_coverage_xy->Fill(0.0,0.0);
		kaon_coverage_xy->Fill(0.0,0.0);

	}
	
	printf("Found %d pion points on the target\n", (int) hit_points_pion.size());
	double x,y;
	for (unsigned int i = 0; i < hit_points_pion.size(); i++)
	{
		x = hit_points_pion[i].x;
		y = hit_points_pion[i].y;
// 		if (hit_points_pion[i].t > 0) continue;
		pion_dist_xy->Fill(x,y);
	}
	
	printf("Found %d kaon points on the target\n", (int) hit_points_kaon.size());
	for (unsigned int i = 0; i < hit_points_kaon.size(); i++)
	{
		x = hit_points_kaon[i].x;
		y = hit_points_kaon[i].y;
// 		if (hit_points_pion[i].t > 0) continue;
		kaon_dist_xy->Fill(x,y);
	}
	
	std::vector<double> dist_traveled = dirc_model->get_dist_traveled();
	for (unsigned int i = 0; i < dist_traveled.size(); i++)
	{
		liquid_dist->Fill(dist_traveled[i]);
	}
	
	if (refraction_sim_n > 0)
	{
	  std::vector<std::pair<double, double> > pion_theta_cphi;
	  std::vector<std::pair<double, double> > kaon_theta_cphi;
	  
	  
	  printf("Creating Refraction histograms...\n");
	  std::vector<double> before_fill;
	  std::vector<double> after_fill;
	  
	  pion_theta_cphi = dirc_model->get_refraction_rand_phi(\
		before_fill,\
		after_fill,\
		refraction_sim_n,\
		pion_angle ,
		particle_x,\
		particle_y,\
		particle_theta,\
		particle_phi,\
		tracking_unc,\
		ckov_unc,\
		0,\
		pion_beta);
	  
	  
	  for (unsigned int i = 0; i < before_fill.size(); i++)
	  {
	    ref_pion_before->Fill(before_fill[i]*rad_to_deg);
	    ref_pion_after->Fill(after_fill[i]*rad_to_deg);
	  }
	  for (unsigned int i = 0; i < pion_theta_cphi.size(); i++)
	  {
	    ref_theta_cphi_pion->Fill(pion_theta_cphi[i].first*rad_to_deg,pion_theta_cphi[i].second*rad_to_deg);
	  }
	  
	  kaon_theta_cphi = dirc_model->get_refraction_rand_phi(\
		before_fill,\
		after_fill,\
		refraction_sim_n,\
		kaon_angle ,
		particle_x,\
		particle_y,\
		particle_theta,\
		particle_phi,\
		tracking_unc,\
		ckov_unc,\
		0,\
		kaon_beta);
	  
	  
	  for (unsigned int i = 0; i < before_fill.size(); i++)
	  {
	    ref_kaon_before->Fill(before_fill[i]*rad_to_deg);
	    ref_kaon_after->Fill(after_fill[i]*rad_to_deg);
	  }
	  for (unsigned int i = 0; i < kaon_theta_cphi.size(); i++)
	  {
	    ref_theta_cphi_kaon->Fill(kaon_theta_cphi[i].first*rad_to_deg,kaon_theta_cphi[i].second*rad_to_deg);
	  }
	  printf("Refraction histograms created\n");
	}
	
	
	pion_dist_xy->Write();
	kaon_dist_xy->Write();
	ll_diff_pion->Write();
	ll_diff_kaon->Write();
	phot_found_pion->Write();
	phot_found_kaon->Write();
	pion_cerenkov->Write();
	kaon_cerenkov->Write();
	pion_coverage_xy->Write();
	kaon_coverage_xy->Write();
	liquid_dist->Write();
	
	ref_pion_before->Write();
	ref_pion_after->Write();
	ref_kaon_before->Write();
	ref_kaon_after->Write();
	
	
	ref_theta_cphi_pion->Write();
	ref_theta_cphi_kaon->Write();
	tfile->Close();
	
	int status = 0;
	return status;
}
