#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

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

// 	Math::Transform3 trans;  
// 	Math::Vector3 startz(0,0,1);
// 	Math::VectorPair3 pair(0,0,0,0,1,0);
// 	
// 	startz = pair.origin();
// 	printf("pair1 %5.03f %5.03f %5.03f\n",pair.x0(),pair.y0(),pair.z0());
// 	printf("startz %5.03f %5.03f %5.03f\n",startz.x(),startz.y(),startz.z());
// 	
// 	trans.reset();
// 	trans.linear_rotation(Math::Vector3(90,0,0));
// 	trans.linear_rotation(Math::Vector3(0,0,45));
// 	startz = trans.transform_linear(startz);
// 	
// 	printf("startz %5.03f %5.03f %5.03f\n",startz.x(),startz.y(),startz.z());

	double in_num = 0;
	

	if (nargs > 1)
	{
		in_num = atof(argv[1]);
		printf("in number: %12.04f\n",in_num);
	}
	
	double resx = in_num;
	double resy = in_num;
	double minx = -1500;
	double maxx = -minx;
	double miny = -500;
	double maxy = -miny;
	
	double rad_to_deg = 57.2958;
	
	double res_enhance = 1;

	double tracking_unc = .0000*57.3; //mrad
// 	double ckov_unc = .0077*57.3; //chromatic + optical aberation = 7.7mrad
	double ckov_unc = .003*57.3; //transport = 3mrad
		
	double energy = 4.5;
	double kmass = .4937;
	double pimass = .1396;
	
	double particle_x = 0;
	double particle_y = 0;
	double particle_theta = 0;
	double particle_phi = 0;
	
	int num_runs = 4000;
	
	int n_sim_phots = 40;
	
	int n_phi_phots =120000;
	int n_z_phots = 4;	
	double sfunc_m = 25;
	double sfunc_r = 35;
	double sfunc_sig = 6;
	
	//20k,Fold_x, 6, and 1.5 result in .8935 - 3.3 mrad or better
	
	//Good perpendicular version needs up_douwn_sep=false
	
	double pdf_unc_red_fac = 1;
	bool fold_dist = false;
	double wedge_uncertainty = 0/57.3;
	double mirror_angle_change = 0;
	double mirror_angle_change_unc = 0;
	double box_rot = 0;
	double box_rot_unc = 0;
	double bar_box_box_angle = 0/57.3;
	double mirror_r_difference = -10000000000;
	double wedge_non_uniformity = 0;
	
	
	double liquid_absorbtion = 0*-log(.7)/1000;
	double liquid_index = 1.33;
	
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
		pion_cerenkov->Fill(dirc_model->get_cerenkov_angle_rand(pion_beta,ckov_unc));
		kaon_cerenkov->Fill(dirc_model->get_cerenkov_angle_rand(kaon_beta,ckov_unc));
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
	
	for (int i = 0; i < 1; i++)
	{
// 		dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
// 		dirc_model->set_bar_box_angle(bar_box_box_angle);
		hit_points_pion = dirc_model->sim_reg_n_photons(\
			n_phi_phots,\
			n_z_phots,\
			false,\
			false,\
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
			false,\
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

	if (fold_dist == true)
	{
		hit_points_pion = fold_x(hit_points_pion);
		hit_points_kaon = fold_x(hit_points_kaon);
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
	

	
	for (int i = 0; i < num_runs; i++)
	{
// 		dirc_model->set_pmt_angle(spread_ang.Gaus(47.87,box_rot_unc));
		dirc_model->set_focus_mirror_angle(spread_ang.Gaus(74.11,mirror_angle_change_unc));
// 		dirc_model->set_upper_wedge_angle_diff(spread_ang.Gaus(0,wedge_uncertainty));
// 		dirc_model->set_bar_box_angle(spread_ang.Gaus(0,bar_box_box_angle));
		
		printf("\r                                                    ");
		printf("\rrunning iter %d/%d  ",i+1,num_runs);
		
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
		
		digitizer.digitize_points(sim_points);
		
		if (fold_dist == true) sim_points = fold_x(sim_points);
		
		llc = pdf_pion->get_log_likelihood(sim_points);
		llf = pdf_kaon->get_log_likelihood(sim_points);

		ll_diff_pion->Fill(llc-llf);
		
// 		ll_diff_pion->Fill(sep_pdfs->get_log_likelihood_spread_diff(sim_points));
		
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
		
		digitizer.digitize_points(sim_points);
		
		if (fold_dist == true) sim_points = fold_x(sim_points);
		
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
	
	tfile->Close();
	
	int status = 0;
	return status;
}
