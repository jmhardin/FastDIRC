#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "include/dirc_optical_sim.h"
#include "include/dirc_point.h"
#include "include/dirc_probability_spread.h"
#include "include/dirc_probability_separation.h"
#include "include/dirc_spread_radius.h"
#include "include/dirc_spread_relative.h"
#include "include/dirc_spread_linear_soft.h"
#include "include/dirc_spread_gaussian.h"
#include "include/dirc_digitizer.h"
#include "include/dirc_progressive_separation.h"
#include <TFile.h>
#include <TTree.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TRandom3.h>

std::vector<dirc_point> fold_x(std::vector<dirc_point> inpoints) {
	std::vector<dirc_point> outvec;
	for (unsigned int i = 0; i < inpoints.size(); i++) {
		dirc_point tpoint = inpoints[i];
		tpoint.x = fabs(tpoint.x);
		outvec.push_back(tpoint);
	}
	return outvec;
}
int main(int nargs, char* argv[])
{  
	clock_t timing_clock;

	const char* in_str;
	bool inputfile = false;
	bool inputrootfile = false;
	bool out_csv = false;
	bool slac_run = false;
	double time_window=-1;//time window for confounded pmt hits, in ns	

	double energy = 5.0;
	double energy_mean = energy;
	double energy_spread = 0;
	double kmass = .4937;
	double pimass = .1396;
	double mumass = .1057;

	double particle_x = 0;
	double particle_y = 0;
	double particle_x_mean = particle_x;
	double particle_y_mean = particle_y;
	double particle_x_spread = 0;
	double particle_y_spread = 0;
	double particle_theta = 4;
	double particle_theta_mean = particle_theta;
	double particle_theta_spread = 0;
	double particle_phi = 40;
	double particle_phi_mean = particle_phi;
	double const_track_off = 0;

	double particle_flight_distance = 0;

	int box_check_n = -1;
	double box_check_theta = 0;
	double box_check_phi = 0;
	double box_check_theta_unc = 0;
	double box_check_phi_unc = 0;
	double box_check_overall_theta = -15;
	int phot_check_n = -1;
	double phot_check_max_theta = 4;


	bool force_kinematics = false;
	bool use_prog_sep = false;
	bool kaleidoscope_plot = false;
	bool monochrome_plot = false;
	bool flatten_time = false;
	bool sep_updown = false;

	int num_runs = 1000;
	int max_particles = 1000000;
	int phi_phots_reduce = 1;
	int refraction_sim_n = -1;

	double mean_n_phot = 40;
	double spread_n_phot = 0;

	double wedge_uncertainty = 0/57.3;
	double refrac_index=1.47;
	double mirror_angle_change = 0;
	double mirror_angle_change_unc = 0;
	double mirror_angle_change_yunc = 0;
	double box_rot = 0;
	double box_rot_unc = 0;
	double bar_box_box_angle = 0/57.3;
	double mirror_r_difference = -400;
	//	double mirror_r_difference = 0;
	double wedge_non_uniformity = 0;
	double pmt_offset = 0;
	double main_mirror_nonuniformity = 0;


	double pmt_min_z = -1000;
	double pmt_max_z = 1000;
	double large_mirror_min_z = -1000;
	double large_mirror_max_z = 1000;

	double upper_wedge_yang_spread = 0;
	int rseed = 1337;

	int broaden_events = 0;

	double tracking_unc = .0000*57.3; //mrad
	// 	double ckov_unc = .0077*57.3; //chromatic + optical aberation = 7.7mrad
	double ckov_unc = .003*57.3; //transport = 3mrad

	double resx = 6;
	double resy = 6;
	double rest = 1;
	//	double minx = -8000;
	//	double maxx = -minx;
	//	double miny = -800;
	//	double maxy = -miny;
	double minx = -1500;
	double maxx = 1500;
	double miny = -1500;
	double maxy = 1500;
	double mint = 0;
	double maxt = 1000;
	double t_unc = .27;
	double t_bin_size = 1;

	double digit_miny = -50;
	double digit_maxy = 300;

	digit_miny = miny;
	digit_maxy = maxy;



	//Sets the side boundarys of the distributions
	double sm_xl = -10000000;
	double sm_xr = -sm_xl;

	double s_func_x = 6;
	double s_func_y = s_func_x;
	double s_func_t = 2;
	double sfunc_sig = 1;

	int n_sim_phots = 40;

	int n_phi_phots = 150000;
	int n_z_phots = 4;
	int n_step_phots = 1000;
	// 	n_step_phots = n_z_phots*n_phi_phots;

/*
	sm_xl = -50;
	sm_xr = sm_xl + 440;
*/
//	sm_xl = 0;
//	sm_xr = sm_xl + 1100;

//	double overlap_x = -1;

	bool three_seg_mirror = true;
	bool fill_distributions	= false;
	bool fill_kinematics_yields = false;
	double kinematics_yields_min_theta = 0;
	double kinematics_yields_max_theta = 12;
	double kinematics_yields_step_theta = .25;
	double kinematics_yields_min_phi = 0;
	double kinematics_yields_max_phi = 45;
	double kinematics_yields_step_phi = 1;
	int  kinematics_n_phots = 100000;


	double liquid_absorbtion = 0*-log(.7)/1000;
	double liquid_index = 1.33;

	bool coverage_plot = false;
	int num_cov = 100000;

	char* rootfilename = new char[256];
	char* inputrootfilename = new char[256];
	sprintf(rootfilename,"fitdirc.root");	

	printf("Arguments Passed=%d\n",nargs);

	if(nargs==2){
		in_str = argv[1];
		printf("%s\n",in_str);
		printf("nargs=2\n");
		inputfile = true;
	}
	else{
		for (int i = 1; i < nargs; i++)
		{
			if (strcmp(argv[i], "-if") == 0)
			{
				i++;
				in_str = argv[i];
				inputfile = true;
			}
			else if (strcmp(argv[i], "-of") == 0)
			{
				i++;
				sprintf(rootfilename,"%s",argv[i]);	
			}
			else if (strcmp(argv[i], "-fill_dist") == 0)
			{
				fill_distributions = true;	
			}
			else if (strcmp(argv[i], "-cylindrical_mirror") == 0)
			{
				three_seg_mirror = false;	
			}
			else if (strcmp(argv[i], "-updown") == 0)
			{
				printf("Up-Down Histogram filling only implemented in loop mode, option currently does nothing in other modes\n");
				sep_updown = true;	
			}
			else if (strcmp(argv[i], "-root_input") == 0)
			{
				i++;
				inputrootfile = true;
				sprintf(inputrootfilename,"%s",argv[i]);
			}
			else if (strcmp(argv[i], "-particle_phi") == 0)
			{
				i++;
				particle_phi = atof(argv[i]);
				particle_phi_mean = particle_phi;
			}
			else if (strcmp(argv[i], "-particle_flight_distance") == 0)
			{
				//meters
				i++;
				particle_flight_distance = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-tracking_unc") == 0)
			{
				i++;
				tracking_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-const_track_off") == 0)
			{
				i++;
				const_track_off = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-ckov_unc") == 0)
			{
				i++;
				ckov_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-refraction_sim_n") == 0)
			{
				i++;
				refraction_sim_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_theta") == 0)
			{
				i++;
				particle_theta = atof(argv[i]);
				particle_theta_mean = particle_theta;
			}
			else if (strcmp(argv[i], "-force_kinematics") == 0)
			{
				force_kinematics = true;
			}
			else if (strcmp(argv[i], "-fill_kinematics_yields") == 0)
			{
				fill_kinematics_yields = true;
			}
			else if (strcmp(argv[i], "-coverage_plot") == 0)
			{
				coverage_plot = true;
			}
			else if (strcmp(argv[i], "-kaleidoscope_plot") == 0)
			{
				kaleidoscope_plot = true;
			}
			else if (strcmp(argv[i], "-monochrome_light") == 0)
			{
				monochrome_plot = true;
			}
			else if (strcmp(argv[i], "-prog_sep") == 0)
			{
				use_prog_sep = true;
			}
			else if (strcmp(argv[i], "-slac_run") == 0)
			{
				slac_run = true;
				three_seg_mirror = false;
				mirror_r_difference = 0;
				//mean_n_phot = 31.1;
				mean_n_phot = 32.4;
				spread_n_phot = 6;
				liquid_index = 1.47;
				phi_phots_reduce = 10;
			}
			else if (strcmp(argv[i], "-slac_geometry") == 0)
			{
				three_seg_mirror = false;
				mirror_r_difference = 0;
				//mean_n_phot = 31.1;
				mean_n_phot = 32.4;
				spread_n_phot = 6;
				liquid_index = 1.47;
				phi_phots_reduce = 10;
				//	sm_xl = -300;
				//	sm_xr = sm_xl + 1000;
				miny = -1500;
				maxy = 1500;
				digit_miny = miny;
				digit_maxy = maxy;
			}
			else if (strcmp(argv[i], "-flatten_time") == 0)
			{
				flatten_time = true;
			}
			else if (strcmp(argv[i], "-open_image_plane") == 0)
			{
				pmt_min_z = -559;
				pmt_max_z = -329;
				large_mirror_min_z = -559;
				large_mirror_max_z = -130;
			}
			else if (strcmp(argv[i], "-mean_n_phot") == 0)
			{
				i++;
				mean_n_phot = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-s_func_t") == 0)
			{
				i++;
				s_func_t = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-t_unc") == 0)
			{
				i++;
				t_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-t_bin_size") == 0)
			{
				i++;
				t_bin_size = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-n_phi_phots") == 0)
			{
				i++;
				n_phi_phots = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-n_z_phots") == 0)
			{
				i++;
				n_z_phots = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-spread_n_phot") == 0)
			{
				i++;
				spread_n_phot = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-t") == 0)
			{
				i++;
				time_window = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_x_spread") == 0)
			{
				i++;
				particle_x_spread = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_y_spread") == 0)
			{
				i++;
				particle_y_spread = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-particle_theta_spread") == 0)
			{
				i++;
				particle_theta_spread = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-energy_spread") == 0)
			{
				i++;
				energy_spread = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-E") == 0)
			{
				i++;
				energy = atof(argv[i]);
				energy_mean = energy;
			}
			else if (strcmp(argv[i], "-n") == 0)
			{
				i++;
				num_runs = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-max_particles") == 0)
			{
				i++;
				max_particles = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-pmt_res") == 0)
			{
				i++;
				resx = atof(argv[i]);
				resy = resx;
			}
			else if (strcmp(argv[i], "-liquid_index") == 0)
			{
				i++;
				liquid_index = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-broaden_events") == 0)
			{
				i++;
				broaden_events = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-wedge_uncertainty") == 0)
			{
				i++;
				wedge_uncertainty = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-refrac_index") == 0)
			{
				i++;
				refrac_index = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-mirror_angle_change") == 0)
			{
				i++;
				mirror_angle_change = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-mirror_angle_change_unc") == 0)
			{
				i++;
				mirror_angle_change_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-mirror_angle_change_yunc") == 0)
			{
				i++;
				mirror_angle_change_yunc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_rot") == 0)
			{
				i++;
				box_rot = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_rot_unc") == 0)
			{
				i++;
				box_rot_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-bar_box_box_angle") == 0)
			{
				i++;
				bar_box_box_angle = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-mirror_r_difference") == 0)
			{
				i++;
				mirror_r_difference = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-wedge_non_uniformity") == 0)
			{
				i++;
				wedge_non_uniformity = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-pmt_offset") == 0)
			{
				i++;
				pmt_offset = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-main_mirror_nonuniformity") == 0)
			{
				i++;
				main_mirror_nonuniformity = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-upper_wedge_yang_spread") == 0)
			{
				i++;
				upper_wedge_yang_spread = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-rseed") == 0)
			{
				i++;
				rseed = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-upper_wedge_yang_spread") == 0)
			{
				i++;
				upper_wedge_yang_spread = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_check_n") == 0)
			{
				i++;
				box_check_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-box_check_theta") == 0)
			{
				i++;
				box_check_theta = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_check_phi") == 0)
			{
				i++;
				box_check_phi = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_check_theta_unc") == 0)
			{
				i++;
				box_check_theta_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_check_phi_unc") == 0)
			{
				i++;
				box_check_phi_unc = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-box_check_overall_theta") == 0)
			{
				i++;
				box_check_overall_theta = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-phot_check_n") == 0)
			{
				i++;
				phot_check_n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-phot_check_max_theta") == 0)
			{
				i++;
				phot_check_max_theta = atof(argv[i]);
			}
			else
			{
				printf("Unrecognized argument: %s\n",argv[i]);
			}
		}
	}


	double main_mirror_angle = 74.11+mirror_angle_change;



	// 	double t_unc = .05;


	double rad_to_deg = 57.2958;

	double res_enhance = 1;







	double prog_thresh = 500;

	if (flatten_time == true)
	{
		s_func_t = 100000000;
	}	

	double outcsv_x,outcsv_y,outcsv_t;
	outcsv_x = 0*35;//bars are 35mm wide
	outcsv_y = 0;//mm
	outcsv_t = 0;//ns


	if(out_csv){
		n_phi_phots = 4000;
		n_z_phots = 8;
		particle_y += outcsv_y;}
	double pdf_unc_red_fac = 1;


	TRandom3 spread_ang(rseed+3);


	DircOpticalSim *dirc_model = new DircOpticalSim(\
			rseed,\
			-1200 + mirror_r_difference,\
			300.38,\
			main_mirror_angle,\
			600,\
			47.87 + box_rot + mirror_angle_change);

	dirc_model->set_store_traveled(false);// uses LOTS of memory if set to true.
	dirc_model->set_wedge_mirror_rand(wedge_non_uniformity);
	dirc_model->set_three_seg_mirror(three_seg_mirror);
	dirc_model->set_kaleidoscope_plot(kaleidoscope_plot);	
	dirc_model->set_pmt_plane_zs(pmt_min_z,pmt_max_z);
	dirc_model->set_large_mirror_zs(large_mirror_min_z,large_mirror_max_z);


	double muon_beta, pion_beta, kaon_beta/*, electron_beta:=1*/;
	muon_beta=pion_beta=kaon_beta=-1;
	double muon_angle, pion_angle, kaon_angle;
	muon_angle=pion_angle=kaon_angle = -1;


	TFile* tfile = new TFile(rootfilename,"RECREATE");

	TH1F *ll_diff_pion = new TH1F("ll_diff_pion","Difference of log likelihood real = pion",200000,-200,200);
	TH1F *ll_diff_kaon = new TH1F("ll_diff_kaon","Difference of log likelihood real = kaon",200000,-200,200);
	TH1F *ll_diff_pion_up = new TH1F("ll_diff_pion_up","Difference of log likelihood real = pion \"up\" going photons",200000,-200,200);
	TH1F *ll_diff_kaon_up = new TH1F("ll_diff_kaon_up","Difference of log likelihood real = kaon \"up\" going photons",200000,-200,200);
	TH1F *ll_diff_pion_down = new TH1F("ll_diff_pion_down","Difference of log likelihood real = pion \"down\" going photons",200000,-200,200);
	TH1F *ll_diff_kaon_down = new TH1F("ll_diff_kaon_down","Difference of log likelihood real = kaon \"down\" going photons",200000,-200,200);
	TH1F *phot_found_pion = new TH1F("phot_found_pion","number of photons found on pion angle", 1001,-.5,1000.5);
	TH1F *phot_found_kaon = new TH1F("phot_found_kaon","number of photons found on kaon angle", 1001,-.5,1000.5);
	TH1F *phot_found_pion_up = new TH1F("phot_found_pion_up","number of photons found on pion angle \"up\" going photons", 1001,-.5,1000.5);
	TH1F *phot_found_kaon_up = new TH1F("phot_found_kaon_up","number of photons found on kaon angle \"up\" going photons", 1001,-.5,1000.5);
	TH1F *phot_found_pion_down = new TH1F("phot_found_pion_down","number of photons found on pion angle \"down\" going photons", 1001,-.5,1000.5);
	TH1F *phot_found_kaon_down = new TH1F("phot_found_kaon_down","number of photons found on kaon angle \"down\" going photons", 1001,-.5,1000.5);

	TH1F *ref_pion_before = new TH1F("ref_pion_before","Angle of Pion Photons going into interface", 9000,0,90);
	TH1F *ref_kaon_before = new TH1F("ref_kaon_before","Angle of Kaon Photons going into interface", 9000,0,90);
	TH1F *ref_pion_after = new TH1F("ref_pion_after","Angle of Pion Photons going out of interface", 9000,0,90);
	TH1F *ref_kaon_after = new TH1F("ref_kaon_after","Angle of Kaon Photons going out of interface", 9000,0,90);
	TH1F *ref_pion_sens_plane = new TH1F("ref_pion_sens_plane","Angle of Pion Photons going into window to PMTs", 9000,0,90);
	TH1F *ref_kaon_sens_plane = new TH1F("ref_kaon_sens_plane","Angle of Kaon Photons going into window to PMTs", 9000,0,90);



	TH1F *pion_dist_x = new TH1F("pion_dist_x","x val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx);
	TH1F *pion_dist_y = new TH1F("pion_dist_y","y val of intercepted points - pion",(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH1F *pion_dist_t = new TH1F("pion_dist_t","t val of intercepted points - pion",(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH1F *kaon_dist_x = new TH1F("kaon_dist_x","x val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx);
	TH1F *kaon_dist_y = new TH1F("kaon_dist_y","y val of intercepted points - kaon",(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH1F *kaon_dist_t = new TH1F("kaon_dist_t","t val of intercepted points - kaon",(maxt-mint)/(res_enhance*rest),mint,maxt);

	TH1F *box_check_x = new TH1F("box_check_x","Check Box X",(maxx-minx)/(.1*res_enhance*resx),minx,maxx);
	TH1F *box_check_y = new TH1F("box_check_y","Check Box Y",(maxy-miny)/(.1*res_enhance*resy),miny,maxy);
	TH1F *box_check_d = new TH1F("box_check_d","Check Box Dist",800000,0,8000);

	TH1F *pion_dist_geant_x = new TH1F("pion_dist_geant_x","x val of intercepted points (from geant) - pion",(maxx-minx)/(res_enhance*resx),minx,maxx);
	TH1F *pion_dist_geant_y = new TH1F("pion_dist_geant_y","y val of intercepted points (from geant) - pion",(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH1F *pion_dist_geant_t = new TH1F("pion_dist_geant_t","t val of intercepted points (from geant) - pion",(maxt-mint)/(res_enhance*rest),mint,maxt);

	TH2F *ref_theta_cphi_pion = new TH2F("ref_theta_cphi_pion","Emitted angle of Pion Photons versus angle into interface", 7200,0,360,1800,0,90);
	TH2F *ref_theta_cphi_kaon = new TH2F("ref_theta_cphi_kaon","Emitted angle of Kaon Photons versus angle into interface", 7200,0,360,1800,0,90);

	TH2F *pion_dist_xy = new TH2F("pion_dist_xy","xy val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH2F *kaon_dist_xy = new TH2F("kaon_dist_xy","xy val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

	TH2F *box_check_xy = new TH2F("box_check_xy","xy val of points from top of wedge",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

	TH2F *pion_dist_geant_xy = new TH2F("pion_dist_geant_xy","xy val of intercepted points (from geant) - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

	TH2F *pion_dist_xt = new TH2F("pion_dist_xt","xt val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *kaon_dist_xt = new TH2F("kaon_dist_xt","xt val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);

	TH2F *pion_dist_geant_xt = new TH2F("pion_dist_geant_xt","xt val of intercepted points (from geant) - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);

	TH2F *pion_dist_yt = new TH2F("pion_dist_yt","yt val of intercepted points - pion",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);
	TH2F *kaon_dist_yt = new TH2F("kaon_dist_yt","yt val of intercepted points - kaon",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);

	TH2F *pion_dist_geant_yt = new TH2F("pion_dist_geant_yt","yt val of intercepted points (from geant) - pion",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);


	maxy *= 5;

	TH2F *pion_coverage_xy = new TH2F("pion_coverage_xy","xy val of generated points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
	TH2F *kaon_coverage_xy = new TH2F("kaon_coverage_xy","xy val of generated points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

	TH2F *kinematics_yields = new TH2F("kinematics_yields","Fractional Yields versus angles",(kinematics_yields_max_theta-kinematics_yields_min_theta)/kinematics_yields_step_theta*1.00001,\
			kinematics_yields_min_theta,\
			kinematics_yields_max_theta,\
			(kinematics_yields_max_phi-kinematics_yields_min_phi)/kinematics_yields_step_phi*1.00001,\
			kinematics_yields_min_phi,\
			kinematics_yields_max_phi);

	TH1F *pion_cerenkov = new TH1F("pion_cerenkov","pion cerenkov angle distribution",300,45,48);
	TH1F *kaon_cerenkov = new TH1F("kaon_cerenkov","kaon cerenkov angle distribution",300,45,48);
	TH1F *pion_lambda = new TH1F("pion_lambda","pion wavelength distribution",450,250,700);
	TH1F *kaon_lambda = new TH1F("kaon_lambda","kaon wavelength distribution",450,250,700);

	TH1F *liquid_dist = new TH1F("liquid_dist","distance travelled in liquid (mm)",1500,0,1500);

	TH1F *simulation_time = new TH1F("simulation_time","Simulation time per particle",4001,-.5,4000.5);

	TH1F *zero_pion_id = new TH1F("zero_pion_id","Correct pion ID with cut at 0",2,-0.5,1.5);
	TH1F *zero_kaon_id = new TH1F("zero_kaon_id","Correct kaon ID with cut at 0",2,-0.5,1.5);


	DircDigitizer digitizer(\
			minx,\
			maxx,\
			resx,\
			digit_miny,\
			digit_maxy,\
			resy,\
			t_unc,\
			t_bin_size);

	//  	
	printf("Beginning Run\n");
	double llc, llf, ll_diff;
	llc=llf=ll_diff=0;
	int pion_kaon_diff = 0;
	std::vector<dirc_point> sim_points;
	std::vector<dirc_point> confound_points;
	dirc_model->set_focmirror_nonuniformity(main_mirror_nonuniformity);
	/*   read in file for events and 
	 *   initialize an array 
	 *    */
	if (inputrootfile)
	{
		printf("Opening root file: %s\n",inputrootfilename);
		//Gross... I got root types in my loops...
		Float_t px[10];
		Float_t py[10];
		Float_t pz[10];
		Float_t particle_x[10];
		Float_t particle_y[10];
		Float_t particle_z[10];
		Float_t gen_t[10];
		Float_t gen_E[10];
		Float_t hitw[2000];
		Float_t hity[2000];
		Float_t hit_t[2000];
		int nHit;
		Int_t particle_id[10];

		TFile* inputrootfile = new TFile(inputrootfilename,"READONLY");
		TTree *in_tree = (TTree*) inputrootfile->Get("dircTreePID");

		int nentries = in_tree->GetEntries();

		in_tree->SetBranchAddress("ParticleID",&particle_id);
		in_tree->SetBranchAddress("genTrk_px",&px);
		in_tree->SetBranchAddress("genTrk_py",&py);
		in_tree->SetBranchAddress("genTrk_pz",&pz);
		in_tree->SetBranchAddress("genTrk_x",&particle_x);
		in_tree->SetBranchAddress("genTrk_y",&particle_y);
		in_tree->SetBranchAddress("genTrk_z",&particle_z);
		in_tree->SetBranchAddress("genTrk_t",&gen_t);
		in_tree->SetBranchAddress("genTrk_E",&gen_E);

		in_tree->SetBranchAddress("hitpixel_w",&hitw);
		in_tree->SetBranchAddress("hitpixel_y",&hity);
		in_tree->SetBranchAddress("hit_t",&hit_t);

		in_tree->SetBranchAddress("nHit",&nHit);

		std::vector<std::vector<dirc_point> > root_hits_array;

		unsigned int r = 0;

		int iPID, iBAR;
		//int ievent_index;
		double ix,iy,it,itheta,iphi,iE;
		std::vector<int> PID,BAR;
		std::vector<int> event_index;
		std::vector<double> x,y,t,theta,phi,E;

		double min_t = 1000;
		double avg_x = 0;
		double avg_y = 0;
		double avg_t = 0;
		int tot_nhit = 0;
		for (int i = 0; i < nentries; i++) {

			if (r >= (unsigned int) max_particles) break;

			in_tree->GetEntry(i);
			// 			if (iE > 4) continue;	
			//TODO read more than 1 particle per event
			iPID = particle_id[0];
			iBAR = dirc_model->get_bar_from_x(particle_x[0]*10);

			if (iBAR != 6) {continue;}
			//			if (iPID == 8) {continue;}
			ix = particle_x[0]*10 - dirc_model->get_bar_offset(iBAR);
			iy = particle_y[0]*10 + 700;//Figure out that geometry
			it = gen_t[0] + 2*r;//2ns spacing of events at least.
			iE = gen_E[0];
			double p_trans = sqrt(py[0]*py[0]+px[0]*px[0]);
			itheta = atan(p_trans/pz[0])*180/3.14159;
			iphi = atan2(py[0],px[0])*180/3.14159;

			if (force_kinematics == true)
			{
				itheta = particle_theta;
				iphi = particle_phi;
				iE = energy;
			}
			event_index.push_back(i);
			PID.push_back(iPID);
			BAR.push_back(iBAR);
			x.push_back(ix);
			y.push_back(iy);
			t.push_back(it);
			theta.push_back(itheta);
			phi.push_back(iphi);
			E.push_back(iE);

			if (t[r] < 0) printf("low t: %12.04f\n",t[r]);

			std::vector<dirc_point> root_hits_add;
			if (iPID == 8 || iPID == 9){
				tot_nhit += nHit;
			}
			for (int j = 0; j < nHit; j++)
			{
				dirc_point addpt;
				//these are swapped knowingly
				addpt.x = hity[j]*10 + 151.4 + 22;
				addpt.y = hitw[j]*10 - 28.85 + 25.1;
				//Strange extra hits at 149.15 for some reason
				if (addpt.y > 149.1 && addpt.y < 149.2) continue;
				//addpt.t = t[r] + hit_t[j] - 11;
				addpt.t = hit_t[j] + 2*r;

				if (addpt.t - 2*r < 0)
				{//hmmmm can't be - ignore it
					printf("low addt: %12.04f %d %d\n", addpt.t - 2*r, i, j);
					if (iPID == 8 || iPID == 9)
					{
						tot_nhit--;
					}
					continue;
				}

				root_hits_add.push_back(addpt);



				if (fill_distributions == true)
				{
					addpt.t -= 2*r;
					pion_dist_geant_x->Fill(addpt.x);
					pion_dist_geant_y->Fill(addpt.y);
					pion_dist_geant_t->Fill(addpt.t);
					pion_dist_geant_xy->Fill(addpt.x,addpt.y);
					pion_dist_geant_xt->Fill(addpt.x,addpt.t);
					pion_dist_geant_yt->Fill(addpt.y,addpt.t);
					addpt.t += 2*r;
				}

				if (addpt.t < 0) printf("low addpt: %12.04f\n",addpt.t);
				if ((iPID == 8 || iPID == 9))
				{
					min_t = std::min(min_t,addpt.t - 2*r);
					avg_x += addpt.x;
					avg_y += addpt.y;
					avg_t += addpt.t - 2*r;
				}
			}
			root_hits_array.push_back(root_hits_add);

			r++;
		}
		double high_avg = 0;
		double low_avg = 0;
		double high_min = 10000000;
		double low_min = 10000000;
		int high_count = 0;
		int low_count = 0;
		double tmp_t = -1;
		for (unsigned int j = 0; j < root_hits_array.size(); j++)
		{
			for (unsigned int k = 0; k < root_hits_array[j].size(); k++)
			{
				tmp_t = root_hits_array[j][k].t - 2*j;
				if (tmp_t < 0) printf("low t: %12.04f %d %d\n",tmp_t,j,k);
				if (tmp_t > avg_t/tot_nhit)
				{
					high_count++;
					high_avg += tmp_t;
					high_min = std::min(high_min, tmp_t);
				}
				else
				{
					low_count++;
					low_avg += tmp_t;
					low_min = std::min(low_min, tmp_t);
				}
			}
		}


		printf("Avg hit x y t: %12.04f %12.04f %12.04f\n",avg_x/tot_nhit,avg_y/tot_nhit, avg_t/tot_nhit);
		printf("min input t delta: %12.04f\n",min_t);
		printf("T abv avg : below avg : diff : %12.04f %12.04f %12.04f\n", high_avg/high_count, low_avg/low_count, high_avg/high_count - low_avg/low_count);
		printf("T abv min : below min : diff : %12.04f %12.04f %12.04f\n", high_min, low_min, high_min - low_min);
		inputrootfile->Close("R");

		int mc_tally=0;
//		int confounded_tally=0;
		num_runs=0;
		coverage_plot=false;

		//declare memory addresses for inputs
		//units: mm, ns, deg, GeV 

		//declare for beta and angle based on PID, and constructors for hitpoints and pdfs
		double pion_mc_beta,kaon_mc_beta,pion_mc_angle,kaon_mc_angle;
		std::vector<dirc_point> hits_trk_is_pion;
		std::vector<dirc_point> hits_trk_is_kaon;

		double min_pval = 1e-2;
		DircSpreadGaussian * pdf_as_pion = new DircSpreadGaussian(\
				sfunc_sig,\
				hits_trk_is_pion,\
				s_func_x,\
				s_func_y,\
				s_func_t,\
				min_pval);
		DircSpreadGaussian * pdf_as_kaon = new DircSpreadGaussian(\
				sfunc_sig,\
				hits_trk_is_kaon,\
				s_func_x,\
				s_func_y,\
				s_func_t,\
				min_pval);//could make this an array that is filled...
		//Might be worth making a new copy each time
		DircProbabilitySeparation * sep_pdfs_mc;

		DircProgressiveSeparation *progressive_separation = \
								    new DircProgressiveSeparation(\
										    dirc_model,\
										    n_phi_phots*n_z_phots,\
										    n_step_phots,\
										    sfunc_sig,\
										    s_func_x,\
										    s_func_y,\
										    s_func_t,\
										    kmass,\
										    pimass,\
										    prog_thresh);

		//		r = 100;

		for(unsigned int n=0;n < r; n++){

			min_t = 1000;

			if(fabs(PID[n])>7 && fabs(PID[n])<13)//8 and 9 are pions, 11 and 12 are kaons
			{

				mc_tally++;

				//
				//and replace acos......
				//cdd replace 1.47 with double refrac_index=1.47;
				pion_mc_beta = dirc_model->get_beta(E[n],pimass);
				kaon_mc_beta = dirc_model->get_beta(E[n],kmass);
				pion_mc_angle = rad_to_deg*acos(1/(1.47*pion_mc_beta));
				kaon_mc_angle = rad_to_deg*acos(1/(1.47*kaon_mc_beta));

				//sep_pdfs_mc = new DircProbabilitySeparation(pdf_as_kaon,pdf_as_pion);

				dirc_model->set_focus_mirror_angle(\
						spread_ang.Gaus(main_mirror_angle,mirror_angle_change_unc),\
						spread_ang.Gaus(0,mirror_angle_change_yunc));
				dirc_model->set_upper_wedge_angle_diff(\
						spread_ang.Gaus(0,wedge_uncertainty),\
						spread_ang.Gaus(0,upper_wedge_yang_spread));
				printf("\r                                                                                                      ");
				printf("\rRunning monte carlo event  %i/%i. Found %i pions/kaons", n,r-1,mc_tally);

				fflush(stdout);

				//simulate pmt hits from primary mc event to compare with preconstructed pdfs
				//


				if(abs(PID[n])==8 || PID[n]==9 ){
					dirc_model->sim_rand_n_photons(\
							sim_points,\
							n_sim_phots,\
							pion_mc_angle,\
							BAR[n],\
							x[n],\
							y[n],\
							0,\
							theta[n],\
							phi[n],\
							tracking_unc,\
							ckov_unc,\
							pion_mc_beta);
				//TODO: Check timing
					pion_kaon_diff = 1;
				}
				else{
					dirc_model->sim_rand_n_photons(\
							sim_points,\
							n_sim_phots,\
							kaon_mc_angle,\
							BAR[n],\
							x[n],\
							y[n],\
							0,\
							theta[n],\
							phi[n],\
							tracking_unc,\
							ckov_unc,\
							kaon_mc_beta);
				//TODO: Check timing
					pion_kaon_diff = -1;
				}


				for (unsigned int i = 0; i < sim_points.size(); i++)
				{
					sim_points[i].t += t[n];
					// 					printf("%12.04f\n",sim_points[i].t);
				}
				//simulate pmt hits from background events that occur within t ns of the primary event

				/*Ignore confounding for now
				  confounded_tally=0;
				  for (unsigned int j =0; j < r;j++){

				  if (abs(event_index[n] - event_index[j]) > broaden_events)
				  {
				  continue;
				  }

				  if(fabs(t[n]-t[j])<time_window && j!=n){
				  confounded_tally++; 
				  if(abs(PID[j])==2 ||PID[j]==3){

				  dirc_model->sim_rand_n_photons(\
				  confound_points,\
				  n_sim_phots,\
				  47.135,\
				  BAR[j],\
				  x[j],\
				  y[j],\
				  theta[j],\
				  phi[j],\
				  tracking_unc,\
				  ckov_unc,\
				  1);
				  }
				  else if(abs(PID[j])==8||PID[j]==9){
				  pion_beta = dirc_model->get_beta(E[j],pimass);
				  pion_angle = rad_to_deg*acos(1/(refrac_index*pion_beta));

				  dirc_model->sim_rand_n_photons(\
				  confound_points,\
				  n_sim_phots,\
				  pion_angle,\
				  BAR[j],\
				  x[j],\
				  y[j],\
				  theta[j],\
				  phi[j],\
				  tracking_unc,\
				  ckov_unc,\
				  pion_beta);

				  }
				  else if(abs(PID[j])==11 || PID[j]==12){
				  kaon_beta = dirc_model->get_beta(E[j],kmass);
				  kaon_angle = rad_to_deg*acos(1/(refrac_index*kaon_beta));

				  dirc_model->sim_rand_n_photons(\
				  confound_points,\
				  n_sim_phots,\
				  kaon_angle,\
				  BAR[j],\
				  x[j],\
				  y[j],\
				  theta[j],\
				  phi[j],\
				  tracking_unc,\
				  ckov_unc,\
				  kaon_beta);

				  }
				  else if(abs(PID[j])==5 || PID[j]==6){
				  muon_beta = dirc_model->get_beta(E[j],mumass);
				  muon_angle = rad_to_deg*acos(1/(refrac_index*muon_beta));

				  dirc_model->sim_rand_n_photons(\
				  confound_points,\
				  n_sim_phots,\
				  muon_angle,\
				  BAR[j],\
				x[j],\
					y[j],\
					theta[j],\
					phi[j],\
					tracking_unc,\
					ckov_unc,\
					muon_beta);
			}
			}
			for (unsigned int i = 0; i < confound_points.size(); i++)
			{
				confound_points[i].t += t[j];
				// 						confound_points[i].x += 150;
				sim_points.push_back(confound_points[i]);
			}
			confound_points.clear();
			}//end confounded point generation
			*/
				// 				printf("Found %i confounding events for track %i. sim_points.size()=%lu\n",confounded_tally,n,sim_points.size());


				//Begin analysis

				//Zero out uncertainties (assume we are starting from calibrated versions

				//				digitizer.digitize_points(sim_points);

				//Ignore all the simpoints - use the root file
				sim_points = root_hits_array[n];

			timing_clock = clock();

			if (use_prog_sep == true)
			{
				dirc_model->set_focus_mirror_angle(\
						spread_ang.Gaus(main_mirror_angle,0),\
						spread_ang.Gaus(0,0));
				dirc_model->set_upper_wedge_angle_diff(\
						spread_ang.Gaus(0,0),\
						spread_ang.Gaus(0,0));
				//Make Prog Sep with proper points
				/*
				   ll_diff = progressive_separation->get_ll_progressive(\
				   sim_points,\
				   BAR[n],\
				   E[n],\
				   x[n],\
				   y[n],\
				   theta[n],\
				   phi[n],\
				   tracking_unc,\
				   ckov_unc);
				 */
			}
			else
			{
				dirc_model->set_focus_mirror_angle(\
						spread_ang.Gaus(main_mirror_angle,0),\
						spread_ang.Gaus(0,0));
				dirc_model->set_upper_wedge_angle_diff(\
						spread_ang.Gaus(0,0),\
						spread_ang.Gaus(0,0));

				dirc_model->sim_reg_n_photons(\
						hits_trk_is_pion,\
						n_phi_phots,\
						n_z_phots,\
						pion_mc_angle,\
						BAR[n],\
						x[n],\
						y[n],\
						0,
						theta[n],\
						phi[n],\
						0,\
						ckov_unc/pdf_unc_red_fac,\
						pion_mc_beta);
				//TODO: Check timing

				dirc_model->sim_reg_n_photons(\
						hits_trk_is_kaon,\
						n_phi_phots,\
						n_z_phots,\
						kaon_mc_angle,
						BAR[n],\
						x[n],\
						y[n],\
						0,
						theta[n],\
						phi[n],\
						0,\
						ckov_unc/pdf_unc_red_fac,\
						kaon_mc_beta);
				//TODO: Check timing
				double avg_x = 0;	
				double avg_y = 0;	
				double avg_t = 0;	
				for (unsigned int k = 0; k < hits_trk_is_pion.size(); k++)
				{
					hits_trk_is_pion[k].t += t[n];
					min_t = std::min(min_t,hits_trk_is_pion[k].t - 2*n);
					avg_x += hits_trk_is_pion[k].x;
					avg_y += hits_trk_is_pion[k].y;
					avg_t += hits_trk_is_pion[k].t - 2*n;
				}

				double high_avg = 0;
				double low_avg = 0;
				double high_min = 10000000;
				double low_min = 10000000;
				int high_count = 0;
				int low_count = 0;
				double tmp_t = -1;
				double tmp_x = -1;
				double tmp_y = -1;
				for (unsigned int k = 0; k < hits_trk_is_pion.size(); k++)
				{
					tmp_t = hits_trk_is_pion[k].t - 2*n;
					tmp_x = hits_trk_is_pion[k].x;
					tmp_y = hits_trk_is_pion[k].y;
					if (n == 0 && fill_distributions == true)
					{
						pion_dist_x->Fill(tmp_x);
						pion_dist_y->Fill(tmp_y);
						pion_dist_t->Fill(tmp_t);
						pion_dist_xy->Fill(tmp_x, tmp_y);
						pion_dist_xt->Fill(tmp_x, tmp_t);
						pion_dist_yt->Fill(tmp_t, tmp_t);
					}

					if (tmp_t > avg_t/hits_trk_is_pion.size())
					{
						high_count++;
						high_avg += tmp_t;
						high_min = std::min(high_min,tmp_t);
					}
					else
					{
						low_count++;
						low_avg += tmp_t;
						low_min = std::min(low_min,tmp_t);
					}
				}
				printf("\nAvg t High : low : diff : %12.04f : %12.04f : % 12.04f\n",high_avg/high_count,low_avg/low_count,high_avg/high_count-low_avg/low_count);
				printf("Min t High : low : diff : %12.04f : %12.04f : % 12.04f\n",high_min,low_min,high_min - low_min);

				for (unsigned int k = 0; k < hits_trk_is_kaon.size(); k++)
				{
					hits_trk_is_kaon[k].t += t[n];
				}
				for (unsigned int k = 0; k < hits_trk_is_kaon.size(); k++)
				{
					tmp_t = hits_trk_is_kaon[k].t - 2*n;
					tmp_x = hits_trk_is_kaon[k].x;
					tmp_y = hits_trk_is_kaon[k].y;
					if (n == 0 && fill_distributions == true)
					{
						kaon_dist_x->Fill(tmp_x);
						kaon_dist_y->Fill(tmp_y);
						kaon_dist_t->Fill(tmp_t);
						kaon_dist_xy->Fill(tmp_x, tmp_y);
						kaon_dist_xt->Fill(tmp_x, tmp_t);
						kaon_dist_yt->Fill(tmp_t, tmp_t);
					}

					if (tmp_t > avg_t/hits_trk_is_kaon.size())
					{
						high_count++;
						high_avg += tmp_t;
						high_min = std::min(high_min,tmp_t);
					}
					else
					{
						low_count++;
						low_avg += tmp_t;
						low_min = std::min(low_min,tmp_t);
					}
				}
				//printf("Pion Avg x y t: %12.04f %12.04f %12.04f\n",avg_x/hits_trk_is_kaon.size(),avg_y/hits_trk_is_kaon.size(),avg_t/hits_trk_is_kaon.size());
				//printf("kaon min_t: %12.04f\n",min_t);
				pdf_as_pion->set_support(hits_trk_is_pion);
				pdf_as_kaon->set_support(hits_trk_is_kaon);


				// 					printf("betas: %12.04f %12.04f\n",pion_mc_beta,kaon_mc_beta);

				llc = pdf_as_pion->get_log_likelihood(sim_points);
				llf = pdf_as_kaon->get_log_likelihood(sim_points);

				// 					printf("%06d\n",sim_points.size());

				ll_diff = llc - llf;

				// 					printf("\n%12.04f %12.04f %12.04f\n",llc,llf,ll_diff);
				int correct_id = 1;
				if (ll_diff*pion_kaon_diff < 0)
				{
					correct_id = 0;
					for (unsigned int k = 0; k < sim_points.size(); k++)
					{
						//  							printf("%12.04f %12.04f %12.04f ll: %12.04f %12.04f, %12.04f\n",sim_points[k].x,sim_points[k].y,sim_points[k].t,\
						log(pdf_as_pion->get_single_likelihood(sim_points[k])),\
							log(pdf_as_kaon->get_single_likelihood(sim_points[k])),\
							log(pdf_as_pion->get_single_likelihood(sim_points[k]))-log(pdf_as_kaon->get_single_likelihood(sim_points[k])));
					}
				}
				for (unsigned int k = 0; k < sim_points.size(); k++)
				{
					//  						printf("%12.04f %12.04f %12.04f\n",sim_points[k].x,sim_points[k].y,sim_points[k].t);
				}
				if (ll_diff == ll_diff)
				{
					// 						printf("\n%12.04f %12.04f %12.04f %d\n",llc,llf,ll_diff,correct_id);
					for (unsigned int k = 0; k < sim_points.size(); k++)
					{
						// 							printf("%12.04f %12.04f %12.04f ll: %12.04f %12.04f, %12.04f\n",sim_points[k].x,sim_points[k].y,sim_points[k].t,\
						log(pdf_as_pion->get_single_likelihood(sim_points[k])),\
							log(pdf_as_kaon->get_single_likelihood(sim_points[k])),\
							log(pdf_as_pion->get_single_likelihood(sim_points[k]))-log(pdf_as_kaon->get_single_likelihood(sim_points[k])));

					}
				}
			}

			timing_clock = clock() - timing_clock;

			simulation_time->Fill(((float)timing_clock)/(CLOCKS_PER_SEC)*1000);


			// 				printf("\nPID[n]=%i loglikehood difference(pion-kaon)=%f\n",PID[n],ll_diff);
			if(abs(PID[n])==8||PID[n]==9){
				ll_diff_pion->Fill(ll_diff);

				if (ll_diff > 0)
				{
					zero_pion_id->Fill(1);
				}
				else
				{
					zero_pion_id->Fill(0);
				}

				phot_found_pion->Fill(sim_points.size());
			}
			else{

				ll_diff_kaon->Fill(ll_diff);
				if (ll_diff < 0)
				{
					zero_kaon_id->Fill(1);
				}
				else
				{
					zero_kaon_id->Fill(0);
				}
				phot_found_kaon->Fill(sim_points.size());
			}
			}
		}
		printf("\n");//ensure we're on a new line at the end
		printf("Kaon missID: %12.04f%% \n",(double) 100.0*zero_kaon_id->GetBinContent(1)/(zero_kaon_id->GetBinContent(1)+zero_kaon_id->GetBinContent(2)));
		printf("Pion missID: %12.04f%% \n",(double) 100.0*zero_pion_id->GetBinContent(1)/(zero_pion_id->GetBinContent(1)+zero_pion_id->GetBinContent(2)));
	}
	else if(inputfile){

		printf("Csv file with name %s selected.\n",in_str);
		printf("Opening %s with time window of %fns\n",in_str,time_window);

		int mc_tally=0;
		int confounded_tally=0;
		num_runs=0;
		coverage_plot=false;
		std::ifstream f(in_str);
		if(!f.is_open()){printf("%s not found!",in_str);return -1;}	

		int line_buf_size = 1000; 
		char  s[line_buf_size];
		f.getline(s,line_buf_size); //skips first line
		//declare memory addresses for inputs
		//units: mm, ns, deg, GeV 
		int iPID, iBAR,ievent_index;
		double ix,iy,it,itheta,iphi,iE;
		std::vector<int> PID,BAR,event_index;
		std::vector<double> x,y,t,theta,phi,E;

		//declare for beta and angle based on PID, and constructors for hitpoints and pdfs
		double pion_mc_beta,kaon_mc_beta,pion_mc_angle,kaon_mc_angle;
		std::vector<dirc_point> hits_trk_is_pion;
		std::vector<dirc_point> hits_trk_is_kaon;

		double min_pval = 1e-2;
		DircSpreadGaussian * pdf_as_pion = new DircSpreadGaussian(\
				sfunc_sig,\
				hits_trk_is_pion,\
				s_func_x,\
				s_func_y,\
				s_func_t,\
				min_pval);
		DircSpreadGaussian * pdf_as_kaon = new DircSpreadGaussian(\
				sfunc_sig,\
				hits_trk_is_kaon,\
				s_func_x,\
				s_func_y,\
				s_func_t,\
				min_pval);//could make this an array that is filled...
		//Might be worth making a new copy each time
		DircProbabilitySeparation * sep_pdfs_mc;

		DircProgressiveSeparation *progressive_separation = \
								    new DircProgressiveSeparation(\
										    dirc_model,\
										    n_phi_phots*n_z_phots,\
										    n_step_phots,\
										    sfunc_sig,\
										    s_func_x,\
										    s_func_y,\
										    s_func_t,\
										    kmass,\
										    pimass,\
										    prog_thresh);


		unsigned int r=0;
		while(f>>ievent_index>>iPID>>iBAR>>ix>>iy>>it>>itheta>>iphi>>iE)
		{
			// 			if (iE > 4) continue;
			if (r >= (unsigned int) max_particles) break;

			if (force_kinematics == true)
			{
				itheta = particle_theta;
				iphi = particle_phi;
				iE = energy;
			}

			event_index.push_back(ievent_index);
			PID.push_back(iPID);
			BAR.push_back(iBAR);
			x.push_back(ix);
			y.push_back(iy);
			t.push_back(it);
			theta.push_back(itheta);
			phi.push_back(iphi);
			E.push_back(iE);
			r++;
		}
		f.close();
		for(unsigned int n=0;n<r;n++){
			if(fabs(PID[n])>7 && fabs(PID[n])<13)//8 and 9 are pions, 11 and 12 are kaons
			{

				mc_tally++;

				//
				//and replace acos......
				//cdd replace 1.47 with double refrac_index=1.47;
				pion_mc_beta = dirc_model->get_beta(E[n],pimass);
				kaon_mc_beta = dirc_model->get_beta(E[n],kmass);
				pion_mc_angle = rad_to_deg*acos(1/(1.47*pion_mc_beta));
				kaon_mc_angle = rad_to_deg*acos(1/(1.47*kaon_mc_beta));

				//sep_pdfs_mc = new DircProbabilitySeparation(pdf_as_kaon,pdf_as_pion);

				dirc_model->set_focus_mirror_angle(\
						spread_ang.Gaus(main_mirror_angle,mirror_angle_change_unc),\
						spread_ang.Gaus(0,mirror_angle_change_yunc));
				dirc_model->set_upper_wedge_angle_diff(\
						spread_ang.Gaus(0,wedge_uncertainty),\
						spread_ang.Gaus(0,upper_wedge_yang_spread));
				printf("\r                                                                                                      ");
				printf("\rRunning monte carlo event  %i/%i. Found %i pions/kaons", n,r-1,mc_tally);

				fflush(stdout);

				//simulate pmt hits from primary mc event to compare with preconstructed pdfs
				//


				if(abs(PID[n])==8 || PID[n]==9 ){
					dirc_model->sim_rand_n_photons(\
							sim_points,\
							n_sim_phots,\
							pion_mc_angle,\
							BAR[n],\
							x[n],\
							y[n],\
							0,\
							theta[n],\
							phi[n],\
							0,\
							ckov_unc,\
							pion_mc_beta);
				//TODO: Check timing
					pion_kaon_diff = 1;
				}
				else{
					dirc_model->sim_rand_n_photons(\
							sim_points,\
							n_sim_phots,\
							kaon_mc_angle,\
							BAR[n],\
							x[n],\
							y[n],\
							0,\
							theta[n],\
							phi[n],\
							0,\
							ckov_unc,\
							kaon_mc_beta);
				//TODO: Check timing
					pion_kaon_diff = -1;
				}
				for (unsigned int i = 0; i < sim_points.size(); i++)
				{
					sim_points[i].t += t[n];
					// 					sim_points[i].x += 150;
					// 					printf("%12.04f\n",sim_points[i].t);
				}
				//simulate pmt hits from background events that occur within t ns of the primary event
				confounded_tally=0;
				for (unsigned int j =0; j < r;j++){

					if (abs(event_index[n] - event_index[j]) > broaden_events)
					{
						continue;
					}

					if(fabs(t[n]-t[j])<time_window && j!=n){
						confounded_tally++; 
						if(abs(PID[j])==2 ||PID[j]==3){

							dirc_model->sim_rand_n_photons(\
									confound_points,\
									n_sim_phots,\
									47.135,\
									BAR[j],\
									x[j],\
									y[j],\
									0,\
									theta[j],\
									phi[j],\
									tracking_unc,\
									ckov_unc,\
									1);
						//TODO: Check timing
						}
						else if(abs(PID[j])==8||PID[j]==9){
							pion_beta = dirc_model->get_beta(E[j],pimass);
							pion_angle = rad_to_deg*acos(1/(refrac_index*pion_beta));

							// 							dirc_model->sim_rand_n_photons(\
							confound_points,\
								n_sim_phots,\
								pion_angle,\
								BAR[j],\
								x[j],\
								y[j],\
								theta[j],\
								phi[j],\
								tracking_unc,\
								ckov_unc,\
								pion_beta);

						}
						else if(abs(PID[j])==11 || PID[j]==12){
							kaon_beta = dirc_model->get_beta(E[j],kmass);
							kaon_angle = rad_to_deg*acos(1/(refrac_index*kaon_beta));

							// 							dirc_model->sim_rand_n_photons(\
							confound_points,\
								n_sim_phots,\
								kaon_angle,\
								BAR[j],\
								x[j],\
								y[j],\
								theta[j],\
								phi[j],\
								tracking_unc,\
								ckov_unc,\
								kaon_beta);

						}
						else if(abs(PID[j])==5 || PID[j]==6){
							muon_beta = dirc_model->get_beta(E[j],mumass);
							muon_angle = rad_to_deg*acos(1/(refrac_index*muon_beta));

							dirc_model->sim_rand_n_photons(\
									confound_points,\
									n_sim_phots,\
									muon_angle,\
									BAR[j],\
									x[j],\
									y[j],\
									0,\
									theta[j],\
									phi[j],\
									tracking_unc,\
									ckov_unc,\
									muon_beta);
							//TODO: Check timing
						}
					}
					for (unsigned int i = 0; i < confound_points.size(); i++)
					{
						confound_points[i].t += t[j];
						// 						confound_points[i].x += 150;
						sim_points.push_back(confound_points[i]);
					}
					confound_points.clear();
				}//end confounded point generation

				//Begin analysis

				digitizer.digitize_points(sim_points);



				timing_clock = clock();

				if (use_prog_sep == true)
				{
					dirc_model->set_focus_mirror_angle(\
							spread_ang.Gaus(main_mirror_angle,0),\
							spread_ang.Gaus(0,0));
					dirc_model->set_upper_wedge_angle_diff(\
							spread_ang.Gaus(0,0),\
							spread_ang.Gaus(0,0));

					ll_diff = progressive_separation->get_ll_progressive(\
							sim_points,\
							BAR[n],\
							E[n],\
							x[n],\
							y[n],\
							theta[n],\
							phi[n],\
							tracking_unc,\
							ckov_unc);
				}
				else
				{
					dirc_model->set_focus_mirror_angle(\
							spread_ang.Gaus(main_mirror_angle,0),\
							spread_ang.Gaus(0,0));
					dirc_model->set_upper_wedge_angle_diff(\
							spread_ang.Gaus(0,0),\
							spread_ang.Gaus(0,0));

					dirc_model->sim_reg_n_photons(\
							hits_trk_is_pion,\
							n_phi_phots,\
							n_z_phots,\
							pion_mc_angle,\
							BAR[n],\
							x[n],\
							y[n],\
							0,
							theta[n],\
							phi[n],\
							0,\
							ckov_unc/pdf_unc_red_fac,\
							pion_mc_beta);
						//TODO: Check timing

					dirc_model->sim_reg_n_photons(\
							hits_trk_is_kaon,\
							n_phi_phots,\
							n_z_phots,\
							kaon_mc_angle,
							BAR[n],\
							x[n],\
							y[n],\
							0,
							theta[n],\
							phi[n],\
							0,\
							ckov_unc/pdf_unc_red_fac,\
							kaon_mc_beta);
						//TODO: Check timing

					for (unsigned int k = 0; k < hits_trk_is_pion.size(); k++)
					{
						hits_trk_is_pion[k].t += t[n];
					}
					for (unsigned int k = 0; k < hits_trk_is_kaon.size(); k++)
					{
						hits_trk_is_kaon[k].t += t[n];
					}
					pdf_as_pion->set_support(hits_trk_is_pion);
					pdf_as_kaon->set_support(hits_trk_is_kaon);


					// 					printf("betas: %12.04f %12.04f\n",pion_mc_beta,kaon_mc_beta);

					llc = pdf_as_pion->get_log_likelihood(sim_points);
					llf = pdf_as_kaon->get_log_likelihood(sim_points);

					// 					printf("%06d\n",sim_points.size());

					ll_diff = llc - llf;

					// 					printf("\n%12.04f %12.04f %12.04f\n",llc,llf,ll_diff);
					int correct_id = 1;
					if (ll_diff*pion_kaon_diff < 0)
					{
						correct_id = 0;

						printf("wrong: %12.04f %d\n",ll_diff,pion_kaon_diff);
						for (unsigned int k = 0; k < sim_points.size(); k++)
						{
							printf("%12.04f %12.04f %12.04f ll: %12.04f %12.04f, %12.04f\n",sim_points[k].x,sim_points[k].y,sim_points[k].t,\
									log(pdf_as_pion->get_single_likelihood(sim_points[k])),\
									log(pdf_as_kaon->get_single_likelihood(sim_points[k])),\
									log(pdf_as_pion->get_single_likelihood(sim_points[k]))-log(pdf_as_kaon->get_single_likelihood(sim_points[k])));
						}
					}
					for (unsigned int k = 0; k < sim_points.size(); k++)
					{
						//  						printf("%12.04f %12.04f %12.04f\n",sim_points[k].x,sim_points[k].y,sim_points[k].t);
					}
					if (ll_diff == ll_diff)
					{
						// 						printf("\n%12.04f %12.04f %12.04f %d\n",llc,llf,ll_diff,correct_id);
						for (unsigned int k = 0; k < sim_points.size(); k++)
						{
							// 							printf("%12.04f %12.04f %12.04f ll: %12.04f %12.04f, %12.04f\n",sim_points[k].x,sim_points[k].y,sim_points[k].t,\
							log(pdf_as_pion->get_single_likelihood(sim_points[k])),\
								log(pdf_as_kaon->get_single_likelihood(sim_points[k])),\
								log(pdf_as_pion->get_single_likelihood(sim_points[k]))-log(pdf_as_kaon->get_single_likelihood(sim_points[k])));

						}
					}
				}

				timing_clock = clock() - timing_clock;

				simulation_time->Fill(((float)timing_clock)/(CLOCKS_PER_SEC)*1000);


				// 				printf("\nPID[n]=%i loglikehood difference(pion-kaon)=%f\n",PID[n],ll_diff);
				if(abs(PID[n])==8||PID[n]==9){
					ll_diff_pion->Fill(ll_diff);

					if (ll_diff > 0)
					{
						zero_pion_id->Fill(1);
					}
					else
					{
						zero_pion_id->Fill(0);
					}

					phot_found_pion->Fill(sim_points.size());
				}
				else{

					ll_diff_kaon->Fill(ll_diff);
					if (ll_diff < 0)
					{
						zero_kaon_id->Fill(1);
					}
					else
					{
						zero_kaon_id->Fill(0);
					}
					phot_found_kaon->Fill(sim_points.size());
				}
			}
		}
		printf("\n");//ensure we're on a new line at the end
		printf("Kaon missID: %12.04f%%\n",100.0*zero_kaon_id->GetBinContent(1)/(zero_kaon_id->GetBinContent(1)+zero_kaon_id->GetBinContent(2)));
		printf("Pion missID: %12.04f%%\n",100.0*zero_pion_id->GetBinContent(1)/(zero_pion_id->GetBinContent(1)+zero_pion_id->GetBinContent(2)));
		//end mc reading script

	}
	else if (slac_run == true)
	{ 
		printf("Testing the SLAC geometry/run\n");

		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);

		std::vector<dirc_point> hit_points_pion;
		std::vector<dirc_point> hit_points_kaon;

		dirc_model->set_liquid_absorbtion(liquid_absorbtion);
		dirc_model->set_liquid_index(liquid_index);
		dirc_model->set_three_seg_mirror(three_seg_mirror);
		dirc_model->set_sidemirror(sm_xr,sm_xl);


		dirc_model->set_pmt_offset(pmt_offset);
		dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
		dirc_model->set_bar_box_angle(bar_box_box_angle);

		//ns
		double pion_time = particle_flight_distance/(pion_beta*.3);
		double kaon_time = particle_flight_distance/(kaon_beta*.3);

		dirc_model->sim_reg_n_photons(\
				hit_points_pion,\
				n_phi_phots/10,\
				n_z_phots,\
				pion_angle,\
				1,\
				particle_x,\
				particle_y,\
				pion_time,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/pdf_unc_red_fac,\
				pion_beta);

		dirc_model->sim_reg_n_photons(\
				hit_points_kaon,\
				n_phi_phots/10,\
				n_z_phots,\
				pion_angle,\
				1,\
				particle_x,\
				particle_y,\
				kaon_time,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/pdf_unc_red_fac,\
				kaon_beta);


		DircSpreadGaussian* pdf_pion = new DircSpreadGaussian(\
				sfunc_sig,\
				hit_points_pion,\
				s_func_x,\
				s_func_y,\
				s_func_t);
		DircSpreadGaussian* pdf_kaon = new DircSpreadGaussian(\
				sfunc_sig,\
				hit_points_kaon,\
				s_func_x,\
				s_func_y,\
				s_func_t);

		for (int i = 0; i < num_runs; i++)
		{
			dirc_model->set_focus_mirror_angle(\
					spread_ang.Gaus(main_mirror_angle,mirror_angle_change_unc),\
					spread_ang.Gaus(0,mirror_angle_change_yunc));
			dirc_model->set_upper_wedge_angle_diff(\
					spread_ang.Gaus(0,wedge_uncertainty),\
					spread_ang.Gaus(0,upper_wedge_yang_spread));
			dirc_model->set_bar_box_angle(spread_ang.Gaus(0,box_rot_unc));

			printf("\r                                                    ");
			printf("\rrunning iter %8d/%d  ",i+1,num_runs);


			fflush(stdout);

			particle_theta = spread_ang.Uniform(0,box_rot_unc);
			particle_phi = spread_ang.Uniform(0,360);
			//Sim new photons
			
			//ns
			double pion_time = particle_flight_distance/(pion_beta*.3);
			double kaon_time = particle_flight_distance/(kaon_beta*.3);

			dirc_model->sim_reg_n_photons(\
					hit_points_pion,\
					n_phi_phots/phi_phots_reduce,\
					n_z_phots,\
					pion_angle,\
					1,\
					particle_x,\
					particle_y,\
					pion_time,\
					particle_theta,\
					particle_phi,\
					0,\
					ckov_unc/pdf_unc_red_fac,\
					pion_beta);

			dirc_model->sim_reg_n_photons(\
					hit_points_kaon,\
					n_phi_phots/phi_phots_reduce,\
					n_z_phots,\
					pion_angle,\
					1,\
					particle_x,\
					particle_y,\
					kaon_time,\
					particle_theta,\
					particle_phi,\
					0,\
					ckov_unc/pdf_unc_red_fac,\
					kaon_beta);

			pdf_pion->set_support(hit_points_pion);
			pdf_kaon->set_support(hit_points_kaon);


			n_sim_phots = spread_ang.Gaus(mean_n_phot,spread_n_phot);

			//assume its a middle bar
			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					pion_angle,\
					1,\
					particle_x,\
					particle_y,\
					pion_time,\
					particle_theta,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					pion_beta);
			digitizer.digitize_points(sim_points);


			llc = pdf_pion->get_log_likelihood(sim_points);
			llf = pdf_kaon->get_log_likelihood(sim_points);


			ll_diff_pion->Fill(1*(llc-llf));

			phot_found_pion->Fill(sim_points.size());

			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					kaon_angle,\
					1,\
					particle_x,\
					particle_y,\
					kaon_time,\
					particle_theta,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					kaon_beta);


			digitizer.digitize_points(sim_points);

			llc = pdf_pion->get_log_likelihood(sim_points);
			llf = pdf_kaon->get_log_likelihood(sim_points);

			ll_diff_kaon->Fill(1*(llc-llf));

			// 		ll_diff_kaon->Fill(sep_pdfs->get_log_likelihood_spread_diff(sim_points));

			phot_found_kaon->Fill(sim_points.size());
		}

		printf("\nRun Completed\n");
		if (coverage_plot == true)
		{
			printf("starting coverage ouput\n");

			for (int i = 0; i < num_cov; i++)
			{
				if ((i+1) % 1000 == 0)
				{
					//printf("\r                                                                                                                           ");
					printf("\rcoverage iter %d/%d",i+1,num_cov);
				}
				fflush(stdout);

				particle_theta = spread_ang.Uniform(0,14);
				particle_phi = spread_ang.Uniform(0,360);
				energy = spread_ang.Uniform(2,4.5);
				pion_beta = dirc_model->get_beta(energy,pimass);
				kaon_beta = dirc_model->get_beta(energy,kmass);
				dirc_model->sim_rand_n_photons(\
						sim_points,\
						n_sim_phots,\
						pion_angle,\
						1,\
						particle_x,\
						particle_y,\
						0,\
						particle_theta,\
						particle_phi,\
						tracking_unc,\
						ckov_unc,\
						pion_beta);

				for (unsigned int j = 0; j < sim_points.size(); j++)
				{
					//printf("%12.04f %12.04f\n",sim_points[j].x,sim_points[j].y);
					pion_coverage_xy->Fill(sim_points[j].x+outcsv_x,sim_points[j].y);
				}

				dirc_model->sim_rand_n_photons(\
						sim_points,\
						n_sim_phots,\
						kaon_angle,\
						1,\
						particle_x,\
						particle_y,\
						0,\
						particle_theta,\
						particle_phi,\
						tracking_unc,\
						ckov_unc,\
						kaon_beta);

				for (unsigned int j = 0; j < sim_points.size(); j++)
				{
					kaon_coverage_xy->Fill(sim_points[j].x+outcsv_x,sim_points[j].y);
				}

			}
			printf("\nfinished output of coverage plots\n");
		}
		else
		{
			pion_coverage_xy->Fill(0.0,0.0);
			kaon_coverage_xy->Fill(0.0,0.0);

		}

		//	printf("Found %d pion points on the target\n", (int) hit_points_pion.size());
		double x,y,t_ns;
		for (unsigned int i = 0; i < hit_points_pion.size(); i++)
		{
			x = hit_points_pion[i].x;
			y = hit_points_pion[i].y;
			t_ns = hit_points_pion[i].t;
			// 		if (hit_points_pion[i].t > 0) continue;
			pion_dist_xy->Fill(x,y);
			pion_dist_xt->Fill(x,t_ns);
			pion_dist_yt->Fill(y,t_ns);
			pion_dist_t->Fill(t_ns);
		}

		//	printf("Found %d kaon points on the target\n", (int) hit_points_kaon.size());
		for (unsigned int i = 0; i < hit_points_kaon.size(); i++)
		{
			x = hit_points_kaon[i].x;
			y = hit_points_kaon[i].y;
			t_ns = hit_points_kaon[i].t;
			// 		if (hit_points_pion[i].t > 0) continue;
			kaon_dist_xy->Fill(x,y);
			kaon_dist_xt->Fill(x,t_ns);
			kaon_dist_yt->Fill(y,t_ns);
			kaon_dist_t->Fill(t_ns);
		}
	}
	else
	{ 
		printf("no input file specified.  Running in loop mode\n");

		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);

		if (monochrome_plot == true)
		{
			pion_angle = rad_to_deg*acos(1/(refrac_index*pion_beta));
			kaon_angle = rad_to_deg*acos(1/(refrac_index*kaon_beta));
			pion_beta = -1;
			kaon_beta = -1;
			ckov_unc = 0;
		}

		std::vector<dirc_point> hit_points_pion;
		std::vector<dirc_point> hit_points_kaon;



		// 	dirc_model->set_upper_wedge_angle_diff(0);

		dirc_model->set_liquid_absorbtion(liquid_absorbtion);
		dirc_model->set_liquid_index(liquid_index);
		dirc_model->set_three_seg_mirror(three_seg_mirror);
		dirc_model->set_sidemirror(sm_xr,sm_xl);


		dirc_model->set_pmt_offset(pmt_offset);
		dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
		dirc_model->set_bar_box_angle(bar_box_box_angle);

		//ns
		double pion_time = particle_flight_distance/(pion_beta*.3);
		double kaon_time = particle_flight_distance/(kaon_beta*.3);

		dirc_model->sim_reg_n_photons(\
				hit_points_pion,\
				n_phi_phots,\
				n_z_phots,\
				pion_angle,\
				1,\
				particle_x,\
				particle_y,\
				pion_time,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/pdf_unc_red_fac,\
				pion_beta);

		dirc_model->sim_reg_n_photons(\
				hit_points_kaon,\
				n_phi_phots,\
				n_z_phots,\
				kaon_angle,\
				1,\
				particle_x,\
				particle_y,\
				kaon_time,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/pdf_unc_red_fac,\
				kaon_beta);

		if (flatten_time == true)
		{
			for (unsigned int i = 0; i < hit_points_pion.size(); i++)
			{
				hit_points_pion[i].t = 0;
			}
			for (unsigned int i = 0; i < hit_points_kaon.size(); i++)
			{
				hit_points_kaon[i].t = 0;
			}
		}


		DircSpreadGaussian* pdf_pion = new DircSpreadGaussian(\
				sfunc_sig,\
				hit_points_pion,\
				s_func_x,\
				s_func_y,\
				s_func_t);
		DircSpreadGaussian* pdf_kaon = new DircSpreadGaussian(\
				sfunc_sig,\
				hit_points_kaon,\
				s_func_x,\
				s_func_y,\
				s_func_t);

		for (int i = 0; i < num_runs; i++)
		{
			dirc_model->set_focus_mirror_angle(\
					spread_ang.Gaus(main_mirror_angle,mirror_angle_change_unc),\
					spread_ang.Gaus(0,mirror_angle_change_yunc));
			dirc_model->set_upper_wedge_angle_diff(\
					spread_ang.Gaus(0,wedge_uncertainty),\
					spread_ang.Gaus(0,upper_wedge_yang_spread));
			dirc_model->set_bar_box_angle(spread_ang.Gaus(0,box_rot_unc));

			if (particle_theta_mean < .01)
			{
				particle_phi = spread_ang.Uniform(0,360);
			}

			particle_theta = spread_ang.Gaus(particle_theta_mean, particle_theta_spread);
			particle_x = spread_ang.Gaus(particle_x_mean, particle_x_spread);
			particle_y = spread_ang.Gaus(particle_y_mean, particle_y_spread);
			energy = spread_ang.Gaus(energy_mean,energy_spread);
			pion_beta = dirc_model->get_beta(energy,pimass);
			kaon_beta = dirc_model->get_beta(energy,kmass);

			if (monochrome_plot == true)
			{
				pion_angle = rad_to_deg*acos(1/(refrac_index*pion_beta));
				kaon_angle = rad_to_deg*acos(1/(refrac_index*kaon_beta));
				pion_beta = -1;
				kaon_beta = -1;
				ckov_unc = 0;
			}
			//			printf("%12.04f %12.04f %12.04f %12.04f\n",particle_x, particle_y, particle_theta, energy);

			printf("\r                                                    ");
			printf("\rrunning iter %8d/%d  ",i+1,num_runs);


			fflush(stdout);

			n_sim_phots = spread_ang.Gaus(mean_n_phot,spread_n_phot);

			//assume its a middle bar
			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					pion_angle,\
					1,\
					particle_x,\
					particle_y,\
					pion_time,\
					particle_theta+const_track_off,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					pion_beta);
			digitizer.digitize_points(sim_points);

			if (flatten_time == true)
			{
				for (unsigned int i = 0; i < sim_points.size(); i++)
				{
					sim_points[i].t = 0;
				}
			}
			llc = pdf_pion->get_log_likelihood(sim_points);
			llf = pdf_kaon->get_log_likelihood(sim_points);

			ll_diff_pion->Fill(1*(llc-llf));
			phot_found_pion->Fill(sim_points.size());

			if (sep_updown == true)
			{
				std::vector<dirc_point> sim_points_up;
				std::vector<dirc_point> sim_points_down;
				for (unsigned int j = 0; j < sim_points.size(); j++)
				{
					if (sim_points[j].updown == 1)
					{
						sim_points_up.push_back(sim_points[j]);
					}
					else if (sim_points[j].updown == -1)
					{
						sim_points_down.push_back(sim_points[j]);
					}
					else
					{
						printf("Found a sim point without an updown flag.  This shouldn't happen.\n");
					}
				}
				//just with the up
				llc = pdf_pion->get_log_likelihood(sim_points_up);
				llf = pdf_kaon->get_log_likelihood(sim_points_up);

				ll_diff_pion_up->Fill(1*(llc-llf));

				//just with the down
				llc = pdf_pion->get_log_likelihood(sim_points_down);
				llf = pdf_kaon->get_log_likelihood(sim_points_down);

				ll_diff_pion_down->Fill(1*(llc-llf));

				phot_found_pion_up->Fill(sim_points_up.size());
				phot_found_pion_down->Fill(sim_points_down.size());
			}




			dirc_model->sim_rand_n_photons(\
					sim_points,\
					n_sim_phots,\
					kaon_angle,\
					1,\
					particle_x,\
					particle_y,\
					kaon_time,\
					particle_theta+const_track_off,\
					particle_phi,\
					tracking_unc,\
					ckov_unc,\
					kaon_beta);

			if (flatten_time == true)
			{
				for (unsigned int i = 0; i < sim_points.size(); i++)
				{
					sim_points[i].t = 0;
				}
			}

			digitizer.digitize_points(sim_points);

			llc = pdf_pion->get_log_likelihood(sim_points);
			llf = pdf_kaon->get_log_likelihood(sim_points);

			ll_diff_kaon->Fill(1*(llc-llf));
			phot_found_kaon->Fill(sim_points.size());

			if (sep_updown == true)
			{
				std::vector<dirc_point> sim_points_up;
				std::vector<dirc_point> sim_points_down;
				for (unsigned int j = 0; j < sim_points.size(); j++)
				{
					if (sim_points[j].updown == 1)
					{
						sim_points_up.push_back(sim_points[j]);
					}
					else if (sim_points[j].updown == -1)
					{
						sim_points_down.push_back(sim_points[j]);
					}
					else
					{
						printf("Found a sim point without an updown flag.  This shouldn't happen.\n");
					}
				}
				//just with the up
				llc = pdf_pion->get_log_likelihood(sim_points_up);
				llf = pdf_kaon->get_log_likelihood(sim_points_up);

				ll_diff_kaon_up->Fill(1*(llc-llf));

				//just with the down
				llc = pdf_pion->get_log_likelihood(sim_points_down);
				llf = pdf_kaon->get_log_likelihood(sim_points_down);

				ll_diff_kaon_down->Fill(1*(llc-llf));
				phot_found_kaon_up->Fill(sim_points_up.size());
				phot_found_kaon_down->Fill(sim_points_down.size());
			}
			// 		ll_diff_kaon->Fill(sep_pdfs->get_log_likelihood_spread_diff(sim_points));

		}

		printf("\nRun Completed\n");
		if (coverage_plot == true)
		{
			printf("starting coverage ouput\n");

			for (int i = 0; i < num_cov; i++)
			{
				if ((i+1) % 1000 == 0)
				{
					//printf("\r                                                                                                                           ");
					printf("\rcoverage iter %d/%d",i+1,num_cov);
				}
				fflush(stdout);

				particle_theta = spread_ang.Uniform(0,14);
				particle_phi = spread_ang.Uniform(0,360);
				energy = spread_ang.Uniform(2,4.5);
				pion_beta = dirc_model->get_beta(energy,pimass);
				kaon_beta = dirc_model->get_beta(energy,kmass);
				dirc_model->sim_rand_n_photons(\
						sim_points,\
						n_sim_phots,\
						pion_angle,\
						1,\
						particle_x,\
						particle_y,\
						0,\
						particle_theta,\
						particle_phi,\
						tracking_unc,\
						ckov_unc,\
						pion_beta);

				for (unsigned int j = 0; j < sim_points.size(); j++)
				{
					//printf("%12.04f %12.04f\n",sim_points[j].x,sim_points[j].y);
					pion_coverage_xy->Fill(sim_points[j].x+outcsv_x,sim_points[j].y);
				}

				dirc_model->sim_rand_n_photons(\
						sim_points,\
						n_sim_phots,\
						kaon_angle,\
						1,\
						particle_x,\
						particle_y,\
						0,\
						particle_theta,\
						particle_phi,\
						tracking_unc,\
						ckov_unc,\
						kaon_beta);

				for (unsigned int j = 0; j < sim_points.size(); j++)
				{
					kaon_coverage_xy->Fill(sim_points[j].x+outcsv_x,sim_points[j].y);
				}

			}
			printf("\nfinished output of coverage plots\n");
		}
		else
		{
			pion_coverage_xy->Fill(0.0,0.0);
			kaon_coverage_xy->Fill(0.0,0.0);

		}
		if (fill_kinematics_yields == true)
		{
			int itheta_bin = 0;
			int iphi_bin = 0;
			printf("starting kinematics yields ouput\n");
			for (double itheta = kinematics_yields_min_theta + kinematics_yields_step_theta/2; itheta < kinematics_yields_max_theta; itheta += kinematics_yields_step_theta)
			{
				itheta_bin++;//I know there's a better way to do this, but I'm lazy
				iphi_bin = 0;
				for (double iphi = kinematics_yields_min_phi + kinematics_yields_step_phi/2; iphi < kinematics_yields_max_phi; iphi += kinematics_yields_step_phi)
				{
					iphi_bin++;
					printf("\r                                                                             ");
					printf("\rPhi: %12.04f    Theta: %12.04f",iphi,itheta);
					fflush(stdout);	
					pion_beta = dirc_model->get_beta(5,pimass); //High energy
					dirc_model->sim_rand_n_photons(
							sim_points,
							kinematics_n_phots,
							pion_angle,
							1,
							0, //particle_x
							0, //particle_y
							0,\
							itheta,
							iphi,
							0, //ignore tracking
							ckov_unc,
							pion_beta);

					kinematics_yields->SetBinContent(itheta_bin,iphi_bin,(1.0*sim_points.size())/kinematics_n_phots);
				}
			}
			printf("\n");
		}


		//Make the distribution - use random	
		pion_beta = dirc_model->get_beta(energy,pimass);
		kaon_beta = dirc_model->get_beta(energy,kmass);
		if (monochrome_plot == true)
		{
			pion_angle = rad_to_deg*acos(1/(refrac_index*pion_beta));
			kaon_angle = rad_to_deg*acos(1/(refrac_index*kaon_beta));
			pion_beta = -1;
			kaon_beta = -1;
			ckov_unc = 0;
		}
		/*
		   dirc_model->sim_rand_n_photons(\
		   hit_points_pion,\
		   n_phi_phots*n_z_phots,\
		   pion_angle,\
		   1,\
		   particle_x,\
		   particle_y,\
		   particle_theta,\
		   particle_phi,\
		   0,\
		   ckov_unc/pdf_unc_red_fac,\
		   pion_beta);

		   dirc_model->sim_reg_n_photons(\
		   hit_points_kaon,\
		   n_phi_phots*n_z_phots,\
		   kaon_angle,\
		   1,\
		   particle_x,\
		   particle_y,\
		   particle_theta,\
		   particle_phi,\
		   0,\
		   ckov_unc/pdf_unc_red_fac,\
		   kaon_beta);
		 */
		//	printf("Found %d pion points on the target\n", (int) hit_points_pion.size());
		double x,y,t_ns;
		for (unsigned int i = 0; i < hit_points_pion.size(); i++)
		{
			x = hit_points_pion[i].x;
			y = hit_points_pion[i].y;
			t_ns = hit_points_pion[i].t;
			// 		if (hit_points_pion[i].t > 0) continue;
			pion_dist_x->Fill(x);
			pion_dist_y->Fill(y);
			pion_dist_t->Fill(t_ns);
			pion_dist_xy->Fill(x,y);
			pion_dist_xt->Fill(x,t_ns);
			pion_dist_yt->Fill(y,t_ns);
			pion_dist_t->Fill(t_ns);
		}

		//	printf("Found %d kaon points on the target\n", (int) hit_points_kaon.size());
		for (unsigned int i = 0; i < hit_points_kaon.size(); i++)
		{
			x = hit_points_kaon[i].x;
			y = hit_points_kaon[i].y;
			t_ns = hit_points_kaon[i].t;
			// 		if (hit_points_pion[i].t > 0) continue;
			kaon_dist_x->Fill(x);
			kaon_dist_y->Fill(y);
			kaon_dist_t->Fill(t_ns);
			kaon_dist_xy->Fill(x,y);
			kaon_dist_xt->Fill(x,t_ns);
			kaon_dist_yt->Fill(y,t_ns);
			kaon_dist_t->Fill(t_ns);
		}
	}
	/*
	   if (out_csv == true)
	   {
	   ofstream pion_csv;
	   ofstream kaon_csv;

	   pion_csv.open("pion_dist.csv");
	   kaon_csv.open("kaon_dist.csv");
	   char format_buffer[256];

	   for (unsigned int i = 0; i < hit_points_pion.size(); i++)
	   {
	   x = hit_points_pion[i].x;
	   y = hit_points_pion[i].y;
	   t_ns = hit_points_pion[i].t;
	   sprintf(format_buffer,"%12.06e %12.06e %12.06e\n",x+outcsv_x,y,t_ns+outcsv_t);

	   pion_csv << format_buffer;
	   }
	   for (unsigned int i = 0; i < hit_points_kaon.size(); i++)
	   {
	   x = hit_points_kaon[i].x;
	   y = hit_points_kaon[i].y;
	   t_ns = hit_points_kaon[i].t;
	   sprintf(format_buffer,"%12.06e %12.06e %12.06e\n",x+outcsv_x,y,t_ns+outcsv_t);

	   kaon_csv << format_buffer;
	   }	
	   pion_csv.close();
	   kaon_csv.close();
	   }//
	   std::vector<double> dist_traveled = dirc_model->get_dist_traveled();
	   for (unsigned int i = 0; i < dist_traveled.size(); i++)
	   {
	   liquid_dist->Fill(dist_traveled[i]);
	   }
	 */
	double dist_mean = 0;
	if (box_check_n > 0)
	{	
		printf("Box Checking Histograms...\n");
		std::vector<dirc_point> box_check_points;

		dirc_model->test_from_wedge_top(\
				box_check_points,\
				box_check_n,\
				1 /*Bar*/,\
				particle_x,\
				box_check_theta,\
				box_check_phi,\
				box_check_theta_unc,\
				box_check_phi_unc,\
				box_check_overall_theta);

		dist_mean = 0;

		for (unsigned int i = 0; i < box_check_points.size(); i++)
		{
			box_check_x->Fill(box_check_points[i].x);	
			box_check_y->Fill(box_check_points[i].y);	
			box_check_d->Fill(box_check_points[i].t);
			dist_mean += box_check_points[i].t;
			box_check_xy->Fill(box_check_points[i].x,box_check_points[i].y);
		}
		dist_mean /= box_check_points.size();

		printf("Box Check Results\n");
		printf("n-phots XSpread YSpread Dist X/Dist Y/Dist\n");
		printf("%07d %12.04f %12.04f %12.04f %12.04e %12.04e\n",(int) box_check_points.size(),box_check_x->GetRMS(),box_check_y->GetRMS(),dist_mean,box_check_x->GetRMS()/dist_mean,box_check_y->GetRMS()/dist_mean);


	}
	if (phot_check_n > 0)
	{	
		printf("Checking Number of Photons Found\n");
		std::vector<dirc_point> phot_check_points;

		int num_phot_points = 100;
		double phot_check_theta = 0;
		double phot_check_phi = 0;
		dist_mean = 0;	

		for (int i = 0; i < phot_check_n; i++)
		{
			pion_beta = dirc_model->get_beta(energy,pimass);
			kaon_beta = dirc_model->get_beta(energy,kmass);
			if (monochrome_plot == true)
			{
				pion_angle = rad_to_deg*acos(1/(refrac_index*pion_beta));
				kaon_angle = rad_to_deg*acos(1/(refrac_index*kaon_beta));
				pion_beta = -1;
				kaon_beta = -1;
				ckov_unc = 0;
			}
			phot_check_theta = spread_ang.Uniform(0,phot_check_max_theta);
			phot_check_phi = spread_ang.Uniform(0,360);
			dirc_model->sim_rand_n_photons(\
					phot_check_points,\
					num_phot_points,\
					pion_angle,\
					1,\
					particle_x,\
					particle_y,\
					0,\
					phot_check_theta,\
					phot_check_phi,\
					tracking_unc,\
					ckov_unc,\
					pion_beta);

			dist_mean += phot_check_points.size();
			//printf("%d\n",phot_check_points.size());
			phot_check_points.clear();
		}

		printf("Relative photons found:\n");
		printf("%12.04f\n",dist_mean/(num_phot_points*phot_check_n));


	}
	if (refraction_sim_n > 0)
	{
		std::vector<std::pair<double, double> > pion_theta_cphi;
		std::vector<std::pair<double, double> > kaon_theta_cphi;

		particle_theta=particle_theta_mean;	
		particle_phi=particle_phi_mean;	

		printf("Creating Refraction histograms...\n");
		std::vector<double> before_fill;
		std::vector<double> after_fill;
		std::vector<double> pmt_incidence;

		pion_theta_cphi = dirc_model->get_refraction_rand_phi(\
				before_fill,\
				after_fill,\
				pmt_incidence,\
				refraction_sim_n,\
				pion_angle ,
				particle_x,\
				particle_y,\
				particle_theta,\
				particle_phi,\
				0,\
				0,\
				pion_beta);
		//Since we're running this all at once, ignore tracking and other uncertainties

		for (unsigned int i = 0; i < before_fill.size(); i++)
		{
			ref_pion_before->Fill(before_fill[i]*rad_to_deg);
		}
		for (unsigned int i = 0; i < after_fill.size(); i++)
		{
			ref_pion_after->Fill(after_fill[i]*rad_to_deg);
		}
		for (unsigned int i = 0; i < pmt_incidence.size(); i++)
		{
			ref_pion_sens_plane->Fill(pmt_incidence[i]*rad_to_deg);
		}
		for (unsigned int i = 0; i < pion_theta_cphi.size(); i++)
		{
			ref_theta_cphi_pion->Fill(pion_theta_cphi[i].first*rad_to_deg,pion_theta_cphi[i].second*rad_to_deg);
		}

		kaon_theta_cphi = dirc_model->get_refraction_rand_phi(\
				before_fill,\
				after_fill,\
				pmt_incidence,\
				refraction_sim_n,\
				kaon_angle ,
				particle_x,\
				particle_y,\
				particle_theta,\
				particle_phi,\
				0,\
				0,\
				kaon_beta);


		for (unsigned int i = 0; i < before_fill.size(); i++)
		{
			ref_kaon_before->Fill(before_fill[i]*rad_to_deg);
		}
		for (unsigned int i = 0; i < after_fill.size(); i++)
		{
			ref_kaon_after->Fill(after_fill[i]*rad_to_deg);
		}
		for (unsigned int i = 0; i < pmt_incidence.size(); i++)
		{
			ref_kaon_sens_plane->Fill(pmt_incidence[i]*rad_to_deg);
		}
		for (unsigned int i = 0; i < kaon_theta_cphi.size(); i++)
		{
			ref_theta_cphi_kaon->Fill(kaon_theta_cphi[i].first*rad_to_deg,kaon_theta_cphi[i].second*rad_to_deg);
		}
		printf("Refraction histograms created\n");
	}


	double t_lambda = 0;
	for (int i = 0; i < 1000000; i++)
	{
		//4GeV pion
		pion_cerenkov->Fill(dirc_model->get_cerenkov_angle_rand(.9994,0,t_lambda));
		pion_lambda->Fill(t_lambda);
		//4GeV kaon
		kaon_cerenkov->Fill(dirc_model->get_cerenkov_angle_rand(.99235,0,t_lambda));
		kaon_lambda->Fill(t_lambda);
	}


	tfile->cd();

	pion_dist_geant_xy->Write();
	pion_dist_geant_xt->Write();
	pion_dist_geant_yt->Write();
	pion_dist_geant_x->Write();
	pion_dist_geant_y->Write();
	pion_dist_geant_t->Write();

	pion_dist_x->Write();
	pion_dist_y->Write();
	kaon_dist_x->Write();
	kaon_dist_y->Write();

	box_check_x->Write();
	box_check_y->Write();
	box_check_d->Write();
	box_check_xy->Write();

	pion_dist_xy->Write();
	kaon_dist_xy->Write();
	pion_dist_xt->Write();
	kaon_dist_xt->Write();
	pion_dist_yt->Write();
	kaon_dist_yt->Write();
	pion_dist_t->Write();
	kaon_dist_t->Write();
	ll_diff_pion->Write();
	ll_diff_kaon->Write();
	ll_diff_pion_up->Write();
	ll_diff_kaon_up->Write();
	ll_diff_pion_down->Write();
	ll_diff_kaon_down->Write();
	phot_found_pion->Write();
	phot_found_kaon->Write();
	phot_found_pion_up->Write();
	phot_found_kaon_up->Write();
	phot_found_pion_down->Write();
	phot_found_kaon_down->Write();
	pion_cerenkov->Write();
	kaon_cerenkov->Write();
	pion_coverage_xy->Write();
	kaon_coverage_xy->Write();
	kinematics_yields->Write();
	liquid_dist->Write();
	pion_lambda->Write();
	kaon_lambda->Write();


	ref_pion_before->Write();
	ref_pion_after->Write();
	ref_kaon_before->Write();
	ref_kaon_after->Write();
	ref_pion_sens_plane->Write();
	ref_kaon_sens_plane->Write();


	ref_theta_cphi_pion->Write();
	ref_theta_cphi_kaon->Write();

	simulation_time->Write();
	tfile->Close();

	int status = 0;
	return status;

}
