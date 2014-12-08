#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

#include <stdlib.h>
#include <math.h>

#include "include/dirc_optical_sim.h"
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
//
//
//
//
//
//cdd work on faster acos calculation, that is getting angle from beta (most of it comes out at 45 deg.....)

//
//
//
//
//cdd move this fold_x(inpoints) to header?
std::vector<dirc_point> fold_x(std::vector<dirc_point> inpoints) {
    std::vector<dirc_point> outvec;
    for (unsigned int i = 0; i < inpoints.size(); i++) {
        dirc_point tpoint = inpoints[i];
        tpoint.x = fabs(tpoint.x);
        outvec.push_back(tpoint);
    }
    return outvec;
}

int main(int nargs, char* argv[]) {
    double in_num = 0;
    const char* in_str;
    bool inputfile =true;
    bool out_csv=false;
    double time_window=10;//time window for compounded pmt hits, in ns
    printf("Arguments Passed=%d\n",nargs);

    if(nargs==1) {
        inputfile=false;
        out_csv=true;
    };
    if(nargs==2) {
        in_str = argv[1];
        printf("%s\n",in_str);
        printf("nargs=1");
    };
    if(nargs==3) {
        in_str = argv[1];
        time_window= atof(argv[2]);
        printf("Opening %s with time window of %fns\n",in_str,time_window);
    };


    double resx = 6;
    double resy = 6;
    double rest = 1;
    double minx = -8000;
    double maxx = -minx;
    double miny = -800;
    double maxy = -miny;
    double mint = 0;
    double maxt = 1000;
    double t_unc = 1;
// 	double t_unc = .05;


    double rad_to_deg = 57.2958;

    double res_enhance = 1;

    double tracking_unc = .0000*57.3; //mrad
// 	double ckov_unc = .0077*57.3; //chromatic + optical aberation = 7.7mrad
    double ckov_unc = .003*57.3; //transport = 3mrad


    double energy = 5.0;
    double kmass = .4937;
    double pimass = .1396;
    double mumass = .1057;

    double particle_x = 0;
    double particle_y = 0;
    double particle_theta = 4;
    double particle_phi = 40;

    int num_runs = 1000;

    int n_sim_phots = 40;

    int refraction_sim_n = 0;

    int n_phi_phots = 100000;
// 	int n_phi_phots =2;
    int n_z_phots = 4;
    double s_func_x = 6;
    double s_func_y = s_func_x;
    double s_func_t = 2;
// 	double s_func_t = .3;
    double sfunc_sig = 1;

    double outcsv_x,outcsv_y,outcsv_t;
    outcsv_x = 0*35;//bars are 35mm wide
    outcsv_y = 0;//mm
    outcsv_t = 0;//ns


    if(out_csv) {
        n_phi_phots = 4000;
        n_z_phots = 8;
        particle_y += outcsv_y;
    }
    double pdf_unc_red_fac = 1;
    double wedge_uncertainty = 0/57.3;
    double mirror_angle_change = 0;
    double mirror_angle_change_unc = in_num;
    double mirror_angle_change_yunc = 0;
    double box_rot = 0;
    double box_rot_unc = 0;
    double bar_box_box_angle = 0/57.3;
    double mirror_r_difference = -400;
    double wedge_non_uniformity = 0;
    double pmt_offset = 0;
    double main_mirror_nonuniformity = 0;

    double upper_wedge_yang_spread = 0;
    int rseed = 0;

    double sm_xl = -10000000;
    double sm_xr = -sm_xl;

    double overlap_x = -1;

    bool three_seg_mirror = true;


    double liquid_absorbtion = 0*-log(.7)/1000;
    double liquid_index = 1.33;

    bool coverage_plot = false;
    int num_cov = 100000;

    bool up_down_sep = false;
    TRandom3 spread_ang(rseed+3);


    DircOpticalSim *dirc_model = new DircOpticalSim(\
            rseed,\
            -1200 + mirror_r_difference,\
            300.38,\
            74.11 + mirror_angle_change,\
            600,\
            47.87 + box_rot);

    double muon_beta, pion_beta, kaon_beta/*, electron_beta:=1*/;
    double muon_angle, pion_angle, kaon_angle;

    char* rootfilename = new char[256];
    sprintf(rootfilename,"fitdirc.root");

    TFile* tfile = new TFile(rootfilename,"RECREATE");

    TH1F *ll_diff_pion = new TH1F("ll_diff_pion","Difference of log likelihood real = pion",12000,-1500,1500);
    TH1F *ll_diff_kaon = new TH1F("ll_diff_kaon","Difference of log likelihood real = kaon",12000,-1500,1500);
    TH1F *phot_found_pion = new TH1F("phot_found_pion","number of photons found on pion angle", 1001,-.5,1000.5);
    TH1F *phot_found_kaon = new TH1F("phot_found_kaon","number of photons found on kaon angle", 1001,-.5,1000.5);

    TH1F *ref_pion_before = new TH1F("ref_pion_before","Angle of Pion Photons going into interface", 9000,0,90);
    TH1F *ref_kaon_before = new TH1F("ref_kaon_before","Angle of Kaon Photons going into interface", 9000,0,90);
    TH1F *ref_pion_after = new TH1F("ref_pion_after","Angle of Pion Photons going out of interface", 9000,0,90);
    TH1F *ref_kaon_after = new TH1F("ref_kaon_after","Angle of Kaon Photons going out of interface", 9000,0,90);

    TH1F *pion_dist_t = new TH1F("pion_dist_t","t val of intercepted points - pion",(maxt-mint)/(res_enhance*rest),mint,maxt);
    TH1F *kaon_dist_t = new TH1F("kaon_dist_t","t val of intercepted points - kaon",(maxt-mint)/(res_enhance*rest),mint,maxt);


    TH2F *ref_theta_cphi_pion = new TH2F("ref_theta_cphi_pion","Emitted angle of Pion Photons versus angle into interface", 7200,0,360,1800,0,90);
    TH2F *ref_theta_cphi_kaon = new TH2F("ref_theta_cphi_kaon","Emitted angle of Kaon Photons versus angle into interface", 7200,0,360,1800,0,90);

    TH2F *pion_dist_xy = new TH2F("pion_dist_xy","xy val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
    TH2F *kaon_dist_xy = new TH2F("kaon_dist_xy","xy val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);

    TH2F *pion_dist_xt = new TH2F("pion_dist_xt","xt val of intercepted points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);
    TH2F *kaon_dist_xt = new TH2F("kaon_dist_xt","xt val of intercepted points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxt-mint)/(res_enhance*rest),mint,maxt);

    TH2F *pion_dist_yt = new TH2F("pion_dist_yt","yt val of intercepted points - pion",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);
    TH2F *kaon_dist_yt = new TH2F("kaon_dist_yt","yt val of intercepted points - kaon",(maxy-miny)/(res_enhance*resy),miny,maxy,(maxt-mint)/(res_enhance*rest),mint,maxt);


    maxy *= 5;

    TH2F *pion_coverage_xy = new TH2F("pion_coverage_xy","xy val of generated points - pion",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);
    TH2F *kaon_coverage_xy = new TH2F("kaon_coverage_xy","xy val of generated points - kaon",(maxx-minx)/(res_enhance*resx),minx,maxx,(maxy-miny)/(res_enhance*resy),miny,maxy);


    TH1F *pion_cerenkov = new TH1F("pion_cerenkov","pion cerenkov angle distribution",300,45,48);
    TH1F *kaon_cerenkov = new TH1F("kaon_cerenkov","kaon cerenkov angle distribution",300,45,48);
    TH1F *pion_lambda = new TH1F("pion_lambda","pion wavelength distribution",450,250,700);
    TH1F *kaon_lambda = new TH1F("kaon_lambda","kaon wavelength distribution",450,250,700);

    TH1F *liquid_dist = new TH1F("liquid_dist","distance travelled in liquid (mm)",1500,0,1500);
    dirc_model ->set_store_traveled(true);

    /*	for (int i = 0; i < 1000000; i++)
    	{
    		double tmp;
    		pion_cerenkov->Fill(dirc_model->get_cerenkov_angle_rand(pion_beta,ckov_unc, tmp));
    		pion_lambda->Fill(tmp);
    		kaon_cerenkov->Fill(dirc_model->get_cerenkov_angle_rand(kaon_beta,ckov_unc, tmp));
    		kaon_lambda->Fill(tmp);
    	}
    //*/



    DircDigitizer digitizer(\
                            minx,\
                            maxx,\
                            resx,\
                            miny,\
                            maxy,\
                            resy,\
                            t_unc);
    /*
    std::vector<dirc_point> hit_points_pion;
    std::vector<dirc_point> hit_points_kaon;

    // 	dirc_model->set_wedge_mirror_rand(wedge_non_uniformity);

    // 	dirc_model->set_upper_wedge_angle_diff(0);

    dirc_model->set_liquid_absorbtion(liquid_absorbtion);
    dirc_model->set_liquid_index(liquid_index);
    dirc_model->set_three_seg_mirror(three_seg_mirror);
    dirc_model->set_sidemirror(sm_xr,sm_xl);


    dirc_model->set_pmt_offset(pmt_offset);
     		dirc_model->set_upper_wedge_angle_diff(wedge_uncertainty);
    dirc_model->set_bar_box_angle(bar_box_box_angle);


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

    DircProbabilitySeparation* sep_pdfs = new DircProbabilitySeparation(pdf_kaon,pdf_pion);
    //  */
    printf("Beginning Run\n");
    double llc, llf;
    std::vector<dirc_point> sim_points;
    std::vector<dirc_point> confound_points;
    dirc_model->set_focmirror_nonuniformity(main_mirror_nonuniformity);
    if(inputfile) {
        int mc_tally=0;
        int confounded_tally=0;
        num_runs=0;
        coverage_plot=false;
        std::ifstream f(in_str);
        if(!f.is_open()) {
            printf("%s not found!",in_str);
            return -1;
        }

        int line_buf_size = 1000;
        char  s[line_buf_size];
        f.getline(s,line_buf_size); //skips first line
        //declare memory addresses for inputs
        //units: mm, ns, deg, GeV
        int iPID, iBAR;
        double ix,iy,it,itheta,iphi,iE;
        std::vector<int> PID,BAR;
        std::vector<double> x,y,t,theta,phi,E;

        //declare for beta and angle based on PID, and constructors for hitpoints and pdfs
        double pion_mc_beta,kaon_mc_beta,pion_mc_angle,kaon_mc_angle;
        std::vector<dirc_point> hits_trk_is_pion;
        std::vector<dirc_point> hits_trk_is_kaon;
        //DircSpreadGaussian  * pdf_as_pion, * pdf_as_kaon;
        DircProbabilitySeparation * sep_pdfs_mc;

        unsigned int r=0;
        while(f>>iPID>>iBAR>>ix>>iy>>it>>itheta>>iphi>>iE) {
            PID.push_back(iPID);
            BAR.push_back(iBAR);
            x.push_back(ix);
            y.push_back(iy);
            t.push_back(it);
            theta.push_back(itheta);
            phi.push_back(iphi);
            E.push_back(iE);
            //std::cout<<PID[r]<<" "<<BAR[r]<<" "<<x[r]<<" "<<y[r]<<" "<<t[r]<<" "<<theta[r]<<" "<<phi[r]<<" "<<E[r]<<std::endl;
            r++;
        }
        f.close();
        for(unsigned int n=0; n<r; n++) {
            if(fabs(PID[n])>7 && fabs(PID[n])<13) { //8 and 9 are pions, 11 and 12 are kaons

                mc_tally++;
                //
                //
                //
                dirc_model->set_focus_mirror_angle(\
                                                   spread_ang.Gaus(74.11,mirror_angle_change_unc),\
                                                   spread_ang.Gaus(0,mirror_angle_change_yunc));
                dirc_model->set_upper_wedge_angle_diff(\
                                                       spread_ang.Gaus(0,0),\
                                                       spread_ang.Gaus(0,0));
                //
                //and replace acos......
                //cdd replace 1.47 with double refrac_index=1.47;
                pion_mc_beta = dirc_model->get_beta(E[n],pimass);
                kaon_mc_beta = dirc_model->get_beta(E[n],kmass);
                pion_mc_angle = rad_to_deg*acos(1/(1.47*pion_mc_beta));
                kaon_mc_angle = rad_to_deg*acos(1/(1.47*kaon_mc_beta));

                hits_trk_is_pion = dirc_model->sim_reg_n_photons(\
                                   n_phi_phots,\
                                   n_z_phots,\
                                   pion_mc_angle,\
                                   BAR[n],\
                                   x[n],\
                                   y[n],\
                                   theta[n],\
                                   phi[n],\
                                   0,\
                                   ckov_unc/pdf_unc_red_fac,\
                                   pion_mc_beta);

                hits_trk_is_kaon = dirc_model->sim_reg_n_photons(\
                                   n_phi_phots,\
                                   n_z_phots,\
                                   kaon_mc_angle,
                                   BAR[n],\
                                   x[n],\
                                   y[n],\
                                   theta[n],\
                                   phi[n],\
                                   0,\
                                   ckov_unc/pdf_unc_red_fac,\
                                   kaon_mc_beta);

                DircSpreadGaussian pdf_as_pion(\
                                               sfunc_sig,\
                                               hits_trk_is_pion,\
                                               s_func_x,\
                                               s_func_y,\
                                               s_func_t);

                DircSpreadGaussian pdf_as_kaon(\
                                               sfunc_sig,\
                                               hits_trk_is_kaon,\
                                               s_func_x,\
                                               s_func_y,\
                                               s_func_t);

                //sep_pdfs_mc = new DircProbabilitySeparation(pdf_as_kaon,pdf_as_pion);

                dirc_model->set_focus_mirror_angle(\
                                                   spread_ang.Gaus(74.11,mirror_angle_change_unc),\
                                                   spread_ang.Gaus(0,mirror_angle_change_yunc));
                dirc_model->set_upper_wedge_angle_diff(\
                                                       spread_ang.Gaus(0,wedge_uncertainty),\
                                                       spread_ang.Gaus(0,upper_wedge_yang_spread));
                printf("\r  ");
                printf("\rRunning monte carlo event  %i/%i. Found %i pions/kaons", n,r-1,mc_tally);

                fflush(stdout);

                //simulate pmt hits from primary mc event to compare with preconstructed pdfs
                //
                if(abs(PID[n])==8 || PID[n]==9 ) {
                    sim_points = dirc_model->sim_rand_n_photons(\
                                 n_sim_phots,\
                                 pion_mc_angle,\
                                 BAR[n],\
                                 x[n],\
                                 y[n],\
                                 theta[n],\
                                 phi[n],\
                                 tracking_unc,\
                                 ckov_unc,\
                                 pion_mc_beta);
                } else {
                    sim_points = dirc_model->sim_rand_n_photons(\
                                 n_sim_phots,\
                                 kaon_mc_angle,\
                                 BAR[n],\
                                 x[n],\
                                 y[n],\
                                 theta[n],\
                                 phi[n],\
                                 tracking_unc,\
                                 ckov_unc,\
                                 kaon_mc_beta);
                }
                //simulate pmt hits from background events that occur within t ns of the primary event
                confounded_tally=0;
                for (unsigned int i =0; i< r; i++) {

                    if(fabs(t[n]-t[i])<time_window && i!=n) {
                        confounded_tally++;
                        if(abs(PID[i])==2 ||PID[i]==3) {

                            confound_points=dirc_model->sim_rand_n_photons(\
                                            n_sim_phots,\
                                            47.135,\
                                            BAR[i],\
                                            x[i],\
                                            y[i],\
                                            theta[i],\
                                            phi[i],\
                                            tracking_unc,\
                                            ckov_unc,\
                                            1);
                        } else if(abs(PID[i])==8||PID[i]==9) {
                            pion_beta = dirc_model->get_beta(E[i],pimass);
                            pion_angle = rad_to_deg*acos(1/(1.47*pion_beta));

                            confound_points=dirc_model->sim_rand_n_photons(\
                                            n_sim_phots,\
                                            pion_angle,\
                                            BAR[i],\
                                            x[i],\
                                            y[i],\
                                            theta[i],\
                                            phi[i],\
                                            tracking_unc,\
                                            ckov_unc,\
                                            pion_beta);

                        } else if(abs(PID[i])==11 || PID[i]==12) {
                            kaon_beta = dirc_model->get_beta(E[i],kmass);
                            kaon_angle = rad_to_deg*acos(1/(1.47*kaon_beta));

                            confound_points=dirc_model->sim_rand_n_photons(\
                                            n_sim_phots,\
                                            kaon_angle,\
                                            BAR[i],\
                                            x[i],\
                                            y[i],\
                                            theta[i],\
                                            phi[i],\
                                            tracking_unc,\
                                            ckov_unc,\
                                            kaon_beta);

                        } else if(abs(PID[i])==5 || PID[i]==6) {
                            muon_beta = dirc_model->get_beta(E[i],mumass);
                            muon_angle = rad_to_deg*acos(1/(1.47*muon_beta));

                            confound_points=dirc_model->sim_rand_n_photons(\
                                            n_sim_phots,\
                                            muon_angle,\
                                            BAR[i],\
                                            x[i],\
                                            y[i],\
                                            theta[i],\
                                            phi[i],\
                                            tracking_unc,\
                                            ckov_unc,\
                                            muon_beta);
                        }
                    }
                    for (unsigned int i = 0; i < confound_points.size(); i++) {
                        sim_points.push_back(confound_points[i]);
                    }
                    confound_points.clear();
                }//end confounded point generation

                digitizer.digitize_points(sim_points);

                llc = pdf_as_pion.get_log_likelihood(sim_points);
                llf = pdf_as_kaon.get_log_likelihood(sim_points);
                //printf("Found %i confounding events for track %i. sim_points.size()=%lu\n",confounded_tally,n,sim_points.size());

            } else {
                continue;
            }


            //printf("\nPID[n]=%i loglikehood difference(pion-kaon)=%f\n",PID[n],llc-llf);
            if(abs(PID[n])==8||PID[n]==9) {
                ll_diff_pion->Fill(llc-llf);

                phot_found_pion->Fill(sim_points.size()/(confounded_tally+1));//divide by tally+1 for an average number of photons per event? or do we want to know the total number of photons in the timewindow?
            } else {

                ll_diff_kaon->Fill(llc-llf);

                phot_found_kaon->Fill(sim_points.size()/(confounded_tally+1));
            }
            //        delete[] pdf_as_kaon;delete[] pdf_as_pion;
        }


    }//end mc reading script
    /*
    	for (int i = 0; i < num_runs; i++)

    	{
    // // 		dirc_model->set_pmt_angle(spread_ang.Gaus(47.87,box_rot_unc));
    		dirc_model->set_focus_mirror_angle(\
    			spread_ang.Gaus(74.11,mirror_angle_change_unc),\
    			spread_ang.Gaus(0,mirror_angle_change_yunc));
    		dirc_model->set_upper_wedge_angle_diff(\
    			spread_ang.Gaus(0,wedge_uncertainty),\
    			spread_ang.Gaus(0,upper_wedge_yang_spread));
    // 		dirc_model->set_bar_box_angle(spread_ang.Gaus(0,bar_box_box_angle));

    		printf("\r                                                    ");
    		printf("\rrunning iter %d/%d  ",i+1,num_runs);

    		if (overlap_x > 0)
    		{
    			confound_points = dirc_model->sim_rand_n_photons(\
    			      n_sim_phots,\
    			      pion_angle,\
    			      particle_x,\
    			      particle_y,\
    			      particle_theta,\
    			      particle_phi,\
    			      tracking_unc,\
    			      ckov_unc,\
    			      1);
    			//assume electrom has v=c

    			for (unsigned int i = 0; i < confound_points.size(); i++)
    			{
    			      confound_points[i].x += overlap_x;
    			}
    		}


    		fflush(stdout);
    		sim_points = dirc_model->sim_rand_n_photons(\
    			n_sim_phots,\
    			pion_angle ,
    			particle_x,\
    			particle_y,\
    			particle_theta,\
    			particle_phi,\
    			tracking_unc,\
    			ckov_unc,\
    			pion_beta);
    		for (unsigned int i = 0; i < confound_points.size(); i++)
    		{
    			sim_points.push_back(confound_points[i]);
    		}//
    		digitizer.digitize_points(sim_points);

    		llc = pdf_pion->get_log_likelihood(sim_points);
    		llf = pdf_kaon->get_log_likelihood(sim_points);


    		ll_diff_pion->Fill(llc-llf);

    		phot_found_pion->Fill(sim_points.size());

    		sim_points = dirc_model->sim_rand_n_photons(\
    			n_sim_phots,\
    			kaon_angle,\
    			particle_x,\
    			particle_y,\
    			particle_theta,\
    			particle_phi,\
    			tracking_unc,\
    			ckov_unc,\
    			kaon_beta);

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
    				pion_angle ,
    				particle_x,\
    				particle_y,\
    				particle_theta,\
    				particle_phi,\
    				tracking_unc,\
    				ckov_unc,\
    				pion_beta);

    			for (unsigned int j = 0; j < sim_points.size(); j++)
    			{
    				pion_coverage_xy->Fill(sim_points[j].x,sim_points[j].y);
    			}

    			sim_points = dirc_model->sim_rand_n_photons(\
    				n_sim_phots,\
    				pion_angle ,
    				particle_x,\
    				particle_y,\
    				particle_theta,\
    				particle_phi,\
    				tracking_unc,\
    				ckov_unc,\
    				kaon_beta);

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

    /*	if (out_csv == true)
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
    */

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
    phot_found_pion->Write();
    phot_found_kaon->Write();
    pion_cerenkov->Write();
    kaon_cerenkov->Write();
    pion_coverage_xy->Write();
    kaon_coverage_xy->Write();
    liquid_dist->Write();
    pion_lambda->Write();
    kaon_lambda->Write();


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
