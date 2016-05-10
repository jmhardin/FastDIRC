#include <vector>

#include "../include/dirc_base_sim.h"
#include "../include/dirc_point.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
//
#include <stdlib.h>
#include <math.h>

#include <TRandom3.h>



DircBaseSim::DircBaseSim(
		int rand_seed /*=4357*/,\
		double ibarLength/*=4900*/, \
		double ibarWidth/*=35*/, \
		double ibarDepth/*=17.25*/,
		double iupperWedgeTop/*=178.6*/) {

	barLength=ibarLength;
	barWidth=ibarWidth;
	barDepth=ibarDepth;
	wedgeWidthOff = 1.75;
	wedgeDepthOff = 10;
	wedgeFarAngle = .006*57.3;
	wedgeCloseAngle = 30;
	wedgeWidth=barWidth - wedgeWidthOff;
	wedgeDepthHigh = 79;
	wedgeHeight = 91;

	upperWedgeAngleStore = false;

	wedgeDepthHigh = wedgeDepthOff+barDepth+wedgeHeight*sin(wedgeCloseAngle/57.296);

	windowThickness = 9.6;


	upperWedgeTop = iupperWedgeTop;
	upperWedgeBottom = wedgeHeight + windowThickness;
	upperWedgeHeight = upperWedgeTop - upperWedgeBottom;
	upperWedgeDepthHigh = wedgeDepthHigh + (upperWedgeTop-wedgeHeight)*sin(wedgeCloseAngle/57.296);

	lowerWedgeExtensionZ = -wedgeDepthHigh\
			       - tan(wedgeCloseAngle/57.296)*(upperWedgeBottom - wedgeHeight);


	//Variables used for plane intersection
	wedgeClosePlaneNx = 0; //Shouldn't be needed
	wedgeClosePlaneNy = sin(wedgeCloseAngle/57.296);
	wedgeClosePlaneNz = cos(wedgeCloseAngle/57.296);

	upperWedgeClosePlaneNx = 0; //Shouldn't be needed
	upperWedgeClosePlaneNy = sin(wedgeCloseAngle/57.296);
	upperWedgeClosePlaneNz = cos(wedgeCloseAngle/57.296);

	upperWedgeNonUniform = false;
	upperWedgeNonUniformSpread = 0;

	wedgeClosePlaneD = barLength/2*wedgeClosePlaneNy - wedgeClosePlaneNz * (barDepth+wedgeDepthOff);

	//	upperWedgeClosePlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeClosePlaneNy + upperWedgeClosePlaneNz*lowerWedgeExtensionZ;

	//upperWedgeClosePlaneD = wedgeClosePlaneD; //should be in the same plane;

	upperWedgeFarPlaneNx = 0; //Shouldn't be needed
	upperWedgeFarPlaneNy = 0;
	upperWedgeFarPlaneNz = -1;

	upperWedgeFarPlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeFarPlaneNy;

	//Take average
	quartzIndex = 1.47;
	liquidIndex = 1.47;
	quartzLiquidY = upperWedgeBottom + 15;

	//Negative to make reflections easier
	wedgeFarPlaneNx = 0; //Shouldn't be needed
	wedgeFarPlaneNy = -sin(wedgeFarAngle/57.296);
	wedgeFarPlaneNz = -cos(wedgeFarAngle/57.296);


	wedgeFarPlaneD = barLength/2*wedgeFarPlaneNy;

	upperWedgeFarZ = 0;//Change based on geometry

	rand_gen = new TRandom3(rand_seed);

	box_angle_off_cval = 1;
	box_angle_off_sval = 0;


	num_transmittance = 36;
	min_transmittance = 300;
	max_transmittance = 660;
	sep_transmittance = (max_transmittance - min_transmittance)/(num_transmittance - 1);

	double tmp_quartz_transmittance[36] = {\
		0.999572036,0.999544661,0.999515062,0.999483019,0.999448285,\
			0.999410586,0.999369611,0.999325013,0.999276402,0.999223336,\
			0.999165317,0.999101778,0.999032079,0.998955488,0.998871172,\
			0.998778177,0.99867541 ,0.998561611,0.998435332,0.998294892,\
			0.998138345,0.997963425,0.997767484,0.997547418,0.99729958 ,\
			0.99701966 ,0.99670255 ,0.996342167,0.995931242,0.995461041,\
			0.994921022,0.994298396,0.993577567,0.992739402,0.991760297,\
			0.990610945\
	};
	for (int i = 0; i < num_transmittance; i++) {
		quartz_transmittance.push_back(tmp_quartz_transmittance[i]);
	}

	midLineMode = false;
	midLineWedgeWallFlip = 1;

	build_system();
}
void DircBaseSim::set_kaleidoscope_plot(bool ikp) {
	kaleidoscope_plot = ikp;
}
void DircBaseSim::rotate_2d(double &x, double &y, double cval, double sval) {
	//Standard rotatation allows precomputation of matrix elements
	//Only store one variable for speed
	//Sin should be actual sin (not negative sin)
	double tx = x;
	x = cval*tx - sval*y;
	y = sval*tx + cval*y;
}

void DircBaseSim::set_bar_box_angle(double ang) {
	//expect degrees
	//sets angle between readout box and bars - rotates angle coming out of bars
	box_angle_off_cval = cos(ang/57.3);
	box_angle_off_sval = sin(ang/57.3);
}
std::vector<double> DircBaseSim::get_dist_traveled() {
	return dist_traveled;
}
void DircBaseSim::set_liquid_index(double li) {
	//Sets the index of the upper quartz wedge
	liquidIndex = li;
	//printf("Liquid index set: %12.04f\n",liquidIndex);

}
void DircBaseSim::set_upper_wedge_angle_diff(double rads, double rads_y) {
	upperWedgeClosePlaneNx = 0; //Shouldn't be needed
	upperWedgeClosePlaneNy = sin(wedgeCloseAngle/57.296 + rads);
	upperWedgeClosePlaneNz = cos(wedgeCloseAngle/57.296 + rads);

	rotate_2d(upperWedgeClosePlaneNx,upperWedgeClosePlaneNz,cos(rads_y),sin(rads_y));

	upperWedgeClosePlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeClosePlaneNy + upperWedgeClosePlaneNz*lowerWedgeExtensionZ;

	upperWedgeFarPlaneNx = 0; //Shouldn't be needed
	upperWedgeFarPlaneNy = -sin(rads);
	upperWedgeFarPlaneNz = -cos(rads);

	upperWedgeFarPlaneD = upperWedgeBottom*upperWedgeFarPlaneNy;
}
void DircBaseSim::set_wedge_mirror_rand(double ispread) {
	if (ispread > 0) {
		upperWedgeNonUniform = true;
		upperWedgeNonUniformSpread = ispread;
	} else {
		upperWedgeNonUniform = false;
		upperWedgeNonUniformSpread = 0;
	}
	spread_wedge_mirror();
}
void DircBaseSim::spread_wedge_mirror() {

	if (upperWedgeNonUniform == true) {
		double spread_ang = wedgeCloseAngle;
		spread_ang += rand_gen->Gaus(0,upperWedgeNonUniformSpread);
		upperWedgeClosePlaneNx = 0; //Shouldn't be needed
		upperWedgeClosePlaneNy = sin(spread_ang/57.296);
		upperWedgeClosePlaneNz = cos(spread_ang/57.296);
		upperWedgeClosePlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeClosePlaneNy + upperWedgeClosePlaneNz*lowerWedgeExtensionZ;
	}
}
void DircBaseSim::build_system() {
	spread_wedge_mirror();
}
double DircBaseSim::get_quartz_n(double lambda) {
	double B1,B2,B3,C1,C2,C3;
	B1 = 0.6961663;             // B1
	B2 = 0.4079426;             // B2
	B3 = 0.8974794;             // B3
	C1 = 0.0046791;             // C1
	C2 = 0.0135121;             // C2
	C3 = 97.9340025;          // C3
	double lam2;
	double n_lam;
	lam2 = lambda*lambda/1000000;

	n_lam = sqrt(1 + B1*lam2/(lam2-C1) + B2*lam2/(lam2-C2) + B3*lam2/(lam2-C3));

	return n_lam;
}
bool DircBaseSim::quartz_transmission_mc(double R, double lambda) {
	int ind_transmittance;
	double above_ind;
	double tmp_transmittance;
	double tmp_transmittance_constant;
	double trans_prob;

	ind_transmittance = (lambda - min_transmittance)/sep_transmittance;

	above_ind = lambda - (min_transmittance + sep_transmittance*ind_transmittance);

	//Simple linear interp between values.  5th order poly fit looked like it worked too
	tmp_transmittance = quartz_transmittance[ind_transmittance]*(sep_transmittance-above_ind)/sep_transmittance\
			    + quartz_transmittance[ind_transmittance+1]*above_ind/sep_transmittance;

	tmp_transmittance_constant = log(tmp_transmittance);



	trans_prob = exp(R/1000.*tmp_transmittance_constant);


	return (rand_gen->Uniform(0,1) < trans_prob);

}
double DircBaseSim::get_beta(double E, double m) {
	double gam = E/m;
	if (gam > 5) {
		//Approximate and possibly save time
		return 1-1/(2*gam*gam);
	} else {
		return sqrt(1-1/(gam*gam));
	}

}
double DircBaseSim::get_bar_offset(int bar)
{
	return fabs(bar)/bar*((150-0.5*(barWidth))+bar*(barWidth+.015));
}
int DircBaseSim::get_bar_from_x(double x)
{

	if (x < -150)
	{
		return (x+150)/(barWidth+.015) - 1;
	}
	else if (x > 150)
	{
		return  (x-150)/(barWidth+.015) + 1;
	}
	else
	{
		return 0;
	}

}
void DircBaseSim::set_upper_wedge_angle_store(bool istore)
{
	upperWedgeAngleStore = istore;
}
std::vector<double> DircBaseSim::get_upper_wedge_incident()
{
	return upper_wedge_incident;
}
void DircBaseSim::sim_rand_n_photons(\
		std::vector<dirc_point> &out_points,\
		int n_photons, \
		double ckov_theta /*= 47*/, \
		double particle_bar /*=0*/, \
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_t /*=0*/,\
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc /*= .0015*57.3*/,\
		double ckov_theta_unc /* = .0055*57.3*/,\
		double beta /* = -1*/) {

	out_points.clear();
	fill_rand_phi(\
			out_points,\
			n_photons,\
			ckov_theta,\
			particle_bar,\
			particle_x,\
			particle_y,\
			particle_t /*=0*/,\
			particle_theta,\
			particle_phi,\
			phi_theta_unc,\
			ckov_theta_unc,\
			beta);
}
void DircBaseSim::sim_reg_n_photons(\
		std::vector<dirc_point> &out_points,\
		int n_photons_phi, \
		int n_photons_z,\
		double ckov_theta /*= 47*/, \
		double particle_bar /*=0*/, \
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_t /*=0*/,\
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc /*= 0*/,\
		double ckov_theta_unc /* = 0*/,\
		double beta /* = -1*/) {
	//     std::vector<dirc_point> out_points;
	out_points.clear();
	fill_reg_phi(\
			out_points,\
			n_photons_phi,\
			n_photons_z,\
			ckov_theta,\
			particle_bar,\
			particle_x,\
			particle_y,\
			particle_t,\
			particle_theta,\
			particle_phi,\
			phi_theta_unc,\
			ckov_theta_unc,\
			beta);
	//Reflection inside loops
	//    sidemirror_reflect_points(out_points);
	//     return out_points;
}
void DircBaseSim::test_from_wedge_top(\
		std::vector<dirc_point> &ovals,\
		int n_photons, \
		double particle_bar /*= 1*/, \
		double particle_x /*= 0*/, \
		double phot_theta /*= 0*/, \
		double phot_phi /*= 0*/,\
		double theta_unc, /*= 0*/
		double phi_unc /* = 0*/,\
		double overall_theta) {

	ovals.clear();
	//Note that Theta and Phi are defined along the bar axis, not as they are elsewhere
	//negative bar number is - x axis?
	int numPhots = n_photons;

	double sourcez = -barDepth/2;
	sourcez = -wedgeDepthHigh/2;
	sourcez = 0;
	double sourcey = barLength/2 + upperWedgeTop + .1;
	double sourcex = particle_x;

	double temit, randPhi;

	double x,y,z,dx,dy,dz;

	double mm_index = 0;

	for (int i = 0; i < numPhots; i++) {
		mm_index = 0;
		randPhi = phot_phi + rand_gen->Gaus(0,phi_unc);
		temit = phot_theta + rand_gen->Gaus(0,theta_unc);

		x = sourcex;
		y = sourcey;
		z = sourcez;
		//Note, the geometry change again
		dx = sin(temit/57.3)*sin(randPhi/57.3);
		dy = cos(temit/57.3);
		dz = sin(temit/57.3)*cos(randPhi/57.3);

		rotate_2d(dy,dz,cos(overall_theta/57.3),sin(overall_theta/57.3));

		optical_interface_z(quartzIndex,liquidIndex,dx,dz,dy);

		dirc_point out_val;
		warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);

		if (out_val.t < 0)
		{
			continue;
		}	

		ovals.push_back(out_val);
	}
}
void DircBaseSim::sim_lut_points(\
		std::vector<dirc_point> &ovals,\
		std::vector<double> &phis,\
		std::vector<double> &thetas,\
		int n_photons, \
		double particle_bar /*= 0*/){
	
	ovals.clear();
	phis.clear();
	thetas.clear();
	

	double x,y,z,dx,dy,dz;

	double randPhi;
	double randTheta;	

	double mm_index = 0;

	double c_mm_ns = 300;

	for (int i = 0; i < n_photons; i++) {
		randPhi = rand_gen->Uniform(0,2*3.14159265);
		randTheta = acos(2*rand_gen->Uniform(.5,1) - 1);
		//This timing won't be correct
		mm_index = 0;

		x = 0;
		y = barLength/2;
		z = -barDepth/2;

		dx = sin(randTheta)*sin(randPhi);
		dy = cos(randTheta);
		dz = sin(randTheta)*cos(randPhi); 

		mm_index += warp_wedge(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);
		//account (quickly) for the bar box having a different angle than the readout
		rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);

		if (z > 0) {
			continue;
		}
		dirc_point out_val;


		warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);

		if (out_val.t < 0)
		{
			continue;
		}

		//Timing is hard...
		//out_val.t = mm_index/c_mm_ns;
		out_val.updown = 0;
		ovals.push_back(out_val);
		phis.push_back(randPhi);
		thetas.push_back(randTheta);
	}
}
void DircBaseSim::fill_rand_phi(\
		std::vector<dirc_point> &ovals,\
		int n_photons, \
		double ckov_theta /*= 47*/, \
		double particle_bar /*= 0*/, \
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_t /*= 0*/, \
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc, /*= .0015*57.3*/
		double ckov_theta_unc /* = .0055*57.3*/,\
		double beta/* = -1*/) {


	double emitAngle = ckov_theta;
	double particleTheta = particle_theta + rand_gen->Gaus(0,phi_theta_unc);
	double particlePhi = particle_phi + rand_gen->Gaus(0,phi_theta_unc);
	int numPhots = n_photons/cos(particle_theta/57.3);

	double sourceOff,randPhi;

	double temit, rand_add;
	double wavelength = 400;

	double x,y,z,dx,dy,dz;

	double cos_ptheta = cos(particleTheta/57.3);
	double sin_ptheta = sin(particleTheta/57.3);
	double cos_pphi = cos(particlePhi/57.3);
	double sin_pphi = sin(particlePhi/57.3);

	double mm_index = 0;

	int tmp_updown = 0;

	for (int i = 0; i < numPhots; i++) {
		randPhi = rand_gen->Uniform(0,2*3.14159265);
		sourceOff = -rand_gen->Uniform(0,barDepth);
		if (kaleidoscope_plot == true)
		{	
			sourceOff = -barDepth/2;
		}

		if (beta < 0) {
			rand_add = rand_gen->Gaus(0,ckov_theta_unc);
			temit = emitAngle + rand_add;
		} else {
			temit = get_cerenkov_angle_rand(beta,ckov_theta_unc,wavelength);
		}
		mm_index = (sourceOff - barDepth)*1.47/cos(particleTheta/57.3);

		x = 0;
		y = 0;
		z = sourceOff;

		dx = sin(temit/57.3)*cos(randPhi);
		dy = sin(temit/57.3)*sin(randPhi);
		dz = cos(temit/57.3);

		rotate_2d(z,y,cos_ptheta,sin_ptheta);
		rotate_2d(x,y,cos_pphi,sin_pphi);

		rotate_2d(dz,dy,cos_ptheta,sin_ptheta);
		rotate_2d(dy,dx,cos_pphi,sin_pphi);

		z -= barDepth;
		x += particle_x;
		y += particle_y;



		//photon is now defined as up or down
		if (dy > 0)
		{
			tmp_updown = 1;
		}
		else
		{//Assume that zero doesn't count
			tmp_updown = -1;
		}

		mm_index += warp_ray(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz,\
				sqrt(1-1/(1.47*1.47)));

		if (z > 0)
		{
			continue;
		}

		spread_wedge_mirror();

		mm_index += warp_wedge(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);
		//account (quickly) for the bar box having a different angle than the readout
		rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);

		if (z > 0) {
			continue;
		}
		dirc_point out_val;


		warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);

		if (out_val.t < 0)
		{
			continue;
		}

		out_val.t += particle_t;
		out_val.updown = tmp_updown;
		ovals.push_back(out_val);
	}
}
std::vector<std::pair<double,double> > DircBaseSim::get_refraction_rand_phi(\
		std::vector<double> &before_interface,\
		std::vector<double> &after_interface,\
		std::vector<double> &pmt_incidence,\
		int n_photons, \
		double ckov_theta /*= 47*/, \
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc /*= .0015*57.3*/,\
		double ckov_theta_unc /* = .0055*57.3*/,\
		double beta/* = -1*/) {
	//returns theta versus cerenkov phi_theta_unc

	std::vector<std::pair<double,double> > rval;

	before_interface.clear();
	after_interface.clear();
	refraction_before.clear();
	refraction_after.clear();
	pmt_incidence.clear();
	store_refraction = true;

	double emitAngle = ckov_theta;
	double particleTheta = particle_theta + rand_gen->Gaus(0,phi_theta_unc);
	double particlePhi = particle_phi + rand_gen->Gaus(0,phi_theta_unc);

	int numPhots = n_photons/cos(particle_theta/57.3);

	// 	double sourcez = -sDepth;
	double sourcey = particle_y-barDepth*tan(particleTheta/57.3);
	double sourcex = particle_x;
	double tsy = sourcey;
	double tsx = sourcex;
	sourcey = tsy*cos(particlePhi/57.3)-tsx*sin(particlePhi/57.3);
	sourcex = tsy*sin(particlePhi/57.3)+tsx*cos(particlePhi/57.3);


	double sourceOff,randPhi;

	double temit, rand_add;
	double wavelength = 400;

	double x,y,z,dx,dy,dz;

	double cos_ptheta = cos(particle_theta/57.3);
	double sin_ptheta = sin(particle_theta/57.3);
	double cos_pphi = cos(particle_phi/57.3);
	double sin_pphi = sin(particle_phi/57.3);

	double mm_index = 0;
	double c_mm_ns = 300;

	for (int i = 0; i < numPhots; i++) {
		randPhi = rand_gen->Uniform(0,2*3.14159265);
		sourceOff = -rand_gen->Uniform(0,barDepth);

		if (beta < 0) {
			rand_add = rand_gen->Gaus(0,ckov_theta_unc);
			temit = emitAngle + rand_add;
		} else {
			temit = get_cerenkov_angle_rand(beta,ckov_theta_unc,wavelength);
		}
		mm_index = (sourceOff - barDepth)*1.47;

		x = 0;
		y = 0;
		z = sourceOff;

		dx = sin(temit/57.3)*cos(randPhi);
		dy = sin(temit/57.3)*sin(randPhi);
		dz = cos(temit/57.3);

		rotate_2d(z,y,cos_ptheta,sin_ptheta);
		rotate_2d(x,y,cos_pphi,sin_pphi);

		rotate_2d(dz,dy,cos_ptheta,sin_ptheta);
		rotate_2d(dy,dx,cos_pphi,sin_pphi);

		z -= barDepth;
		x += particle_x;
		y += particle_y;


		mm_index += warp_ray(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz,\
				asin(1/1.47));

		if (z > 0) {
			continue;
		}


		spread_wedge_mirror();

		mm_index += warp_wedge(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);

		//account (quickly) for the bar box having a different angle than the readout
		rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);

		if (z > 0) {
			continue;
		}
		dirc_point out_val;

		warp_readout_box(\
				out_val,\
				0,\
				mm_index,\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);
		if (out_val.t < 0)
		{
			continue;
		}

		//should be threading time information into this soon
		out_val.t = mm_index/(c_mm_ns);
		std::pair<double, double> add_val;

		add_val.first = randPhi;
		add_val.second = refraction_before[refraction_before.size() - 1];
		rval.push_back(add_val);
	}

	for (unsigned int i = 0; i < refraction_before.size(); i++) {
		before_interface.push_back(refraction_before[i]);
		after_interface.push_back(refraction_after[i]);
	}

	store_refraction = false;
	refraction_before.clear();
	refraction_after.clear();

	return rval;
}

void DircBaseSim::fill_reg_phi(\
		std::vector<dirc_point> &ovals,\
		int n_photons_phi, \
		int n_photons_z,\
		double ckov_theta /*= 47*/, \
		double particle_bar /*= 0*/,\
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_t /*= 0*/, \
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc /*= 0*/,\
		double ckov_theta_unc /* = 0*/,\
		double beta /* = -1*/)
{
	double sDepth = .95*barDepth;
	double emitAngle = ckov_theta;

	int tmp_updown = 0;

	double sourceOff,regPhi;

	double temit;
	double rand_add;
	double wavelength = 400;

	double x,y,z,dx,dy,dz;

	double cos_ptheta = cos(particle_theta/57.3);
	double sin_ptheta = sin(particle_theta/57.3);
	double cos_pphi = cos(particle_phi/57.3);
	double sin_pphi = sin(particle_phi/57.3);

	double mm_index = 0;

	double sin_emit;
	double cos_emit;
	double sin_regphi;
	double cos_regphi;

	int adj_n_photons_phi = n_photons_phi/cos(particle_theta/57.3);

	for (int i = 0; i < n_photons_z; i++) {
		sourceOff = (i+.5)*sDepth/(n_photons_z);

		if (kaleidoscope_plot == true)
		{
			sourceOff = -sDepth/2;
		}

		for (int j = 0; j < adj_n_photons_phi; j++) {
			regPhi = j*2*3.14159265357/(adj_n_photons_phi);

			if (beta < 0) {
				rand_add = rand_gen->Gaus(0,ckov_theta_unc);
				temit = emitAngle + rand_add;
			} else {
				temit = get_cerenkov_angle_rand(beta,ckov_theta_unc,wavelength);
			}

			mm_index = (sourceOff - barDepth)*1.47;

			x = 0;
			y = 0;
			z = sourceOff;

			//save some time ~30ms per 400k
			//could compute sin even faster with a taylor series
			sin_emit = sin(temit/57.2957795131);
			cos_emit = sqrt(1-sin_emit*sin_emit);
			cos_regphi = cos(regPhi);
			sin_regphi = sgn(3.14159265359 - regPhi)*sqrt(1-cos_regphi*cos_regphi);

			dx = sin_emit*cos_regphi;
			dy = sin_emit*sin_regphi;
			dz = cos_emit;

			rotate_2d(z,y,cos_ptheta,sin_ptheta);
			rotate_2d(x,y,cos_pphi,sin_pphi);

			rotate_2d(dz,dy,cos_ptheta,sin_ptheta);
			rotate_2d(dy,dx,cos_pphi,sin_pphi);

			z -= barDepth;
			x += particle_x;
			y += particle_y;

			//photon is now defined as up or down
			if (dy > 0)
			{
				tmp_updown = 1;
			}
			else
			{//Assume that zero doesn't count
				tmp_updown = -1;
			}

			mm_index += warp_ray(\
					x,\
					y,\
					z,\
					dx,\
					dy,\
					dz,\
					sqrt(1-1/(1.47*1.47)));

			if (z > 0)
			{
				continue;
			}

			spread_wedge_mirror();

			mm_index += warp_wedge(\
					x,\
					y,\
					z,\
					dx,\
					dy,\
					dz);

			//account (quickly) for the bar box having a different angle than the readout
			rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);

			if (z > 0)
			{
				continue;
			}

			dirc_point out_val;
			warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);

			if (out_val.t < 0)
			{
				continue;
			}

			out_val.t += particle_t;
			out_val.updown = tmp_updown;
			ovals.push_back(out_val);
		}
	}

}
bool DircBaseSim::track_single_photon(\
		dirc_point &out_val,\
		double emit_theta,\
		double emit_phi,\
		double particle_theta,\
		double particle_phi,\
		double particle_x,\
		double particle_y,\
		double particle_z,\
		double particle_t,\
		int particle_bar)
{
	double x,y,z,dx,dy,dz;
	double mm_index;

	double temit = emit_theta;
	double regPhi = emit_phi;

	double cos_emit, sin_emit;
	double cos_regphi, sin_regphi;
	int tmp_updown = 0;

	double cos_ptheta = cos(particle_theta/57.3);
	double sin_ptheta = sin(particle_theta/57.3);
	double cos_pphi = cos(particle_phi/57.3);
	double sin_pphi = sin(particle_phi/57.3);

	double sourceOff = particle_z + barDepth;

	mm_index = (sourceOff - barDepth)*1.47;

	x = 0;
	y = 0;
	z = sourceOff;

	//save some time ~30ms per 400k
	//could compute sin even faster with a taylor series
	sin_emit = sin(temit/57.2957795131);
	cos_emit = sqrt(1-sin_emit*sin_emit);
	cos_regphi = cos(regPhi);
	sin_regphi = sgn(3.14159265359 - regPhi)*sqrt(1-cos_regphi*cos_regphi);

	dx = sin_emit*cos_regphi;
	dy = sin_emit*sin_regphi;
	dz = cos_emit;

	rotate_2d(z,y,cos_ptheta,sin_ptheta);
	rotate_2d(x,y,cos_pphi,sin_pphi);

	rotate_2d(dz,dy,cos_ptheta,sin_ptheta);
	rotate_2d(dy,dx,cos_pphi,sin_pphi);

	z -= barDepth;
	x += particle_x;
	y += particle_y;

	//photon is now defined as up or down
	if (dy > 0)
	{
		tmp_updown = 1;
	}
	else
	{//Assume that zero doesn't count
		tmp_updown = -1;
	}

	mm_index += warp_ray(\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz,\
			sqrt(1-1/(1.47*1.47)));


	if (z > 0)
	{
		return false;
	}

	spread_wedge_mirror();


	if (dx < 0)
	{
		lastWallX = 1;
	}
	else
	{
		lastWallX = -1;
	}
	int beforeWedgeLastWall = lastWallX;

	mm_index += warp_wedge(\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);

	if (dx < 0)	
	{
		lastWallX = 1;
	}
	else
	{
		lastWallX = -1;
	}
	out_val.last_wall_x = lastWallX;
	//	printf("lastwall agrees: %d\n",beforeWedgeLastWall==lastWallX);
	//account (quickly) for the bar box having a different angle than the readout
	rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);

	if (z > 0)
	{
		return false;
	}

	out_val.wedge_before_interface = -1;
	warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);

	out_val.wedge_before_interface = wedgeBeforeInterface;

	if (out_val.t < 0)
	{
		return false;
	}

	out_val.t += particle_t;
	out_val.updown = tmp_updown;

	return true;
}
bool DircBaseSim::track_line_photon(\
		dirc_point &out_val,\
		double emit_theta,\
		double emit_phi,\
		double particle_theta,\
		double particle_phi,\
		double particle_x,\
		double particle_y,\
		double particle_z,\
		double particle_t,\
		int particle_bar,\
		double z_at_top /*=1*/)
{
	//I'm being bad and breaking encapsulation here, but it's for the greater good
	//Think of the children
	double saveBarWidth = barWidth;
	double saveWedgeWidthOff = wedgeWidthOff;

	midLineMode = true;
	midLineMode = false;
	double width_fraction = 1;
	barWidth *= width_fraction;
	wedgeWidthOff *= width_fraction;
	wedgeWidthOff = 0;
	build_system();

	double x,y,z,dx,dy,dz;
	double mm_index;

	double temit = emit_theta;
	double regPhi = emit_phi;

	double cos_emit, sin_emit;
	double cos_regphi, sin_regphi;
	int tmp_updown = 0;

	double cos_ptheta = cos(particle_theta/57.3);
	double sin_ptheta = sin(particle_theta/57.3);
	double cos_pphi = cos(particle_phi/57.3);
	double sin_pphi = sin(particle_phi/57.3);

	double sourceOff = particle_z + barDepth;

	mm_index = (sourceOff - barDepth)*1.47;

	x = 0;
	y = 0;
	z = sourceOff;

	//save some time ~30ms per 400k
	//could compute sin even faster with a taylor series
	sin_emit = sin(temit/57.2957795131);
	cos_emit = sqrt(1-sin_emit*sin_emit);
	cos_regphi = cos(regPhi);
	sin_regphi = sgn(3.14159265359 - regPhi)*sqrt(1-cos_regphi*cos_regphi);

	dx = sin_emit*cos_regphi;
	dy = sin_emit*sin_regphi;
	dz = cos_emit;

	rotate_2d(z,y,cos_ptheta,sin_ptheta);
	rotate_2d(x,y,cos_pphi,sin_pphi);

	rotate_2d(dz,dy,cos_ptheta,sin_ptheta);
	rotate_2d(dy,dx,cos_pphi,sin_pphi);

	z -= barDepth;
	x += particle_x;
	y += particle_y;

	//photon is now defined as up or down
	if (dy > 0)
	{
		tmp_updown = 1;
	}
	else
	{//Assume that zero doesn't count
		tmp_updown = -1;
	}

	mm_index += warp_ray(\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz,\
			sqrt(1-1/(1.47*1.47)));


	if (z > 0)
	{
		barWidth = saveBarWidth;
		wedgeWidthOff = saveWedgeWidthOff;
		build_system();
		midLineMode = false;
		return false;
	}

	spread_wedge_mirror();


	if (dx < 0)
	{
		lastWallX = 1;
	}
	else
	{
		lastWallX = -1;
	}
	int beforeWedgeLastWall = lastWallX;


	double target_z = -17.25/2;

	if (z_at_top < 0)
	{
		target_z = z_at_top;
	}

	//z = target_z;
	if (dz > 0)
	{
		//y += dy/dz*(-z);
	//	dz = -dz;
	}


	mm_index += warp_wedge(\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);

	if (dx < 0)	
	{
		lastWallX = 1;
	}
	else
	{
		lastWallX = -1;
	}
	out_val.last_wall_x = lastWallX;
	//	printf("lastwall agrees: %d\n",beforeWedgeLastWall==lastWallX);
	//account (quickly) for the bar box having a different angle than the readout
	rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);

	if (z > 0)
	{
		barWidth = saveBarWidth;
		wedgeWidthOff = saveWedgeWidthOff;
		build_system();
		midLineMode = false;
		return false;
	}
	if (lastWallX == 1)
	{
	//	x += saveBarWidth/2 - 3;
	}
	else if (lastWallX == -1)
	{
		//x += saveBarWidth/2 - saveWedgeWidthOff;
	//	x += saveBarWidth/2 - 3;
	}	
	out_val.wedge_before_interface = -1;
	warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);

	out_val.wedge_before_interface = wedgeBeforeInterface;

	if (out_val.t < 0)
	{
		barWidth = saveBarWidth;
		wedgeWidthOff = saveWedgeWidthOff;
		build_system();
		midLineMode = false;
		return false;
	}

	out_val.t += particle_t;
	out_val.updown = tmp_updown;

	barWidth = saveBarWidth;
	wedgeWidthOff = saveWedgeWidthOff;
	build_system();
	midLineMode = false;
	return true;
}
bool DircBaseSim::track_single_photon_beta(\
		dirc_point &out_val,\
		double particle_beta,\
		double emit_phi,\
		double particle_theta,\
		double particle_phi,\
		double particle_x,\
		double particle_y,\
		double particle_z,\
		double particle_t,\
		int particle_bar)
{
	double x,y,z,dx,dy,dz;
	double mm_index;

	double wavelength = 400;

	double temit = get_cerenkov_angle_rand(particle_beta,0,wavelength);
	double regPhi = emit_phi;

	double cos_emit, sin_emit;
	double cos_regphi, sin_regphi;
	int tmp_updown = 0;

	double cos_ptheta = cos(particle_theta/57.3);
	double sin_ptheta = sin(particle_theta/57.3);
	double cos_pphi = cos(particle_phi/57.3);
	double sin_pphi = sin(particle_phi/57.3);

	double sourceOff = particle_z + barDepth;

	mm_index = (sourceOff - barDepth)*1.47;

	x = 0;
	y = 0;
	z = sourceOff;

	//save some time ~30ms per 400k
	//could compute sin even faster with a taylor series
	sin_emit = sin(temit/57.2957795131);
	cos_emit = sqrt(1-sin_emit*sin_emit);
	cos_regphi = cos(regPhi);
	sin_regphi = sgn(3.14159265359 - regPhi)*sqrt(1-cos_regphi*cos_regphi);

	dx = sin_emit*cos_regphi;
	dy = sin_emit*sin_regphi;
	dz = cos_emit;

	rotate_2d(z,y,cos_ptheta,sin_ptheta);
	rotate_2d(x,y,cos_pphi,sin_pphi);

	rotate_2d(dz,dy,cos_ptheta,sin_ptheta);
	rotate_2d(dy,dx,cos_pphi,sin_pphi);

	z -= barDepth;
	x += particle_x;
	y += particle_y;

	//photon is now defined as up or down
	if (dy > 0)
	{
		tmp_updown = 1;
	}
	else
	{//Assume that zero doesn't count
		tmp_updown = -1;
	}

	mm_index += warp_ray(\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz,\
			sqrt(1-1/(1.47*1.47)));


	if (z > 0)
	{
		return false;
	}

	spread_wedge_mirror();


	if (dx < 0)
	{
		lastWallX = 1;
	}
	else
	{
		lastWallX = -1;
	}
	int beforeWedgeLastWall = lastWallX;

	mm_index += warp_wedge(\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);

	if (dx < 0)	
	{
		lastWallX = 1;
	}
	else
	{
		lastWallX = -1;
	}
	out_val.last_wall_x = lastWallX;
	//	printf("lastwall agrees: %d\n",beforeWedgeLastWall==lastWallX);
	//account (quickly) for the bar box having a different angle than the readout
	rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);

	if (z > 0)
	{
		return false;
	}

	out_val.wedge_before_interface = -1;
	warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);

	out_val.wedge_before_interface = wedgeBeforeInterface;

	if (out_val.t < 0)
	{
		return false;
	}

	out_val.t += particle_t;
	out_val.updown = tmp_updown;

	return true;
}
bool DircBaseSim::track_all_line_photons(\
		std::vector<dirc_point> &left_vals,\
		std::vector<dirc_point> &right_vals,\
		int points_per_side,\
		double emit_theta,\
		double particle_theta,\
		double particle_phi,\
		double particle_x,\
		double particle_y,\
		double particle_z,\
		double particle_t,\
		int particle_bar,\
		double z_at_top /*=1*/)
{
	//I'm being bad and breaking encapsulation here, but it's for the greater good
	//Think of the children
	double saveBarWidth = barWidth;
	double saveWedgeWidthOff = wedgeWidthOff;

	left_vals.clear();
	right_vals.clear();

	midLineMode = true;
	double width_fraction = .001;
	barWidth *= width_fraction;
	wedgeWidthOff *= width_fraction;
	wedgeWidthOff = 0;
	build_system();

	std::vector<double> mustache_phi;

	double x,y,z,dx,dy,dz;
	//righthand side versions
	double x_r,y_r,z_r,dx_r,dy_r,dz_r;
	double mm_index;
	double mm_index_r;

	double temit = emit_theta;
	double regPhi;

	double cos_emit, sin_emit;
	double cos_regphi, sin_regphi;
	int tmp_updown = 0;

	double cos_ptheta = cos(particle_theta/57.3);
	double sin_ptheta = sin(particle_theta/57.3);
	double cos_pphi = cos(particle_phi/57.3);
	double sin_pphi = sin(particle_phi/57.3);

	double sourceOff = particle_z + barDepth;

	double phi_inc = (2*3.14159265359)/points_per_side;
	regPhi = 0;
	for (int i = 0; i < points_per_side; i++)
	{
		x = 0;
		y = 0;
		z = sourceOff;
		mm_index = (sourceOff - barDepth)*1.47;
		regPhi += phi_inc;

		//save some time ~31ms per 400k
		//could compute sin even faster with a taylor series
		sin_emit = sin(temit/57.2957795131);
		cos_emit = sqrt(1-sin_emit*sin_emit);
		cos_regphi = cos(regPhi);
		sin_regphi = sgn(3.14159265359 - regPhi)*sqrt(1-cos_regphi*cos_regphi);

		dx = sin_emit*cos_regphi;
		dy = sin_emit*sin_regphi;
		dz = cos_emit;

		rotate_2d(z,y,cos_ptheta,sin_ptheta);
		rotate_2d(x,y,cos_pphi,sin_pphi);

		rotate_2d(dz,dy,cos_ptheta,sin_ptheta);
		rotate_2d(dy,dx,cos_pphi,sin_pphi);

		z -= barDepth;
		x += particle_x;
		y += particle_y;

		//photon is now defined as up or down
		if (dy > 0)
		{
			tmp_updown = 1;
		}
		else
		{//Assume that zero doesn't count
			tmp_updown = -1;
		}



		mm_index += warp_ray(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz,\
				sqrt(1-1/(1.47*1.47)));


		if (z > 0)
		{
			continue;
		}

		spread_wedge_mirror();


		int beforeWedgeLastWall = lastWallX;


		double target_z = -17.25/2;

		if (z_at_top < 0)
		{
			target_z = z_at_top;
		}

		z = target_z;
		if (dz > 0)
		{
			dz = -dz;
		}

		dirc_point out_val;

		//branch here
	
		if (upperWedgeBottom/(upperWedgeBottom*tan(wedgeCloseAngle/57.3)+1.5*barDepth+wedgeDepthOff) < -dy/dz)
		{
			mustache_phi.push_back(regPhi);
		//	printf("%12.04f\n",regPhi);
		}
	
		mm_index += warp_wedge(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);
		
		mm_index_r = mm_index;

		int dx_lr = sgn(dx);

		out_val.last_wall_x = -1;
		rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);

		x += saveBarWidth/2-3;
	
		//And this set of ifs is why we try and do one photon at a time normally
		if (z < 0)
		{
			out_val.wedge_before_interface = -1;
			warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);
			out_val.wedge_before_interface = wedgeBeforeInterface;
			if (out_val.t > 0)
			{
				out_val.t += particle_t;
				out_val.updown = tmp_updown;
				out_val.init_phi = regPhi;
				if (dx_lr < 0)
				{
					left_vals.push_back(out_val);
				}
				else
				{
					right_vals.push_back(out_val);
				}
			}
		}
	}
	for (unsigned int i = 0; i < mustache_phi.size(); i++)
	{
		x = 0;
		y = 0;
		z = sourceOff;
		mm_index = (sourceOff - barDepth)*1.47;
		regPhi = mustache_phi[i];
		

		//save some time ~31ms per 400k
		//could compute sin even faster with a taylor series
		sin_emit = sin(temit/57.2957795131);
		cos_emit = sqrt(1-sin_emit*sin_emit);
		cos_regphi = cos(regPhi);
		sin_regphi = sgn(3.14159265359 - regPhi)*sqrt(1-cos_regphi*cos_regphi);

		dx = sin_emit*cos_regphi;
		dy = sin_emit*sin_regphi;
		dz = cos_emit;

		rotate_2d(z,y,cos_ptheta,sin_ptheta);
		rotate_2d(x,y,cos_pphi,sin_pphi);

		rotate_2d(dz,dy,cos_ptheta,sin_ptheta);
		rotate_2d(dy,dx,cos_pphi,sin_pphi);

		z -= barDepth;
		x += particle_x;
		y += particle_y;

		//photon is now defined as up or down
		if (dy > 0)
		{
			tmp_updown = 1;
		}
		else
		{//Assume that zero doesn't count
			tmp_updown = -1;
		}



		mm_index += warp_ray(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz,\
				sqrt(1-1/(1.47*1.47)));


		if (z > 0)
		{
			continue;
		}

		spread_wedge_mirror();


		int beforeWedgeLastWall = lastWallX;

		//now in the mustache => at wall
		//And going up
		double target_z = 0;

		if (z_at_top < 0)
		{
			target_z = z_at_top;
		}

		z = target_z;
		if (dz > 0)
		{
			dz = -dz;
		}
		//Goes backwards first and reflects off of wall.  Could account for 6mrad here as well, consider later.
		y += -dy/dz*barDepth;
		//perform proper reflection - factor of 2 for the reflection, and an additional 1 for the corrected wedge.
		double enhance_rot_fac = 3;
		rotate_2d(dy,dz,cos(enhance_rot_fac*wedgeFarAngle/57.3),-sin(enhance_rot_fac*wedgeFarAngle/57.3));

		dirc_point out_val;

		//branch here
	
		mm_index += warp_wedge(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);
		
		mm_index_r = mm_index;

		int dx_lr = sgn(dx);

		out_val.last_wall_x = -1;
		rotate_2d(dy,dz,box_angle_off_cval,box_angle_off_sval);

		x += saveBarWidth/2-3;
	
		//And this set of ifs is why we try and do one photon at a time normally
		if (z < 0)
		{
			out_val.wedge_before_interface = -1;
			warp_readout_box(out_val,particle_bar,mm_index,x,y,z,dx,dy,dz);
			out_val.wedge_before_interface = wedgeBeforeInterface;
			if (out_val.t > 0)
			{
				out_val.t += particle_t;
				out_val.updown = tmp_updown;
				out_val.init_phi = regPhi;
				if (dx_lr < 0)
				{
					left_vals.push_back(out_val);
				}
				else
				{
					right_vals.push_back(out_val);
				}
			}
		}
	}
	barWidth = saveBarWidth;
	wedgeWidthOff = saveWedgeWidthOff;
	build_system();
	midLineMode = false;
	return true;
}
double DircBaseSim::warp_ray(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz,\
		double cos_critical_angle) {
	//Implemented to avoid total internal reflection computations
	//Expects ray to be "prerotated" by particle theta and phi - just propagates the vectors
	//Uses bar geometry.
	//Modifys the ray.
	//Be careful about x,y,and z - can't straight rotate that
	//returns distance traveled (mm) times index

	double rval = 0; //stores total distance;

	double delx, dely,delz;
	// 	double xzR;
	//Start the rays in the quartz for simulation reasons
	double grace_room = 0;

	if (dy > 0) {
		//Going up
		dely = barLength*(0.5-grace_room) - y;
	} else {
		//going down, flip dy
		dely = barLength*(1.5-grace_room) + y;
		dy = -dy;
	}

	rval += dely*dely;

	y = barLength*(.5-grace_room);

	if (fabs(dz) > cos_critical_angle ||\
			fabs(dx) > cos_critical_angle ||\
			(dy < 0 && fabs(dy) > cos_critical_angle) ||
			dy*dy < 1e-4) {
		//If it's not totally internally reflected, assume it escapes
		//also assume failure if it isn't going up fast enough (absorbed)
		//Positive z origin signals fail
		z = 1337;
		return -1;
	}

	//Del on x and z refers to distance from first reflection
	//I sincerly hope dz is positive - might be worth checking

	int nbouncesx;
	// 	int nbouncesy;
	int nbouncesz;

	double remainderx;
	// 	double remaindery;
	double lrx;
	// 	double lrz = 1;
	double remainderz = 0;
	/*Not used right now, so don't branch if you don't have to
	  if(dy > 0)
	  {
	  nbouncesy = 0;
	  }
	  else
	  {
	  nbouncesy = 1;
	  }*/

	//deterimines if particle is going left or right

	lrx = sgn(dx);

	delx = fabs(dely*dx/dy) - lrx*(barWidth/2-x);
	delz = fabs(dely*dz/dy) + z;

	rval += delx*delx + delz*delz;

	// 	if (dz < 0) printf("dz: %12.04f : Negative and confusing\n",dz);


	nbouncesx = delx/barWidth + 1;
	remainderx = delx - barWidth*(nbouncesx-1);

	if (nbouncesx % 2 == 1) {
		dx = -dx;
		x = lrx*(barWidth/2-remainderx);
	} else {
		x = lrx*(-barWidth/2+remainderx);
	}


	if (x > barWidth/2 - wedgeWidthOff) {
		//Hit the 5mm where the wedge is not.  Assume lost
		z = 1337;
		return -1;
	}
	nbouncesz = delz/barDepth + 1;
	remainderz = delz - barDepth*(nbouncesz-1);
	//	double bar_front = -barDepth/2;
	//	double bar_front = -17.25/2;
	double bar_front = 0;
	if (nbouncesz % 2 == 1) {
		dz = -dz;
		z = bar_front-remainderz;
	} else {
		z = remainderz-barDepth;
	}

	return sqrt(rval)*quartzIndex;
}
double DircBaseSim::warp_wedge(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz) {
	//No Critical angle checking - shouldn't need it
	//I am abusing the x angle being the same in the bar and wedge here
	//may have to rotate after to the mirror
	//Starts y at 0

	wedge_bounces = 0;

	double mm_index = 0;
	bool passed_interface = false;
	// 	double x0 = x;
	double dt = 0; //dummy variable representing parametric line
	double n_dot_v = 0; //dummy variable to save dot product
	double n_dot_v0 = 0;

	//deal with yz first

	//Check for reflection from far wedge plane - max 1 bounce
	if (dz > 0) {
		n_dot_v = -(dy*wedgeFarPlaneNy + dz*wedgeFarPlaneNz);
		n_dot_v0 = -(y*wedgeFarPlaneNy + z*wedgeFarPlaneNz);

		dt = -(wedgeFarPlaneD+n_dot_v0)/n_dot_v;
		//pretending y=0 at bar top for speed
		if (dt*dy < wedgeHeight) {
			//reflect it off bottom wedge
			//Does not pass through optical interface

			//Will always be true - just use for propagation (short circuit var?)

			mm_index += dt*quartzIndex;
			x_wedge_coerce_check(x,y,z,dx,dy,dz,dt);


			dy += 2*n_dot_v*wedgeFarPlaneNy;
			dz += 2*n_dot_v*wedgeFarPlaneNz;
		} else {
			n_dot_v = -(dy*upperWedgeFarPlaneNy + dz*upperWedgeFarPlaneNz);
			n_dot_v0 = -(y*upperWedgeFarPlaneNy + z*upperWedgeFarPlaneNz);

			dt = -(upperWedgeFarPlaneD+n_dot_v0)/n_dot_v;
			//never reached??? probably not  can possibly remove if statement
			//No bottom wedge, try top
			if (dt*dy < upperWedgeTop) { //Should always be true... I hope (remove later?)
				//Does pass through optical interface

				//Following statement performs the propagation if it does not fail
				mm_index += dt*quartzIndex;
				if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
					//Goes out the window, fail and return
					z = 1337;
					return -1;
				}

				//Reflect off top wedge
				dy += 2*n_dot_v*upperWedgeFarPlaneNy;
				dz += 2*n_dot_v*upperWedgeFarPlaneNz;
			}
		}
	}
	//Now dz < 0 or we have a new starting vector.  Either way, we intersect with the "close" wedge now
	n_dot_v = -(dy*wedgeClosePlaneNy + dz*wedgeClosePlaneNz);
	n_dot_v0 = -(y*wedgeClosePlaneNy + z*wedgeClosePlaneNz);

	dt = -(wedgeClosePlaneD+n_dot_v0)/n_dot_v;

	//Assume close enough to determine which wedge it hit
	double ty = dt*dy + y - barLength/2;

	if ((ty < wedgeHeight)) {
		//reflect it off bottom wedge

		//Again, always true
		mm_index += dt*quartzIndex;
		if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
			//Goes out the window, fail and return
			//Shouldn't
			z = 1337;
			return -1;
		}
		dy += 2*n_dot_v*wedgeClosePlaneNy;
		dz += 2*n_dot_v*wedgeClosePlaneNz;

		wedgeBeforeInterface = 1;

	} else if (ty > upperWedgeBottom && ty < upperWedgeTop) {
		//passed interface before reflection
		//get dt to propagate to interface
		dt = (quartzLiquidY + barLength/2 - y)/dy;
		mm_index += dt*quartzIndex;
		if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
			//Goes out the window, fail and return
			z = 1337;
			return -1;
		}
		//correct ordering below - interface is xz plane, so dx and dz go first
		optical_interface_z(quartzIndex,liquidIndex,dx,dz,dy);
		passed_interface = true;

		wedgeBeforeInterface = 0;

		//NOW intersect with upper wedge

		//Signs are weird here, it works, but I think they are twisted up with the initialization too
		n_dot_v = -(dy*upperWedgeClosePlaneNy + dz*upperWedgeClosePlaneNz + dx*upperWedgeClosePlaneNx);
		n_dot_v0 = -(y*upperWedgeClosePlaneNy + z*upperWedgeClosePlaneNz + x*upperWedgeClosePlaneNx);


		dt = -(upperWedgeClosePlaneD+n_dot_v0)/n_dot_v;

		if (dt*dy + y < upperWedgeTop + barLength/2) {
			mm_index += dt*liquidIndex;
			if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
				//Goes out the window, fail and return
				z = 1337;
				return -1;
			}

			if (y > wedgeHeight + barLength/2 && y < upperWedgeBottom + upperWedgeGap + barLength/2)
			{
				//We let it go out the gap on the slanted side, but not on the left and right.
				z = 1337;
				return -1;
			}
			
			if (upperWedgeAngleStore == true)
			{
				upper_wedge_incident.push_back(57.3*acos(n_dot_v));
			}
			
			dx += 2*n_dot_v*upperWedgeClosePlaneNx;
			dy += 2*n_dot_v*upperWedgeClosePlaneNy;
			dz += 2*n_dot_v*upperWedgeClosePlaneNz;

		} else {
			//refracted such that it no longer bounces
			dt = (upperWedgeTop + barLength/2 - y)/dy;
			mm_index += dt*liquidIndex;
			if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
				//Goes out the window, fail and return
				z = 1337;
				return -1;
			}
		}

	} else if (ty < upperWedgeTop) {
		//out the window
		z = 1337;
		return -1;
	}

	//Have now performed the maximum number of reflections (besides x direction)
	//Finish by taking it to the top
	if (!passed_interface == true) {
		//Will go through the interface now before it hits the top
		dt = (quartzLiquidY + barLength/2 - y)/dy;
		mm_index += dt*quartzIndex;
		if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
			//Goes out the window, fail and return
			z = 1337;
			return -1;
		}
		//correct ordering below - interface is xz plane, so dx and dz go first
		optical_interface_z(quartzIndex,liquidIndex,dx,dz,dy);
		passed_interface = true;
	}

	//Now we just finish
	//        printf("liquidIndex: %12.04f  quartzIndex: %12.04f\n",liquidIndex,quartzIndex);
	//	printf("%12.04f %12.04f %12.04f %12.04f %12.04f %12.04f\n",x,y,z,dx,dy,dz);

	if (wedge_bounces > 2)
	{	
		//		printf("Wedge Before Interface   :  bounces: %d  :  %d\n",wedgeBeforeInterface,wedge_bounces);
	}
	dt = (upperWedgeTop + barLength/2 - y)/dy;
	mm_index += dt*liquidIndex;
	if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt))) {
		//I'm not sure if it can go out the window at this point, but just in case
		//Goes out the window, fail and return
		z = 1337;
		return -1;
	}


	//and we're done.  Everything should be correct now with the ray at the top of the wedge
	//return the distance traveled times the index it traveled in

	return mm_index;

}
//Possibly inline these or something for speed, but right now, leave them for sanity
bool DircBaseSim::optical_interface_z(\
		double n1,\
		double n2,\
		double &dx,\
		double &dy,\
		double &dz) {
	//n1 is starting index of refraction, n2 is ending
	//Assume that it's going through the z plane
	//(this will change the ordering when called in our current geometry)

	//dz and dy are flipped, so this is really acos
	double before_ang = acos(dz);
	//	printf("%12.04f\n",before_ang);

	double n12rat = n1/n2;
	double n12rat2 = n12rat*n12rat;

	//takes care of trig
	dz = sqrt(1-n12rat2*(1-dz*dz));

	// 	if (dz != dz) return false;//total internal reflection

	//simpler expression than I expected
	dx *= n12rat;
	dy *= n12rat;

	// 	printf("in optical interface store_refraction: %s\n",store_refraction ? "true" : "false");
	if (store_refraction == true)
	{
		double after_ang = acos(dz);
		refraction_before.push_back(before_ang);
		refraction_after.push_back(after_ang);
	}

	return true;
}
bool DircBaseSim::x_wedge_coerce_check(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz,\
		double dt) {
	//already have dt, be sure to add it in when calling this

	//assumes x starts between the sides of the wedge
	//assumes dy/dx constant over the whole distance
	//also performs the index plane crossing

	//changes x and dx, correctly propagating them
	//propagates y and z as well
	//returns false if it goes out the window on the x step
	//Will not check for y and z directions (assumes dt corresponds to an intersection/reflection point)

	double cdt = 0;
	double tdt = 0;
	double ty = y - barLength/2;

	//This is bad and encapsulation breaking, but much faster and reduces code duplication
	if (midLineMode == true)
	{
		if (ty < wedgeHeight)
		{
			tdt = (wedgeHeight - ty)/dy;

			if (tdt > dt)
			{
				//does not go above wedge
				tdt = dt;
			}
			//only randomize if bouncing in infinitely thin wall
			//dx = midLineWedgeWallFlip*fabs(dx);
			//dx = midLineWedgeWallFlip*fabs(dx);
			//x starts at 0
			x = 0;
			midLineWedgeWallFlip *= -1;
		}
		else
		{
			//does not bounce
			tdt = dt;
			x += dx*tdt;
		}

		x += dx*(dt - tdt);
		y += dy*dt;
		z += dz*dt;
		//ignores possibility of bouncing out - only care about that below window anyway
		return true;
	}

	while (true) {
		//Get time of this bounce based on x direction
		if (dx < 0) {
			tdt = x + barWidth/2;
			tdt /= -dx;
			//coming from right wall
		} else {
			tdt = -x + barWidth/2 - wedgeWidthOff;
			tdt /= dx;
			//coming from left wall
		}
		if (tdt + cdt < dt) {
			//bounced
			//add to total time taken
			cdt += tdt;
			ty += dy*tdt;
			if (ty > wedgeHeight && ty < upperWedgeBottom) {
				//out the window on an edge
				//this is ok, finish propagation without bouncing
				//return false;
				x += tdt*dx;
				tdt = dt - cdt;
				break;
			}

			//not out the window, propagate x for the next iteration
			x += tdt*dx;
			if (ty < upperWedgeBottom) {
				//Must be less than Wedge height at this point
				//bounce off lower wedge
				dx = -dx;
			} else {
				//if in upper wedge, don't bounce and propagate like normal
				tdt = dt - cdt;
				break;
			}
		} else {
			//does not bounce - record remaining dt and break;
			tdt = dt - cdt;
			break;
		}
		wedge_bounces++;
	}
	//Finish last leg of trip:
	x += dx*tdt;
	y = barLength/2 + ty + dy*tdt;//Note that y starts at zero (bar middle), but ty starts at bar top
	z += dt*dz;//Implemented for convenience - z direction irrelevant to window calculation

	return true;//not out the window for the whole trip
}
void DircBaseSim::plane_reflect(\
		double Nx,\
		double Ny,\
		double Nz,\
		double D,\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz,\
		double &dt,\
		double offang /*=0*/) {
	//Propagate to the intercept and reflects off a plane
	//Could use this in the Wedge, but would loose some speed

	double n_dot_v = (Nx*dx + Ny*dy + Nz*dz);
	double n_dot_v0 = (Nx*x + Ny*y + Nz*z);

	dt = (D - n_dot_v0)/n_dot_v;

	x += dt*dx;
	y += dt*dy;
	z += dt*dz;

	if (fabs(offang) > 1e-8) {
		//remove for speed obviously
		rotate_2d(Nz,Ny,cos(offang),sin(offang));
	}

	dx -= 2*n_dot_v*Nx;
	dy -= 2*n_dot_v*Ny;
	dz -= 2*n_dot_v*Nz;
}
double DircBaseSim::get_z_intercept(\
		double Nx,\
		double Ny,\
		double Nz,\
		double D,\
		double x,\
		double y,\
		double z,\
		double dx,\
		double dy,\
		double dz) {
	//finds z intercept (without modifying x y and z)
	//Could use this in the Wedge, but would loose some speed

	double n_dot_v = (Nx*dx + Ny*dy + Nz*dz);
	double n_dot_v0 = (Nx*x + Ny*y + Nz*z);
	double dt;//optimized by math compression?
	dt = (D - n_dot_v0)/n_dot_v;

	return z + dt*dz;
}
double DircBaseSim::get_intercept_plane(\
		double Nx,\
		double Ny,\
		double Nz,\
		double D,\
		double &x,\
		double &y,\
		double &z,\
		double dx,\
		double dy,\
		double dz) {
	//Just intersects a plane (modifying x,y,z)
	//Could use this in the Wedge, but would loose some speed
	//Non returning verion for speed when needed?

	double n_dot_v = (Nx*dx + Ny*dy + Nz*dz);
	double n_dot_v0 = (Nx*x + Ny*y + Nz*z);
	double dt;//optimized by math compression?
	dt = (D - n_dot_v0)/n_dot_v;

	x += dx*dt;
	y += dy*dt;
	z += dz*dt;

	return dt;
}
double DircBaseSim::sgn(double val) {
	//stolen from http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
	return (0 < val) - (val < 0);
}

