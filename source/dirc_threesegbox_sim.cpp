#include <vector>

#include "../include/dirc_threesegbox_sim.h"
#include "../include/dirc_point.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
//
#include <stdlib.h>
#include <math.h>

#include <TRandom3.h>



DircThreeSegBoxSim::DircThreeSegBoxSim(
		int rand_seed /*=4357*/,\
		double ifoc_r/*=540.66*/, \
		double ifoc_mirror_size/*=300.38*/, \
		double ifoc_rot/*=-74.11*/, \
		double isens_size/*=600*/, \
		double isens_rot/*=90*/,\
		double ibar_length/*=4900*/,\
		double ibar_width/*=35*/,\
		double ibar_depth/*=17.25*/,
		double iupper_wedge_top/*=178.6*/) : 
			DircBaseSim(
				rand_seed,\
				ibar_length,\
				ibar_width,\
				ibar_depth) {

	printf("BarLWD: %12.04f %12.04f %12.04f\n",barLength,barWidth,barDepth);
	foc_r = ifoc_r;
	foc_mirror_size = ifoc_mirror_size;
	foc_rot = ifoc_rot;
	foc_yrot = 0;
	sens_size = isens_size;
	sens_rot = isens_rot;

	storeOpticalAngles = false;

	//only used for checking collision
	largePlanarMirrorNx = 0; //Shouldn't be needed
	largePlanarMirrorNy = 1;
	largePlanarMirrorNz = 0;
	largePlanarMirrorD = upperWedgeTop + barLength/2;
	largePlanarMirrorMinZ = -559;
	largePlanarMirrorMaxZ = -130;

	pmtPlaneMinZ = -559;
	pmtPlaneMaxZ = -329;

	focMirrorBottom = 139 + upperWedgeTop + barLength/2;


	sidemirror_xl = -1000000;
	sidemirror_xr = 1000000;
	sidemirror_reflectivity = .9;


	//Take average
	quartzIndex = 1.47;
	liquidIndex = 1.47;
	quartzLiquidY = upperWedgeBottom;


	boxCloseZ = -614;

	reflOff = 9;

	three_seg_mirror = false;

	num_QE = 31;
	min_QE = 300;
	max_QE = 600;
	sep_QE = (max_QE - min_QE)/(num_QE - 1);

	double t_QE[31] = {\
		0.016415, 0.074064, 0.141658, 0.184219, 0.20634,  0.217191, 0.223244,
	       0.222296, 0.215232, 0.206302, 0.195574, 0.183007, 0.169403, 0.155447,
	       0.140014, 0.127696, 0.115716, 0.104086, 0.092256, 0.084083, 0.075418,
	       0.067311, 0.060243, 0.053588, 0.047765, 0.04344,  0.037999, 0.034177,
	       0.030869, 0.027848, 0.024141
	};


	// Transmittance of quartz per 1m


	for (int i = 0; i < num_QE; i++) {
		vals_QE.push_back(t_QE[i]);
	}

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
	set_focmirror_nonuniformity(0);
	build_readout_box();
}
double DircThreeSegBoxSim::get_cerenkov_angle_rand(double beta, double additional_spread, double &wavelength) {
        //May be slow enough to consider approximating in distribution generation
        double out_ang = 0;
        double tmp_lam = 0;
        double tmp_QE_val;
        double above_ind;
        int ind_QE;
        double n_lam;

        while (true) {
                tmp_lam = rand_gen->Uniform(min_QE,max_QE);
                wavelength = tmp_lam;

                //Ignoring that the QE points are right on multiples of 10.  Assuming they are for speed.
                //This may not be neccessary, but I doubt it matters.
                ind_QE = (tmp_lam - min_QE)/sep_QE;

                above_ind = tmp_lam - (min_QE + sep_QE*ind_QE);

                //Simple linear interp between values.  5th order poly fit looked like it worked too
                tmp_QE_val = vals_QE[ind_QE]*(sep_QE-above_ind)/sep_QE + vals_QE[ind_QE+1]*above_ind/sep_QE;

                //Max QE val is ~.23, this saves lot of loops
                if (rand_gen->Uniform(0,.25) > tmp_QE_val) continue;

                //Test emission distribution, second b/c it's a less stringent cut
                if (rand_gen->Uniform(0,1/(min_QE*min_QE)) > 1/(tmp_lam*tmp_lam)) continue;


                n_lam = get_quartz_n(tmp_lam);


                out_ang = 57.3*acos(1/(beta*n_lam));
                break;
        }

        out_ang += rand_gen->Gaus(0,additional_spread);

        return out_ang;
}

void DircThreeSegBoxSim::build_readout_box()
{
	build_system();
	fill_threeseg_plane_vecs();
	fill_foc_mirror_vecs();
	fill_sens_plane_vecs();
	//still rebuild wedge and bars when this is called:
}
void DircThreeSegBoxSim::set_pmt_plane_zs(double imin, double imax)
{
	pmtPlaneMinZ = imin;
	pmtPlaneMaxZ = imax;
}
void DircThreeSegBoxSim::set_large_mirror_zs(double imin, double imax)
{
	largePlanarMirrorMinZ = imin;
	largePlanarMirrorMaxZ = imax;
}
void DircThreeSegBoxSim::fill_sens_plane_vecs() {
	double adjusted_sens_size = 312;

	sensPlaneNx = 0;
	sensPlaneNy = sin(sens_rot/57.3);
	sensPlaneNz = cos(sens_rot/57.3);

	sensPlaneY = -adjusted_sens_size*sin(sens_rot/57.3)/2-reflOff+barLength/2;
	sensPlaneZ = boxCloseZ + sens_size*cos(sens_rot/57.3)/2;

	sensPlaneD = sensPlaneNy*sensPlaneY + sensPlaneNz*sensPlaneZ;

	sensPlaneYdistConversion = 1/cos(sens_rot/57.3);
	sensPlaneZdistConversion = 1/sin(sens_rot/57.3);
}
void DircThreeSegBoxSim::set_sidemirror_reflectivity(double isr) {
	sidemirror_reflectivity = isr;
}
void DircThreeSegBoxSim::set_foc_mirror_r(double ifoc_r) {
	foc_r = ifoc_r;
	build_readout_box();
}
void DircThreeSegBoxSim::fill_foc_mirror_vecs() {
	//is off by pi/2 to reduce rounding errors and such
	double foc_center_ang = foc_rot/57.3 + acos(-foc_mirror_size/(2*foc_r));
	focMirrorY = focMirrorBottom - foc_r*cos(foc_center_ang);
	focMirrorZ = foc_r*sin(foc_center_ang);
}
void DircThreeSegBoxSim::fill_threeseg_plane_vecs() {
	focMirrorTop = focMirrorBottom + foc_mirror_size*cos(foc_rot/57.3);
	focMirrorZDim = foc_mirror_size*sin(foc_rot/57.3);
	//If we ever go to more than 3 segments, use a loop

	double theta_m = foc_rot/57.3;//radians of rotation to go through
	double theta_c = fabs(2*asin(focMirrorZDim/(2*foc_r)));//angle subtended by the mirror
	double seg_h = fabs(2*foc_r*sin(theta_c/6));//length of each segment;

	//I had to do some geometry and algebra to get these numbers, but in hindsight, it's obvious.  Always that way.
	double theta_1 = theta_m - theta_c/3;
	double theta_2 = theta_m;
	double theta_3 = theta_m + theta_c/3;

	threeSeg1Y = focMirrorBottom;
	threeSeg1Z = 0;

	threeSeg2Y = threeSeg1Y + seg_h*cos(theta_1);
	threeSeg2Z = threeSeg1Z - seg_h*sin(theta_1);

	threeSeg3Y = threeSeg2Y + seg_h*cos(theta_2);
	threeSeg3Z = threeSeg2Z - seg_h*sin(theta_2);


	threeSeg1Nx = 0;
	threeSeg1Ny = sin(theta_1);
	threeSeg1Nz = cos(theta_1);
	rotate_2d(threeSeg1Nx,threeSeg1Nz,cos(foc_yrot/57.3),sin(foc_yrot/57.3));//I think this is a slightly wrong rotation if the mirrors are carved out of a solid block, but it should be good enough at small angles
	rotate_2d(threeSeg1Nx,threeSeg1Ny,cos(foc_zrot/57.3),sin(foc_zrot/57.3));//I think this is a slightly wrong rotation if the mirrors are carved out of a solid block, but it should be good enough at small angles
	threeSeg1D = threeSeg1Ny*threeSeg1Y + threeSeg1Nz*threeSeg1Z;//Use point x=0 as reference

	threeSeg2Nx = 0;
	threeSeg2Ny = sin(theta_2);
	threeSeg2Nz = cos(theta_2);
	rotate_2d(threeSeg2Nx,threeSeg2Nz,cos(foc_yrot/57.3),sin(foc_yrot/57.3));
	rotate_2d(threeSeg2Nx,threeSeg2Ny,cos(foc_zrot/57.3),sin(foc_zrot/57.3));//I think this is a slightly wrong rotation if the mirrors are carved out of a solid block, but it should be good enough at small angles
	threeSeg2D = threeSeg2Ny*threeSeg2Y + threeSeg2Nz*threeSeg2Z;//Use point x=0 as reference

	threeSeg3Nx = 0;
	threeSeg3Ny = sin(theta_3);
	threeSeg3Nz = cos(theta_3);
	rotate_2d(threeSeg3Nx,threeSeg3Nz,cos(foc_yrot/57.3),sin(foc_yrot/57.3));
	rotate_2d(threeSeg3Nx,threeSeg3Ny,cos(foc_zrot/57.3),sin(foc_zrot/57.3));//I think this is a slightly wrong rotation if the mirrors are carved out of a solid block, but it should be good enough at small angles
	threeSeg3D = threeSeg3Ny*threeSeg3Y + threeSeg3Nz*threeSeg3Z;//Use point x=0 as reference

}
void DircThreeSegBoxSim::set_focmirror_nonuniformity(double nonuni_deg) {
	foc_mirror_nonuni = nonuni_deg;
	nonUniformFocMirror = (fabs(nonuni_deg) > .001);
}
void DircThreeSegBoxSim::sidemirror_reflect_points(std::vector<dirc_point> &points) {
	double tmpx = 0;
	for (unsigned int i = 0; i < points.size(); i++) {
		tmpx = points[i].x;
		while (tmpx < sidemirror_xl || tmpx > sidemirror_xr) {
			if (tmpx < sidemirror_xl) {
				tmpx = 2*sidemirror_xl - tmpx;
			}
			if (tmpx > sidemirror_xr) {
				tmpx = 2*sidemirror_xr - tmpx;
			}
		}
		points[i].x = tmpx;
	}
}
void DircThreeSegBoxSim::sidemirror_reflect_point(dirc_point &point) {
	double tmpx = 0;
	tmpx = point.x;
	while (tmpx < sidemirror_xl || tmpx > sidemirror_xr) {
		//printf("%12.04f ",tmpx);
		if (rand_gen->Uniform(0,1) > sidemirror_reflectivity)
		{
			point.t = -1337;
			break;
		}
		//printf("%12.04f ",tmpx);
		if (tmpx < sidemirror_xl) {
			tmpx = 2*sidemirror_xl - tmpx;
		}
		if (tmpx > sidemirror_xr) {
			tmpx = 2*sidemirror_xr - tmpx;
		}
	}
	point.x = tmpx;
	//printf("%12.04f \n",tmpx);
}
void DircThreeSegBoxSim::set_sidemirror(double ixr, double ixl) {
	sidemirror_xr = ixr;
	sidemirror_xl = ixl;
}
void DircThreeSegBoxSim::set_three_seg_mirror(bool itsm) {
	if (three_seg_mirror != itsm) {
		three_seg_mirror = itsm;
		build_readout_box();
	}
}
void DircThreeSegBoxSim::set_pmt_offset(double r) {
	//positive r increases path length
	boxCloseZ = -614 - r*cos(sens_rot/57.3);
	reflOff = 9 + r*sin(sens_rot/57.3);
	build_readout_box();
}
void DircThreeSegBoxSim::set_liquid_absorbtion(double iabs) {
	liquidAbsorbtion = iabs;
}
std::vector<double> DircThreeSegBoxSim::get_dist_traveled() {
	return dist_traveled;
}
void DircThreeSegBoxSim::set_store_traveled(bool sst/* = true*/) {
	store_traveled = sst;
}
//Liquid Index is the same as the upper quartz wedge - the volumes are connected
void DircThreeSegBoxSim::set_focus_mirror_angle(double ang,double yang, double zang) {
	foc_rot = ang;
	foc_yrot = yang;
	foc_zrot = zang;
	build_readout_box();
}
void DircThreeSegBoxSim::set_pmt_angle(double ang) {
	sens_rot = ang;
	build_readout_box();
}
void DircThreeSegBoxSim::set_store_optical_angles(bool istore)
{
	storeOpticalAngles = istore;
}
std::vector<double> DircThreeSegBoxSim::get_focus_photon_angles()
{
	return focus_photon_angles;
}
std::vector<double> DircThreeSegBoxSim::get_large_flat_photon_angles()
{
	return large_flat_photon_angles;
}
std::vector<double> DircThreeSegBoxSim::get_side_photon_angles()
{
	return side_photon_angles;
}
void DircThreeSegBoxSim::warp_readout_box(
	dirc_point &out_val,\
	int particle_bar,\
	double &mm_index,\
	double &x,\
	double &y,\
	double &z,\
	double &dx,\
	double &dy,\
	double &dz)
{
	//printf("inside warp_box: %12.04f %12.04f\n",quartzIndex,liquidIndex);
	double c_mm_ns = 300;
	mm_index += warp_box(\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);

	if (z > 0) {
		out_val.t = -1337;
		return;
	}

	//check absorbtion
	if (!(absorbtion_mc(dx,dy))) {
		out_val.t = -1337;
		return;
	}

	mm_index += warp_sens_plane(\
			out_val,\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);

	if (z > 0) {
		out_val.t = -1337;
		return;
	}

	out_val.t = mm_index/(c_mm_ns);
	out_val.x += get_bar_offset(particle_bar);
	
	//Must reflect after offset....
	sidemirror_reflect_point(out_val);
}
bool DircThreeSegBoxSim::absorbtion_mc(double dx, double dy) {
	//True if the particle makes it
	//expects vector pointing after bounce
	//Magic number corresponding to minimal distancel traveled
	double min_dist = 650+56;

	//approximating here - assume yz distance is independent of where on mirror it hit
	//measuring side plots with a ruler seems to mean this is good to ~5-10%
	double approx_dist = min_dist*sqrt(dx*dx+dy*dy)/dy;

	if (store_traveled == true) {
		dist_traveled.push_back(approx_dist);
	}

	double prob_transmitted = exp(-liquidAbsorbtion*approx_dist);
	if (rand_gen->Uniform(0,1) < prob_transmitted) {
		return true;
	} else {
		return false;
	}
}
double DircThreeSegBoxSim::warp_box(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz) {
	//propagates light from the top of the wedge until just after bouncing off of the focusing mirrors
	//Consider branch profiling or some such other optimizations
	//Also, Fast Math, but don't tell the people in the penthouse

	//does not currently include x bouncing off the side - need to add that later
	//possibly in the last "warp to plane" bit
	double rval = 0; //distance traveled

	//first reflect off the back of the bar
	double dt;

	if (dz > 0) {
		//Condition may not be needed - try without for speed later
		dt = -z/dz;

		if (y + dy*dt < focMirrorBottom) {
			//reflects off of back
			x += dx*dt;
			y += dy*dt;
			z += dz*dt;

			rval += dt*liquidIndex;

			//Abuse the fact that this vector is normalized and hope it stays that way
			dz = -dz;
			//rval += dt;
		}
	}


	double trval  = 0;

	if (three_seg_mirror == true) {
		trval = three_seg_reflect(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);
	} else {
		trval = cylindrical_reflect(\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz);
	}

	if (trval < 0) {
		z = 1337;
		return -1;
	}

	//short propagate in front of the mirrors:
	//TODO remove later for speed
	x += .1*dx;
	y += .1*dy;
	z += .1*dz;


	rval += trval;

	return rval;
}
double DircThreeSegBoxSim::three_seg_reflect(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz) {
	//Intersects and reflects the three segment plane
	double rval = 0;
	//definitely losing time here - combuting the n_dot_v and n_dot_v0 and dt twice for the chosen plane
	//TODO fix this once the code is correct

	//I hope there's a fast way to do these reflections

	double tz = 0;
	//printf("inside three seg: %12.04f %12.04f\n",quartzIndex,liquidIndex);

	//check first seg (again, loop this if we go more than 3)
	tz = get_z_intercept(\
			threeSeg1Nx,\
			threeSeg1Ny,\
			threeSeg1Nz,\
			threeSeg1D,\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);
	double offang = 0;
	if (tz > threeSeg2Z && tz < 0) {
		//reflect off mirror closest to box
		if (nonUniformFocMirror == true) {
			//obviously not in the real run
			offang = rand_gen->Gaus(0,foc_mirror_nonuni/57.3);
			// 			rotate_2d(threeSeg1Nz,threeSeg1Ny,cos(offang),sin(offang));
		}
		if (storeOpticalAngles == true)
		{
			double ndotv = dx*threeSeg1Nx + dy*threeSeg1Ny + dz*threeSeg1Nz;
			focus_photon_angles.push_back(57.3*acos(ndotv));
		}
		plane_reflect(\
				threeSeg1Nx,\
				threeSeg1Ny,\
				threeSeg1Nz,\
				threeSeg1D,\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz,\
				rval,\
				offang);
		return rval*liquidIndex;
	}

	tz = get_z_intercept(\
			threeSeg2Nx,\
			threeSeg2Ny,\
			threeSeg2Nz,\
			threeSeg2D,\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);
	// 	printf("tz2: %12.04f\n",tz);
	if (tz > threeSeg3Z && tz < 0) {
		//reflect off middle mirror
		if (nonUniformFocMirror == true) {
			//obviously not in the real run
			offang = rand_gen->Gaus(0,foc_mirror_nonuni/57.3);
			// 			rotate_2d(threeSeg2Nz,threeSeg2Ny,cos(offang),sin(offang));
		}
		if (storeOpticalAngles == true)
		{
			double ndotv = dx*threeSeg2Nx + dy*threeSeg2Ny + dz*threeSeg2Nz;
			focus_photon_angles.push_back(57.3*acos(ndotv));
		}
		
		plane_reflect(\
				threeSeg2Nx,\
				threeSeg2Ny,\
				threeSeg2Nz,\
				threeSeg2D,\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz,\
				rval,\
				offang);
		return rval*liquidIndex;
	}

	tz = get_z_intercept(\
			threeSeg3Nx,\
			threeSeg3Ny,\
			threeSeg3Nz,\
			threeSeg3D,\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);
	// 	printf("tz3: %12.04f\n",tz);
	if (tz > -focMirrorZDim && tz < 0) {
		//reflect off mirror closest to box
		if (nonUniformFocMirror == true) {
			//obviously not in the real run
			offang = rand_gen->Gaus(0,foc_mirror_nonuni/57.3);
			// 			printf("%12.04e\n",offang);
			// 			rotate_2d(threeSeg3Nz,threeSeg3Ny,cos(offang),sin(offang));
		}
		if (storeOpticalAngles == true)
		{
			double ndotv = dx*threeSeg3Nx + dy*threeSeg3Ny + dz*threeSeg3Nz;
			focus_photon_angles.push_back(57.3*acos(ndotv));
		}
		plane_reflect(\
				threeSeg3Nx,\
				threeSeg3Ny,\
				threeSeg3Nz,\
				threeSeg3D,\
				x,\
				y,\
				z,\
				dx,\
				dy,\
				dz,\
				rval,\
				offang);
		return rval*liquidIndex;
	}
	//reflects off of nothing :(
	z = 1337;
	return -1;
}
double DircThreeSegBoxSim::cylindrical_reflect(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz) {
	//intersects and reflects the focusing cylinder
	//Pretending the cylinder has no Nx for intersection purposes
	//Should be a valid approximation, and saves a ton of speed
	//could be implemented on the threeseg mirror (at least the branching part

	double rval = 0;

	double dydz_norm = sqrt(dz*dz+dy*dy);
	double dy_norm = dy/dydz_norm;
	double dz_norm = dz/dydz_norm;

	double localNx = 0;
	double localNy = 0;
	double localNz = 0;

	//there's gotta be a faster way to do all of this
	//different than plane intercept D.  Took this algorithm from internet (wolfram Mathworld)
	//Parametric intercept sucks
	double D = (z - focMirrorZ)*dy_norm - (y - focMirrorY)*dz_norm;
	double detD = sqrt(foc_r*foc_r - D*D);

	//dy > 0 always (or should be.  Remove sgn and fabs function calls after confirming)
	double zrel = D*dy_norm + sgn(dy_norm)*dz_norm*detD;
	double yrel = -D*dz_norm + fabs(dy_norm)*detD;

	double newy = yrel + focMirrorY;
	double newz = zrel + focMirrorZ;

	double ydiff = newy - y;
	double zdiff = newz - z;

	rval += sqrt((ydiff*ydiff+zdiff*zdiff))/dydz_norm;

	x += dx*rval;
	y = newy;
	z = newz;

	//combine later to save speed
	localNy = yrel;
	localNz = zrel;

	if (nonUniformFocMirror == true) {
		//obviously not in the real run
		double offang = rand_gen->Gaus(0,foc_mirror_nonuni/57.3);
		rotate_2d(localNz,localNy,cos(offang),sin(offang));
	}


	double norm_loc = sqrt(localNy*localNy + localNz*localNz);//there's gotta be a better way to normalize this

	localNy /= norm_loc;
	localNz /= norm_loc;

	//remove for speed in ideal case
	// 	rotate_2d(localNx,localNz,cos(foc_yrot/57.3),sin(foc_yrot/57.3));

	double n_dot_v = -(dx*localNx + dy*localNy + dz*localNz);

	// 	printf("dx: %8.04f dy: %8.04f dz: %8.04f\n",dx,dy,dz);

	if (storeOpticalAngles == true)
	{
		focus_photon_angles.push_back(57.3*acos(n_dot_v));
	}
	dx += 2*n_dot_v*localNx;
	dy += 2*n_dot_v*localNy;
	dz += 2*n_dot_v*localNz;

	// 	printf("dx: %8.04f dy: %8.04f dz: %8.04f\n",dx,dy,dz);

	//printf("inside cyl: %12.04f %12.04f\n",quartzIndex,liquidIndex);

	return rval*liquidIndex;
}

double DircThreeSegBoxSim::warp_sens_plane(\
		dirc_point &fill_val,\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz) {
	//don't strictly need to modify the z, could be sped up
	//First check to see if it goes through the large planar mirror - it probably doesn't
	
	double tmpz = get_z_intercept(\
			largePlanarMirrorNx,\
			largePlanarMirrorNy,\
			largePlanarMirrorNz,\
			largePlanarMirrorD,\
			x,\
			y,\
			z,\
			dx,\
			dy,\
			dz);

	if (tmpz < largePlanarMirrorMinZ || tmpz > largePlanarMirrorMaxZ)
	{
		z=1337;
		return 100000;
	}
	double rval = get_intercept_plane(\
				     sensPlaneNx,\
				     sensPlaneNy,\
				     sensPlaneNz,\
				     sensPlaneD,\
				     x,\
				     y,\
				     z,\
				     dx,\
				     dy,\
				     dz);
	if (z < pmtPlaneMinZ || z > pmtPlaneMaxZ)
	{
		z=1337;
		return 100000;
	}
	if (storeOpticalAngles == true)
	{
		large_flat_photon_angles.push_back(57.3*acos(fabs(dy)));
	}
	if (storeOpticalAngles == true)
	{
		//Assumes all Photons bounce - not exatly true, but close
		side_photon_angles.push_back(57.3*acos(dx));
	}

	fill_val.x = x;
	//fill_val.y = (y-sensPlaneY)*sensPlaneYdistConversion;
	fill_val.y = (-z-559)*sensPlaneYdistConversion + 220;

	return rval*liquidIndex;
}

