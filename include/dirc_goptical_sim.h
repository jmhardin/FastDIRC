#include <vector>
#include <utility>
#include <TRandom3.h>

#include "dirc_point.h"
#include "dirc_optical_components.h"

#include <Goptical/Sys/System>



#ifndef DIRC_GOPTICAL_SIM
#define DIRC_GOPTICAL_SIM 
using namespace Goptical;
class DircGopticalSim
{
private:
	double foc_r;
	double foc_mirror_size;
	double foc_rot;
	double foc_yrot;
	double sens_size;
	double sens_rot;
	
	double barLength;
	double barWidth;
	double barDepth;
	double wedgeWidthOff;
	double wedgeDepthOff;
	double wedgeFarAngle;
	double wedgeCloseAngle;
	double wedgeWidth;
	double wedgeDepthHigh;
	double wedgeHeight;
	double upperWedgeDepthHigh;
	double upperWedgeTop;
	double upperWedgeHeight;
	double upperWedgeBottom;
	
	double wedgeClosePlaneNx;
	double wedgeClosePlaneNy;
	double wedgeClosePlaneNz;
	double wedgeClosePlaneD;

	double upperWedgeClosePlaneNx;
	double upperWedgeClosePlaneNy;
	double upperWedgeClosePlaneNz;
	double upperWedgeClosePlaneD;
	double lowerWedgeExtensionZ;
	
	
	double focMirrorBottom;
	double focMirrorTop;
	double focMirrorZDim;
	//Multiseg?  probably not.  If it goes up again, use arrays
	double threeSeg1Nx,threeSeg1Ny,threeSeg1Nz,threeSeg1D;
	double threeSeg2Nx,threeSeg2Ny,threeSeg2Nz,threeSeg2D;
	double threeSeg3Nx,threeSeg3Ny,threeSeg3Nz,threeSeg3D;
	
	//Not yet used - implement to branch faster on the threeseg
	double threeSeg1_2dny,threeSeg1_2dnz,threeSeg1_2dd;
	double threeSeg2_2dny,threeSeg2_2dnz,threeSeg2_2dd;
	double threeSeg3_2dny,threeSeg3_2dnz,threeSeg3_2dd;
	
	double threeSeg1Y,threeSeg1Z;
	double threeSeg2Y,threeSeg2Z;
	double threeSeg3Y,threeSeg3Z;
	
	double sensPlaneNx;
	double sensPlaneNy;
	double sensPlaneNz;
	double sensPlaneD;
	double sensPlaneY;
	double sensPlaneZ;
	
	double sensPlaneYdistConversion;
	
	bool upperWedgeNonUniform;
	double upperWedgeNonUniformSpread;
	
	double wedgeFarPlaneNx;
	double wedgeFarPlaneNy;
	double wedgeFarPlaneNz;
	double wedgeFarPlaneD;
	
	double upperWedgeFarPlaneNx;
	double upperWedgeFarPlaneNy;
	double upperWedgeFarPlaneNz;
	double upperWedgeFarPlaneD;
	
	double upperWedgeFarZ;
	
	double boxCloseZ;
	double reflOff;
	
	double focMirrorY;
	double focMirrorZ;
	
	bool three_seg_mirror;
	
	double sidemirror_xr;
	double sidemirror_xl;
	
	double quartzIndex;
	double liquidIndex;
	double quartzLiquidY;
	
	double box_angle_off_cval;
	double box_angle_off_sval;
	
	double liquidAbsorbtion;
	std::vector<double> dist_traveled;
	bool store_traveled;
	
	bool store_refraction;
	std::vector<double> refraction_before;
	std::vector<double> refraction_after;
	
	double min_QE,max_QE,sep_QE;
	int num_QE;
	std::vector<double> vals_QE;
	double min_transmittance,max_transmittance,sep_transmittance;
	int num_transmittance;
	std::vector<double> quartz_transmittance;
	
	Sys::System sys;
// 	Sys::SourceRays srcrays;
	
	ref<qwartzCrystal> qwartz;
	ref<mineralOil> oil;
	ref<Sys::Image> boxImage;
// 	ref<_Goptical::Material::vacuum> goptical_vacuum;
	
	TRandom3 *rand_gen;
	
	void rotate_2d(double &x, double &y, double cval, double sval);
	
	
	void build_system();
	void clear_system();
	void fill_sens_plane_vecs();
	void fill_threeseg_plane_vecs();
	void fill_foc_mirror_vecs();
	void sidemirror_reflect_points(std::vector<dirc_point> &points);
	void spread_wedge_mirror();
	bool quartz_transmission_mc(double R, double lambda);
	bool absorbtion_mc(double dx, double dy);
	std::vector<dirc_point> trace_source_rays(\
		Sys::SourceRays &srcrays, \
		bool outspot, \
		bool outframe);
	void fill_rand_phi(\
		Sys::SourceRays &srcrays,\
		int n_photons, \
		double ckov_theta /*= 47*/, \
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc /*= .0015*57.3*/,\
		double ckov_theta_unc /* = .0055*57.3*/,\
		double beta /* = -1*/,\
		double check_dir /* = 0 */);
	void fill_reg_phi(\
		std::vector<dirc_point> &fill_points,\
		int n_photons_phi, \
		int n_photons_z,\
		double ckov_theta /*= 47*/, \
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc, /*= 0*/
		double ckov_theta_unc /* = 0*/,\
		double beta /* = -1*/,\
		double check_dir /* = 0 */);
	
	double get_quartz_n(double lambda);
	bool optical_interface_z(\
		double n1,\
		double n2,\
		double &dx,\
		double &dy,\
		double &dz);
	double warp_ray(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz,\
		double critical_angle);
	double warp_wedge(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz);
	double warp_box(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz);
	bool x_wedge_coerce_check(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz,\
		double dt);
	double three_seg_reflect(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz);
	//Utility function - should definitely be inlined
	void plane_reflect(\
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
		double &dt);
	//Utility function - should combine this with above for some speed
	double get_z_intercept(\
		double Nx,\
		double Ny,\
		double Nz,\
		double D,\
		double x,\
		double y,\
		double z,\
		double dx,\
		double dy,\
		double dz);
	//yet another inlinable utility function (make a second with not return?)
	double get_intercept_plane(\
		double Nx,\
		double Ny,\
		double Nz,\
		double D,\
		double &x,\
		double &y,\
		double &z,\
		double dx,\
		double dy,\
		double dz);
	//yet another inlinable utility function
	double cylindrical_reflect(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz);
	//The compiler should be inlining this without our help
	double sgn(double val);
	//Technically one of the "warped functions, but can and should be optimized somehow
	double warp_sens_plane(\
		dirc_point &fill_val,\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz);
public:
	void set_sidemirror(double ixr, double ixl);
	void set_three_seg_mirror(bool itsm);
	void set_pmt_offset(double r);
	void set_liquid_absorbtion(double iabs);
	std::vector<double> get_dist_traveled();
	void set_store_traveled(bool sst = true);
	void set_liquid_index(double li);
	void set_bar_box_angle(double ang);
	void set_focus_mirror_angle(double ang,double yang = 0);
	void set_pmt_angle(double ang);
	void set_wedge_mirror_rand(double ispread);
	double get_cerenkov_angle_rand(double beta, double additional_spread, double &wavelength);
	double get_beta(double E, double m);
	void set_upper_wedge_angle_diff(double rads, double radsy_y = 0);	
	DircGopticalSim(\
		int rand_seed = 4357,\
		double ifoc_r = -1200, \
		double ifoc_mirror_size = 300.38, \
		double ifoc_rot = 74.11, \
		double isens_size = 600, \
		double isens_rot = 47.87);
	std::vector<std::pair<double,double> > get_refraction_rand_phi(\
		std::vector<double> &before_interface,\
		std::vector<double> &after_interface,\
		int n_photons, \
		double ckov_theta = 47, \
		double particle_x = 0, \
		double particle_y = 0, \
		double particle_theta = 0, \
		double particle_phi = 0,\
		double phi_theta_unc = .0015*57.3,\
		double ckov_theta_unc = .0055*57.3,\
		double check_dir = 0,\
		double beta = -1);
	std::vector<dirc_point> sim_rand_n_photons(\
		int n_photons,\
		bool outspot = false, \
		bool outframe = false, \
		double ckov_theta = 47, \
		double particle_x = 0, \
		double particle_y = 0, \
		double particle_theta = 0, \
		double particle_phi = 0,\
		double phi_theta_unc = .08594,\
		double ckov_theta_unc = .3151,\
		double beta = -1,\
		bool check_dir = false);	
	std::vector<dirc_point> sim_reg_n_photons(\
		int n_photons_phi,\
		int n_photons_z,\
		bool outfile = false, \
		bool outframe = false, \
		double ckov_theta = 47, \
		double particle_x = 0, \
		double particle_y = 0, \
		double particle_theta = 0, \
		double particle_phi = 0,\
		double phi_theta_unc = 0,\
		double ckov_theta_unc = 0,\
		double beta = -1,\
		bool check_dir = false);	
};
#endif
