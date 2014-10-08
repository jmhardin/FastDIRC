#include <vector>

#include "../include/dirc_goptical_sim.h"
#include "../include/dirc_optical_components.h"
#include "../include/dirc_point.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <Goptical/Analysis/Spot>

#include <Goptical/Math/Vector>
#include <Goptical/Math/VectorPair>

#include <Goptical/Sys/Element>
#include <Goptical/Sys/Container>
#include <Goptical/Sys/System>
#include <Goptical/Sys/SourcePoint>
#include <Goptical/Sys/SourceRays>
#include <Goptical/Sys/Source>
#include <Goptical/Sys/Group>
#include <Goptical/Design/Telescope/Newton>
#include <Goptical/Sys/Image>
#include <Goptical/Sys/Lens>
#include <Goptical/Sys/OpticalSurface>
#include <Goptical/Sys/Mirror>
#include <Goptical/Sys/Mirror>
#include <Goptical/Data/Plot>
#include <Goptical/Data/PlotData>
#include <Goptical/Data/Set>

#include <Goptical/Trace/Tracer>
#include <Goptical/Trace/Result>
#include <Goptical/Trace/Distribution>
#include <Goptical/Trace/Sequence>
#include <Goptical/Trace/Params>
#include <Goptical/Trace/Ray>

#include <Goptical/Light/Ray>
#include <Goptical/Light/SpectralLine>
#include "Goptical/common.hh"

#include <Goptical/Io/RendererSvg>
#include <Goptical/Io/Rgb>
#include <Goptical/Io/RendererViewport>

#include <Goptical/Math/Transform>

#include <stdlib.h>
#include <math.h>

#include <TRandom3.h>


using namespace Goptical;

//Using close simplex 
DircGopticalSim::DircGopticalSim(
	int rand_seed /*=4357*/,\
	double ifoc_r/*=540.66*/, \
	double ifoc_mirror_size/*=300.38*/, \
	double ifoc_rot/*=-74.11*/, \
	double isens_size/*=600*/, \
	double isens_rot/*=90*/)
{
	foc_r = ifoc_r;
	foc_mirror_size = ifoc_mirror_size;
	foc_rot = ifoc_rot;
	foc_yrot = 0;
	sens_size = isens_size;
	sens_rot = isens_rot;

	
	
	barLength=4900;
	barWidth=35;
	barDepth=17.25;
	wedgeWidthOff = 1.75;
	wedgeDepthOff = 10;
	wedgeFarAngle = .006*57.3;
	wedgeCloseAngle = 30;
	wedgeWidth=barWidth - wedgeWidthOff;
	wedgeDepthHigh = 79;
	wedgeHeight = 91;
	upperWedgeDepthHigh = 130;
	upperWedgeTop = 178.6;
	upperWedgeHeight = 78;
	upperWedgeBottom = upperWedgeTop-upperWedgeHeight;
	
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

	upperWedgeClosePlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeClosePlaneNy + upperWedgeClosePlaneNz*lowerWedgeExtensionZ;
	
	upperWedgeFarPlaneNx = 0; //Shouldn't be needed
	upperWedgeFarPlaneNy = 0;
	upperWedgeFarPlaneNz = -1;
	
	upperWedgeFarPlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeFarPlaneNy;
	
	focMirrorBottom = 139 + upperWedgeTop + barLength/2;
	fill_threeseg_plane_vecs();
	fill_foc_mirror_vecs();
	
	
	
	
	
	sidemirror_xl = -1000000;
	sidemirror_xr = 1000000;
	
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

	boxCloseZ = -614;
	
	reflOff = 9;
	
	three_seg_mirror = false;
	
	rand_gen = new TRandom3(rand_seed);
	
	qwartz = ref<qwartzCrystal>::create();
	oil = ref<mineralOil>::create();
	
	box_angle_off_cval = 1;
	box_angle_off_sval = 0;
	
	
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


	
	for (int i = 0; i < num_QE; i++)
	{
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
	for (int i = 0; i < num_transmittance; i++)
	{
		quartz_transmittance.push_back(tmp_quartz_transmittance[i]);
	}
	build_system();
}
void DircGopticalSim::fill_sens_plane_vecs()
{
	sensPlaneNx = 0;
	sensPlaneNy = sin(sens_rot/57.3);
	sensPlaneNz = cos(sens_rot/57.3);
	
	sensPlaneY = -sens_size*sin(sens_rot/57.3)/2-reflOff+barLength/2;
	sensPlaneZ = boxCloseZ + sens_size*cos(sens_rot/57.3)/2;
}
void DircGopticalSim::fill_foc_mirror_vecs()
{
	//is off by pi/2 to reduce rounding errors and such
	double foc_center_ang = foc_rot/57.3 + acos(-foc_mirror_size/(2*foc_r));
	printf("fca: %12.04f f_rot: %12.04f acos: %12.04f ratio: %12.04f\n",foc_center_ang*57.3,foc_rot,acos(-foc_mirror_size/(2*foc_r))*57.3);
	focMirrorY = focMirrorBottom - foc_r*cos(foc_center_ang);
	focMirrorZ = foc_r*sin(foc_center_ang);
	
	printf("focMirrorY: %12.04f foxMirrorZ: %12.04f\n",focMirrorY,focMirrorZ);
}
void DircGopticalSim::fill_threeseg_plane_vecs()
{
	focMirrorTop = focMirrorBottom + foc_mirror_size*cos(foc_rot/57.3);
	focMirrorZDim = foc_mirror_size*sin(foc_rot/57.3);
	printf("MirrorZDim: %12.04f MirrorSize: %12.04f\n",focMirrorZDim,foc_mirror_size);
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
  
// 	printf("ThreeSeg3z: %12.04f  seg_h: %12.04f\n",threeSeg3Z,seg_h);
	
	threeSeg1Ny = sin(theta_1);
	threeSeg1Nz = cos(theta_1);
	rotate_2d(threeSeg1Nx,threeSeg1Nz,cos(foc_yrot/57.3),sin(foc_yrot/57.3));//I think this is a slightly wrong rotation if the mirrors are carved out of a solid block, but it should be good enough at small angles
	threeSeg1D = threeSeg1Ny*threeSeg1Y + threeSeg1Nz*threeSeg1Z;//Use point x=0 as reference
		
	threeSeg2Ny = sin(theta_2);
	threeSeg2Nz = cos(theta_2);
	rotate_2d(threeSeg2Nx,threeSeg2Nz,cos(foc_yrot/57.3),sin(foc_yrot/57.3));
	threeSeg2D = threeSeg2Ny*threeSeg2Y + threeSeg2Nz*threeSeg2Z;//Use point x=0 as reference
	
	threeSeg3Ny = sin(theta_3);
	threeSeg3Nz = cos(theta_3);
	rotate_2d(threeSeg3Nx,threeSeg3Nz,cos(foc_yrot/57.3),sin(foc_yrot/57.3));
	threeSeg3D = threeSeg3Ny*threeSeg3Y + threeSeg3Nz*threeSeg3Z;//Use point x=0 as reference
	
	
}
void DircGopticalSim::sidemirror_reflect_points(std::vector<dirc_point> &points)
{
  double tmpx = 0;
  for (unsigned int i = 0; i < points.size(); i++)
  {
    tmpx = points[i].x;
    while (tmpx < sidemirror_xl || tmpx > sidemirror_xr)
    {
      if (tmpx < sidemirror_xl)
      {
	tmpx = 2*sidemirror_xl - tmpx;
      }
      if (tmpx > sidemirror_xr)
      {
	tmpx = 2*sidemirror_xr - tmpx;
      }
    }
    points[i].x = tmpx;
  }
}
void DircGopticalSim::set_sidemirror(double ixr, double ixl)
{
  sidemirror_xr = ixr;
  sidemirror_xl = ixl;
}
void DircGopticalSim::set_three_seg_mirror(bool itsm)
{
  if (three_seg_mirror != itsm)
  {
    three_seg_mirror = itsm;
    clear_system();
    build_system();
  }
}
void DircGopticalSim::set_pmt_offset(double r)
{
  //positive r increases path length
	clear_system();
    	boxCloseZ = -614 - r*cos(sens_rot/57.3);
	reflOff = 9 + r*sin(sens_rot/57.3);
	build_system();
}
void DircGopticalSim::set_liquid_absorbtion(double iabs)
{
	liquidAbsorbtion = iabs;
}
std::vector<double> DircGopticalSim::get_dist_traveled()
{
	return dist_traveled;
}
void DircGopticalSim::set_store_traveled(bool sst/* = true*/)
{
	store_traveled = sst;
}
void DircGopticalSim::set_liquid_index(double li)
{
	liquidIndex = li;
	double x, y, z;
	x = .5;
	y = 0;
	z = sqrt(3)/2;
	optical_interface_z(quartzIndex,liquidIndex,x,y,z);
	
// 	printf("n12 = %12.04f sin21 = %12.04f\n",quartzIndex/liquidIndex,2*(x/(sqrt(z*z+x*x))));
}
void DircGopticalSim::rotate_2d(double &x, double &y, double cval, double sval)
{
	//Standard rotatation allows precomputation of matrix elements
	//Only store one variable for speed
	//Sin should be actual sin (not negative sin)
	double tx = x;
	x = cval * tx - sval * y;
	y = -sval * tx + cval * y;
}
void DircGopticalSim::set_bar_box_angle(double ang)
{
	//expect radians
	box_angle_off_cval = cos(ang);
	box_angle_off_sval = sin(ang);
}
void DircGopticalSim::set_upper_wedge_angle_diff(double rads, double rads_y)
{
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
void DircGopticalSim::clear_system()
{
	Sys::Container::element_list_t elements = sys.get_element_list();
	Sys::Container::element_list_t::const_iterator iterator;
	for (iterator = elements.begin(); iterator != elements.end(); ++iterator)
	{
		sys.remove(**iterator);
	}
}
void DircGopticalSim::set_focus_mirror_angle(double ang,double yang)
{
	clear_system();
	foc_rot = ang;
	foc_yrot = yang;
	build_system();
}
void DircGopticalSim::set_pmt_angle(double ang)
{
	clear_system();
	sens_rot = ang;
	build_system();
}
void DircGopticalSim::set_wedge_mirror_rand(double ispread)
{
	if (ispread > 0)
	{
		upperWedgeNonUniform = true;
		upperWedgeNonUniformSpread = ispread;
	}
	else
	{
		upperWedgeNonUniform = false;
		upperWedgeNonUniformSpread = 0;
	}
	spread_wedge_mirror();
}
void DircGopticalSim::spread_wedge_mirror()
{
	
	if (upperWedgeNonUniform == true)
	{
		double spread_ang = wedgeCloseAngle;
		spread_ang += rand_gen->Gaus(0,upperWedgeNonUniformSpread);
		upperWedgeClosePlaneNx = 0; //Shouldn't be needed
		upperWedgeClosePlaneNy = sin(spread_ang/57.296);
		upperWedgeClosePlaneNz = cos(spread_ang/57.296);
		upperWedgeClosePlaneD = (barLength/2 + upperWedgeBottom)*upperWedgeClosePlaneNy + upperWedgeClosePlaneNz*lowerWedgeExtensionZ;
	}
}
	
void DircGopticalSim::build_system()
{
	//**********************************************************************
	// Optical system definition

// 	double size = 2200;//overall System scale 13.1 deg = 10+3.1 at 4m
	double size = 4400;
// 	size = 500;
	double testSize = 1800;

	double focR = foc_r;
	double focMirrorSize = foc_mirror_size;
	double focRot = foc_rot;
	double sensSize = sens_size;
	double sensRot = sens_rot;
	
	double testy = barLength/2+1000;
	
	
	
//Must add all elements as "ref<X>", otherwise they will descope and the system won't work after this function call
	ref<Curve::Flat> flatCurve = ref<Curve::Flat>::create();
	
	//Box
	//Focus Mirror
	
	double boxSize = size;
	double mirrorDepth = focMirrorZDim;//I would not really call this the mirror depth any more - more like the head on height
	double focYoff = focMirrorBottom;
	
	if (three_seg_mirror == false)
	{
	  focRot = foc_rot;
	  focMirrorSize = foc_mirror_size;
	  
	  double focz = -focMirrorSize*sin(focRot/57.3)/2;
	  double focy = focMirrorSize*cos(focRot/57.3)/2;
	  
	  
	  
	  double fx,fy,fz;
	  fx = 0;
	  fy = focy + focYoff;
	  fz = focz;
	  
	  ref<Shape::Rectangle> focMirrorShape = ref<Shape::Rectangle>::create(boxSize,focMirrorSize);
	  ref<cylinderY> focCylinder = ref<cylinderY>::create(-focR);
	  ref<Sys::Mirror> focMirror = \
		  ref<Sys::Mirror>::create(\
			  Math::Vector<3>(fx,fy,fz),\
			  focCylinder,\
			  focMirrorShape,\
			  true,\
			  Material::mirror,\
			  oil);
	  focMirror->rotate(focRot,0,0);
	  focMirror->rotate(0,foc_yrot,0);
// 	  sys.add(*focMirror);
	}
	else if (three_seg_mirror == true)//condition not needed obviously
	{
	  double fx,fy,fz; //serve as holding variables for the 3 mirror positions;
	  fx = 0;
	  fy = 0;
	  fz = 0;
	  
	  double theta_m = foc_rot/57.3;//radians of rotation to go through
	  double theta_c = fabs(2*asin(mirrorDepth/(2*focR)));//angle subtended by the mirror
	  double seg_h = fabs(2*focR*sin(theta_c/6));//length of each segment;
	  
	  
	  
	  ref<Shape::Rectangle> segMirrorShape = ref<Shape::Rectangle>::create(boxSize,seg_h);
	  ref<Curve::Flat> segMirrorCurve = ref<Curve::Flat>::create();
	  
	  //I had to do some geometry and algebra to get these numbers, but in hindsight, it's obvious.  Always that way.
	  double theta_1 = theta_m - theta_c/3;
	  double theta_2 = theta_m;
	  double theta_3 = theta_m + theta_c/3;
	  
	  fy = focYoff + seg_h*cos(theta_1)/2;
	  fz = -seg_h*sin(theta_1)/2;
	  
// 	  printf("tm = %12.04f tc = %12.04f fy = %12.04f seg_h = %12.04f focy = %12.04f\n",theta_m,theta_c,fy,seg_h,seg_h*cos(theta_1)/2);
	  
	  ref<Sys::Mirror> segMirror1 = \
		  ref<Sys::Mirror>::create(\
			  Math::Vector<3>(fx,fy,fz),\
			  segMirrorCurve,\
			  segMirrorShape,\
			  true,\
			  Material::mirror,\
			  oil);
	  
	  segMirror1->rotate(0,foc_yrot,0);
	  segMirror1->rotate(theta_1*57.3,0,0);
	  
	  
// 	  sys.add(*segMirror1);

	  fy += seg_h*cos(theta_1)/2 + seg_h*cos(theta_2)/2;//update these by moving the centers
	  fz += -seg_h*sin(theta_1)/2 - seg_h*sin(theta_2)/2;
	  
	  ref<Sys::Mirror> segMirror2 = \
		  ref<Sys::Mirror>::create(\
			  Math::Vector<3>(fx,fy,fz),\
			  segMirrorCurve,\
			  segMirrorShape,\
			  true,\
			  Material::mirror,\
			  oil);
	  
	  segMirror2->rotate(0,foc_yrot,0);
	  segMirror2->rotate(theta_2*57.3,0,0);
	  
	  
// 	  sys.add(*segMirror2);

	  fy += seg_h*cos(theta_2)/2 + seg_h*cos(theta_3)/2;
	  fz += -seg_h*sin(theta_2)/2 - seg_h*sin(theta_3)/2;
	  
	  ref<Sys::Mirror> segMirror3 = \
		  ref<Sys::Mirror>::create(\
			  Math::Vector<3>(fx,fy,fz),\
			  segMirrorCurve,\
			  segMirrorShape,\
			  true,\
			  Material::mirror,\
			  oil);
	  
	  segMirror3->rotate(0,foc_yrot,0);
	  segMirror3->rotate(theta_3*57.3,0,0);
	  
// 	  sys.add(*segMirror3);
	  
	}
	//Image Plane
	
	sensSize = 312;
	//sensSize = 600;
	sensRot = sens_rot;

	double backFlatYSize = sensSize/2;
	
	double sx,sy,sz;
	sx = 0;
	sy = -sensSize*sin(sensRot/57.3)/2-reflOff+barLength/2;
	sz = boxCloseZ + sensSize*cos(sensRot/57.3)/2;
	
	ref<Shape::Rectangle> boxImageShape = ref<Shape::Rectangle>::create(boxSize,3*sensSize);
	ref<Curve::Flat> boxImageCurve = ref<Curve::Flat>::create();
// 	declared in header for calling later
	boxImage = \
		ref<Sys::Image>::create(\
			Math::Vector<3>(sx,sy,sz),\
			boxImageCurve,\
			boxImageShape);
	boxImage->rotate(sensRot,0,0);
	sys.add(*boxImage);
	
	
	//Below is how we'd do it if there was no focusing mirror
	//500 is the 2.5mrad point at the normal parameters
// 	sx = 0;
// 	sy = barLength/2 + upperWedgeTop + 500;
// 	sz = 0;
// 	
// 	ref<Shape::Rectangle> boxImageShape = ref<Shape::Rectangle>::create(boxSize,8*sensSize);
// 	ref<Curve::Flat> boxImageCurve = ref<Curve::Flat>::create();
// 	boxImage = \
// 		ref<Sys::Image>::create(\
// 			Math::Vector<3>(sx,sy,sz),\
// 			boxImageCurve,\
// 			boxImageShape);
// 	boxImage->rotate(90,0,0);
// 	sys.add(*boxImage);
// 	
//End no mirror code
	
	ref<Shape::Rectangle> blockShape = ref<Shape::Rectangle>::create(boxSize,backFlatYSize);
	ref<Sys::Mirror> backFlatMirror = \
		ref<Sys::Mirror>::create(\
			Math::Vector<3>(0,barLength/2+upperWedgeTop + backFlatYSize/2,0),\
			flatCurve,\
			blockShape,\
			true,\
			Material::mirror,\
			oil);
// 	sys.add(*backFlatMirror);
	
	ref<Curve::Flat> imgCurve = ref<Curve::Flat>::create();
	ref<Shape::Rectangle> imgShape = ref<Shape::Rectangle>::create(size/0.5,testSize*2);
	ref<Sys::Image> testImage = \
		ref<Sys::Image>::create(\
			Math::VectorPair<3>(0,testy,0,0,0,0),\
			imgCurve,\
			imgShape);
	testImage->rotate(90,0,0);
// 	sys.add(*testImage);
	
}
double DircGopticalSim::get_quartz_n(double lambda)
{
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
double DircGopticalSim::get_cerenkov_angle_rand(double beta, double additional_spread, double &wavelength)
{
// 	printf("called get ckov\n");
	//May be slow enough to consider approximating in distribution generation
	double out_ang = 0;
	double tmp_lam = 0;
	double tmp_QE_val;
	double above_ind;
	int ind_QE;
	double n_lam;
    
	while (true)
	{
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
bool DircGopticalSim::quartz_transmission_mc(double R, double lambda)
{
	if (R > 10000)
	{
// 	  printf("using R: %12.04f\n",R);
	}
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
	
// 	printf("trans_prob: %12.04f\n",trans_prob);
	
	return (rand_gen->Uniform(0,1) < trans_prob);

}
double DircGopticalSim::get_beta(double E, double m)
{
	double gam = E/m;
	if (gam > 5)
	{
		//Approximate and possibly save 
		return 1-1/(2*gam*gam);
	}
	else
	{
		return sqrt(1-1/(gam*gam));
	}
	
}
std::vector<dirc_point> DircGopticalSim::sim_rand_n_photons(\
	int n_photons, \
	bool outspot /*= false*/, \
	bool outframe /*= false*/, \
	double ckov_theta /*= 47*/, \
	double particle_x /*= 0*/, \
	double particle_y /*= 0*/, \
	double particle_theta /*= 0*/, \
	double particle_phi /*= 0*/,\
	double phi_theta_unc /*= .0015*57.3*/,\
	double ckov_theta_unc /* = .0055*57.3*/,\
	double beta /* = -1*/,\
	bool check_dir /*= false*/)
{
	if (check_dir == false)
	{
		Sys::SourceRays srcrays(Math::Vector<3>(0,0,0));
		
		fill_rand_phi(\
			srcrays,\
			n_photons,\
			ckov_theta,\
			particle_x,\
			particle_y,\
			particle_theta,\
			particle_phi,\
			phi_theta_unc,\
			ckov_theta_unc,\
			0,\
			beta);
		return trace_source_rays(srcrays,outspot,outframe);
	}
	else
	{
		//Slow way for now, faster once goptical is gone
		std::vector<dirc_point> rval;
		std::vector<dirc_point> tval;
		
		Sys::SourceRays upsrcrays(Math::Vector<3>(0,0,0));
		Sys::SourceRays downsrcrays(Math::Vector<3>(0,0,0));
		fill_rand_phi(\
			upsrcrays,\
			n_photons,\
			ckov_theta,\
			particle_x,\
			particle_y,\
			particle_theta,\
			particle_phi,\
			phi_theta_unc,\
			ckov_theta_unc,\
			1,\
			beta);
		fill_rand_phi(\
			downsrcrays,\
			n_photons,\
			ckov_theta,\
			particle_x,\
			particle_y,\
			particle_theta,\
			particle_phi,\
			phi_theta_unc,\
			ckov_theta_unc,\
			-1,\
			beta);
		
		tval = trace_source_rays(upsrcrays,outspot,outframe);
		for (unsigned int i = 0; i < tval.size(); i++)
		{
			dirc_point tpoint;
			tpoint.x = tval[i].x;
			tpoint.y = tval[i].y;
			tpoint.t = 1;
			tpoint.weight = tval[i].weight;
			rval.push_back(tpoint);
		}
		tval = trace_source_rays(downsrcrays,outspot,outframe);
		for (unsigned int i = 0; i < tval.size(); i++)
		{
			dirc_point tpoint;
			tpoint.x = tval[i].x;
			tpoint.y = tval[i].y;
			tpoint.t = -1;
			tpoint.weight = tval[i].weight;
			rval.push_back(tpoint);
		}
		return rval;
	}

}
std::vector<dirc_point> DircGopticalSim::sim_reg_n_photons(\
	int n_photons_phi, \
	int n_photons_z,\
	bool outspot /*= false*/, \
	bool outframe /*= false*/, \
	double ckov_theta /*= 47*/, \
	double particle_x /*= 0*/, \
	double particle_y /*= 0*/, \
	double particle_theta /*= 0*/, \
	double particle_phi /*= 0*/,\
	double phi_theta_unc /*= 0*/,\
	double ckov_theta_unc /* = 0*/,\
	double beta /* = -1*/,\
	bool check_dir /*= false*/)
{
	if (check_dir == false)
	{
		Sys::SourceRays srcrays(Math::Vector<3>(0,0,0));
		fill_reg_phi(\
			srcrays,\
			n_photons_phi,\
			n_photons_z,\
			ckov_theta,\
			particle_x,\
			particle_y,\
			particle_theta,\
			particle_phi,\
			phi_theta_unc,\
			ckov_theta_unc,\
			0,\
			beta);
	
		
		return trace_source_rays(srcrays,outspot,outframe);
	}
	else
	{
		//Slow way for now, faster once goptical is gone
		std::vector<dirc_point> rval;
		std::vector<dirc_point> tval;
		
		Sys::SourceRays upsrcrays(Math::Vector<3>(0,0,0));
		Sys::SourceRays downsrcrays(Math::Vector<3>(0,0,0));
		fill_reg_phi(\
			upsrcrays,\
			n_photons_phi,\
			n_photons_z,\
			ckov_theta,\
			particle_x,\
			particle_y,\
			particle_theta,\
			particle_phi,\
			phi_theta_unc,\
			ckov_theta_unc,\
			1,\
			beta);
		fill_reg_phi(\
			downsrcrays,\
			n_photons_phi,\
			n_photons_z,\
			ckov_theta,\
			particle_x,\
			particle_y,\
			particle_theta,\
			particle_phi,\
			phi_theta_unc,\
			ckov_theta_unc,\
			-1,\
			beta);
		
		
		
		tval = trace_source_rays(upsrcrays,outspot,outframe);
		for (unsigned int i = 0; i < tval.size(); i++)
		{
			dirc_point tpoint;
			tpoint.x = tval[i].x;
			tpoint.y = tval[i].y;
			tpoint.t = 1;
			tpoint.weight = tval[i].weight;
			rval.push_back(tpoint);
		}
		tval = trace_source_rays(downsrcrays,outspot,outframe);
		for (unsigned int i = 0; i < tval.size(); i++)
		{
			dirc_point tpoint;
			tpoint.x = tval[i].x;
			tpoint.y = tval[i].y;
			tpoint.t = -1;
			tpoint.weight = tval[i].weight;
			rval.push_back(tpoint);
		}
		return rval;
	}
		
}
void DircGopticalSim::fill_rand_phi(\
	Sys::SourceRays &srcrays,\
	int n_photons, \
	double ckov_theta /*= 47*/, \
	double particle_x /*= 0*/, \
	double particle_y /*= 0*/, \
	double particle_theta /*= 0*/, \
	double particle_phi /*= 0*/,\
	double phi_theta_unc, /*= .0015*57.3*/
	double ckov_theta_unc /* = .0055*57.3*/,\
	double check_dir /* = 0*/,\
	double beta/* = -1*/)
{
// 	double sDepth = .95*barDepth;
	double emitAngle = ckov_theta;
	double particleTheta = particle_theta + rand_gen->Gaus(0,phi_theta_unc);
	double particlePhi = particle_phi + rand_gen->Gaus(0,phi_theta_unc);
	
	int numPhots = n_photons;
	
// 	double sourcez = -sDepth;
	double sourcey = particle_y-barDepth*tan(particleTheta/57.3);
	double sourcex = particle_x;
	double tsy = sourcey;
	double tsx = sourcex;
	sourcey = tsy*cos(particlePhi/57.3)-tsx*sin(particlePhi/57.3);
	sourcex = tsy*sin(particlePhi/57.3)+tsx*cos(particlePhi/57.3);
	
// 	srcrays.set_position(Math::Vector<3>(sourcex,sourcey,sourcez));
	//Need to hand absolute positions to warp_ray
	srcrays.set_position(Math::Vector<3>(0,0,0));
	srcrays.set_material(oil);
	
	double sourceOff,randPhi;
	
	Math::Transform3 trans;
	trans.reset();
	trans.linear_rotation(Math::Vector3(particleTheta,0,0));
	trans.linear_rotation(Math::Vector3(0,0,particlePhi));
	
	double temit, rand_add;
	double wavelength = 400;
	
	for (int i = 0; i < numPhots; i++)
	{ 
		randPhi = rand_gen->Uniform(0,2*3.14159265);
		sourceOff = -rand_gen->Uniform(0,barDepth);

		if (beta < 0)
		{
			rand_add = rand_gen->Gaus(0,ckov_theta_unc);
			temit = emitAngle + rand_add;
		}
		else
		{
			temit = get_cerenkov_angle_rand(beta,ckov_theta_unc,wavelength);
		}
		Math::Vector3 start_ray(0,0,sourceOff);
		Math::Vector3 dir_ray(\
			sin(temit/57.3)*cos(randPhi),\
			sin(temit/57.3)*sin(randPhi),\
			cos(temit/57.3));
		
		start_ray = trans.transform_linear(start_ray);
		dir_ray = trans.transform_linear(dir_ray);	
		
		start_ray.x() = start_ray.x() + particle_x;
		start_ray.y() = start_ray.y() + particle_y;
		
		Math::VectorPair3 new_ray(\
			start_ray.x(),\
			start_ray.y(),\
			start_ray.z(),\
			dir_ray.x(),\
			dir_ray.y(),\
			dir_ray.z());
		
		if (dir_ray.y() * check_dir < -.00001)
		{
			//Get all rays in one direction or the other.
			//check_dir = 0 sidesteps this
			continue;
		}
		
		
		warp_ray(\
			new_ray.x0(),\
			new_ray.y0(),\
			new_ray.z0(),\
			new_ray.x1(),\
			new_ray.y1(),\
			new_ray.z1(),\
			asin(1/1.47));
		if (new_ray.z0() > 0)
		{
			continue;
		}
		
		spread_wedge_mirror();
			
		warp_wedge(\
			new_ray.x0(),\
			new_ray.y0(),\
			new_ray.z0(),\
			new_ray.x1(),\
			new_ray.y1(),\
			new_ray.z1());
		
		//account (quickly) for the bar box having a different angle than the readout
		rotate_2d(new_ray.y1(),new_ray.z1(),box_angle_off_cval,box_angle_off_sval);
		
		if (new_ray.z0() > 0)
		{
			continue;
		}
		
		//check absorbtion
		if (!(absorbtion_mc(new_ray.x1(),new_ray.y1())))
		{
			continue;
		}
		
// 			printf("%d\n",num_through++);
		
		srcrays.add_rays(new_ray,&srcrays);
	}
}
std::vector<std::pair<double,double> > DircGopticalSim::get_refraction_rand_phi(\
	std::vector<double> &before_interface,\
	std::vector<double> &after_interface,\
	int n_photons, \
	double ckov_theta /*= 47*/, \
	double particle_x /*= 0*/, \
	double particle_y /*= 0*/, \
	double particle_theta /*= 0*/, \
	double particle_phi /*= 0*/,\
	double phi_theta_unc, /*= .0015*57.3*/
	double ckov_theta_unc /* = .0055*57.3*/,\
	double check_dir /* = 0*/,\
	double beta/* = -1*/)
{
	//returns theta versus cerenkov phi_theta_unc
  
	std::vector<std::pair<double,double> > rval;
  
	before_interface.clear();
	after_interface.clear();
	refraction_before.clear();
	refraction_after.clear();
	store_refraction = true;
	
// 	printf("refraction called with n: %d\n",n_photons);
	
	double emitAngle = ckov_theta;
	double particleTheta = particle_theta + rand_gen->Gaus(0,phi_theta_unc);
	double particlePhi = particle_phi + rand_gen->Gaus(0,phi_theta_unc);
	
	int numPhots = n_photons;
	
// 	double sourcez = -sDepth;
	double sourcey = particle_y-barDepth*tan(particleTheta/57.3);
	double sourcex = particle_x;
	double tsy = sourcey;
	double tsx = sourcex;
	sourcey = tsy*cos(particlePhi/57.3)-tsx*sin(particlePhi/57.3);
	sourcex = tsy*sin(particlePhi/57.3)+tsx*cos(particlePhi/57.3);
	
	
	double bar_dist_travelled;
	
	double sourceOff,randPhi;
	
	Math::Transform3 trans;
	trans.reset();
	trans.linear_rotation(Math::Vector3(particleTheta,0,0));
	trans.linear_rotation(Math::Vector3(0,0,particlePhi));
	
	double temit, rand_add;
	double wavelength = 400;
	
	for (int i = 0; i < numPhots; i++)
	{ 
		randPhi = rand_gen->Uniform(0,2*3.14159265);
		sourceOff = -rand_gen->Uniform(0,barDepth);
		

		if (beta < 0)
		{
			rand_add = rand_gen->Gaus(0,ckov_theta_unc);
			temit = emitAngle + rand_add;
		}
		else
		{
			temit = get_cerenkov_angle_rand(beta,ckov_theta_unc,wavelength);
		}
		Math::Vector3 start_ray(0,0,sourceOff);
		Math::Vector3 dir_ray(\
			sin(temit/57.3)*cos(randPhi),\
			sin(temit/57.3)*sin(randPhi),\
			cos(temit/57.3));
		
		start_ray = trans.transform_linear(start_ray);
		dir_ray = trans.transform_linear(dir_ray);	
		
		start_ray.x() = start_ray.x() + particle_x;
		start_ray.y() = start_ray.y() + particle_y;
		
		Math::VectorPair3 new_ray(\
			start_ray.x(),\
			start_ray.y(),\
			start_ray.z(),\
			dir_ray.x(),\
			dir_ray.y(),\
			dir_ray.z());
		
		if (dir_ray.y() * check_dir < -.00001)
		{
			//Get all rays in one direction or the other.
			//check_dir = 0 sidesteps this
			continue;
		}
		

		
		bar_dist_travelled = warp_ray(\
			new_ray.x0(),\
			new_ray.y0(),\
			new_ray.z0(),\
			new_ray.x1(),\
			new_ray.y1(),\
			new_ray.z1(),\
			asin(1/1.47));
		if (new_ray.z0() > 0)
		{
			continue;
		}
		
// 		if (bar_dist_travelled > 70000) continue;
		
		if (!quartz_transmission_mc(bar_dist_travelled,wavelength))
		{
		  //absorbed in bar
		  continue;
		}
		
		
		spread_wedge_mirror();
		
		
		warp_wedge(\
			new_ray.x0(),\
			new_ray.y0(),\
			new_ray.z0(),\
			new_ray.x1(),\
			new_ray.y1(),\
			new_ray.z1());
		
		rotate_2d(new_ray.y1(),new_ray.z1(),box_angle_off_cval,box_angle_off_sval);
		
		if (new_ray.z0() > 0)
		{
			continue;
		}
		
		warp_box(\
			new_ray.x0(),\
			new_ray.y0(),\
			new_ray.z0(),\
			new_ray.x1(),\
			new_ray.y1(),\
			new_ray.z1());
		
		if (new_ray.z0() > 0)
		{
			continue;
		}
		
		std::pair<double, double> add_val;
// 		printf("adding: %12.04f %12.04f\n",randPhi,refraction_before[refraction_before.size() - 1]);
		add_val.first = randPhi;
		add_val.second = refraction_before[refraction_before.size() - 1];
		rval.push_back(add_val);
	}
	
	for (unsigned int i = 0; i < refraction_before.size(); i++)
	{
      // 	printf("filling vector %d\n",i);
		before_interface.push_back(refraction_before[i]);
		after_interface.push_back(refraction_after[i]);
	}
	
	store_refraction = false;
	refraction_before.clear();
	refraction_after.clear();
	
	return rval;
}

void DircGopticalSim::fill_reg_phi(\
	Sys::SourceRays &srcrays,\
	int n_photons_phi, \
	int n_photons_z,\
	double ckov_theta /*= 47*/, \
	double particle_x /*= 0*/, \
	double particle_y /*= 0*/, \
	double particle_theta /*= 0*/, \
	double particle_phi /*= 0*/,\
	double phi_theta_unc /*= 0*/,\
	double ckov_theta_unc /* = 0*/,\
	double check_dir /* =0 */,\
	double beta /* = -1*/)
{
	double sDepth = .95*barDepth;
	double emitAngle = ckov_theta;
	double particleTheta = particle_theta;
	double particlePhi = particle_phi;
	
// 	double sourcez = -sDepth;
	double sourcey = particle_y-barDepth*tan(particleTheta/57.3);
	double sourcex = particle_x;
	double tsy = sourcey;
	double tsx = sourcex;
	sourcey = tsy*cos(particlePhi/57.3)-tsx*sin(particlePhi/57.3);
	sourcex = tsy*sin(particlePhi/57.3)+tsx*cos(particlePhi/57.3);
	
// 	srcrays.set_position(Math::Vector<3>(sourcex,sourcey,sourcez));
	//Need to hand absolute positions to warp_ray
	srcrays.set_position(Math::Vector<3>(0,0,0));
	srcrays.set_material(oil);
	
	double sourceOff,regPhi;
	
	Math::Transform3 trans;
	trans.reset();
	trans.linear_rotation(Math::Vector3(particle_theta,0,0));
	trans.linear_rotation(Math::Vector3(0,0,particle_phi));

	int num_through = 0;
	double temit;
	double rand_add;
	double wavelength = 400;
	for (int i = 0; i < n_photons_z; i++)
	{
		sourceOff = (i+.5)*sDepth/(n_photons_z);
		
		for (int j = 0; j < n_photons_phi; j++)
		{ 
			regPhi = j*2*3.14159265357/(n_photons_phi);
			
			if (beta < 0)
			{
				rand_add = rand_gen->Gaus(0,ckov_theta_unc);
				temit = emitAngle + rand_add;
			}
			else
			{
				temit = get_cerenkov_angle_rand(beta,ckov_theta_unc,wavelength);
			}
				
			
			
			Math::Vector3 start_ray(0,0,sourceOff);
			Math::Vector3 dir_ray(\
				sin(temit/57.3)*cos(regPhi),\
				sin(temit/57.3)*sin(regPhi),\
				cos(temit/57.3));
			
			start_ray = trans.transform_linear(start_ray);
			start_ray.z() -= barDepth;
			dir_ray = trans.transform_linear(dir_ray);	
			
			start_ray.x() = start_ray.x() + particle_x;
			start_ray.y() = start_ray.y() + particle_y;
			
			Math::VectorPair3 new_ray(\
				start_ray.x(),\
				start_ray.y(),\
				start_ray.z(),\
				dir_ray.x(),\
				dir_ray.y(),\
				dir_ray.z());
			
			if (dir_ray.y() * check_dir < -.00001)
			{
			//Get all rays in one direction or the other.
			//check_dir = 0 sidesteps this
				continue;
			}
			
			warp_ray(\
				new_ray.x0(),\
				new_ray.y0(),\
				new_ray.z0(),\
				new_ray.x1(),\
				new_ray.y1(),\
				new_ray.z1(),\
				asin(1/1.47));
			
			if (new_ray.z0() > 0)
			{
				continue;
			}
			
			spread_wedge_mirror();
			
			warp_wedge(\
				new_ray.x0(),\
				new_ray.y0(),\
				new_ray.z0(),\
				new_ray.x1(),\
				new_ray.y1(),\
				new_ray.z1());
			
// 			printf("dz: %12.04f\n",new_ray.z1());
			
			//account (quickly) for the bar box having a different angle than the readout
			rotate_2d(new_ray.y1(),new_ray.z1(),box_angle_off_cval,box_angle_off_sval);
			
// 			printf("%12.04f %12.04f %12.04f :: %12.04f %12.04f %12.04f\n",new_ray.x0(),new_ray.y0(),new_ray.z0(),new_ray.x1(),new_ray.y1(),new_ray.z1());
			
			if (new_ray.z0() > 0)
			{
				continue;
			}
			
			warp_box(\
				new_ray.x0(),\
				new_ray.y0(),\
				new_ray.z0(),\
				new_ray.x1(),\
				new_ray.y1(),\
				new_ray.z1());
			
			if (new_ray.z0() > 0)
			{
				continue;
			}
			
			//check absorbtion
			if (!(absorbtion_mc(new_ray.x1(),new_ray.y1())))
			{
				continue;
			}
// 			printf("%d\n",num_through++);
			srcrays.add_rays(new_ray,&srcrays);
		}
	}
// 	srcrays.rotate(particle_theta,0,0);
// 	srcrays.rotate(0,0,particle_phi);
}
std::vector<dirc_point> DircGopticalSim::trace_source_rays(\
	Sys::SourceRays &srcrays, \
	bool outspot,\
	bool outframe)
{
	sys.add(srcrays);

	Trace::Tracer         tracer(sys);
	tracer.get_params().set_max_bounce(5000);
	tracer.get_trace_result().set_generated_save_state(srcrays);
	tracer.get_trace_result().set_intercepted_save_state(*boxImage);
	
	tracer.trace();
	std::vector<dirc_point> rval;
	
	_Goptical::Trace::rays_queue_t rays = tracer.get_trace_result().get_intercepted(*boxImage);
	
	double x, y;
	for (unsigned int i = 0; i < rays.size(); i++)
	{
		x = ((Math::Vector3)rays[i]->get_intercept_point()).x();
		y = ((Math::Vector3)rays[i]->get_intercept_point()).y();
		
		dirc_point add_point;
		add_point.x = x;
		add_point.y = y;
		add_point.weight = 1;
		
		rval.push_back(add_point);
	}
	
	if (outframe == true)
	{
		printf("Rendering 3d\n");
		
		Io::RendererSvg       svg_renderer("layout.svg", 1000, 700);
		Io::RendererViewport  &renderer = svg_renderer;
		// 3d system layout
		
		renderer.set_feature_size(20);
		
		renderer.set_perspective();
		renderer.set_fov(45);	
		
		
		sys.draw_3d_fit(renderer, 0);
		renderer.set_camera_transform(renderer.get_camera_transform().linear_rotation(Math::Vector3(0,0,0)));
		sys.draw_3d(renderer);
		
		tracer.get_trace_result().draw_3d(renderer);
		
		printf("rendered 3d\n");
	}
	if (outspot == true)
	{
		printf("Rendering Spot\n");
		Analysis::Spot spot(sys);
		spot.get_tracer().get_params().set_max_bounce(5000);
		Io::RendererSvg renderSpot("spot.svg", 1000, 1000, Io::rgb_black);
		
		spot.draw_diagram(renderSpot);
		printf("Rendered spot\n");
	}
	
	sidemirror_reflect_points(rval);
	sys.remove(srcrays);
	
	return rval;
}

// void DircGopticalSim::warp_ray(\
// 	Math::VectorPair3 &ray,\
// 	double critical_angle)
double DircGopticalSim::warp_ray(\
	double &x,\
	double &y,\
	double &z,\
	double &dx,\
	double &dy,\
	double &dz,\
	double critical_angle)
{
//Implemented to avoid total internal reflection computations
//Expects ray to be "prerotated" by particle theta and phi - just propagates the vectors
//Uses bar geometry.
//Modifys the ray.
//Be careful about x,y,and z - can't straight rotate that
//returns distance traveled
  
	double rval = 0; //stores total distance;
  
	double delx, dely,delz;
// 	double xzR;
	//Start the rays in the quartz for simulation reasons
	double grace_room = 0;
	
	if (dy > 0)
	{
		//Going up
		dely = barLength*(0.5-grace_room) - y;
	}
	else
	{
		//going down, flip dy
		dely = barLength*(1.5-grace_room) + y;
		dy = -dy;
	}
	
	rval += dely*dely;
	
	y = barLength*(.5-grace_room);
	
	if (acos(fabs(dz)) < critical_angle ||\
		acos(fabs(dx)) < critical_angle ||\
		(dy < 0 && acos(fabs(dy)) < critical_angle) ||
		dy*dy < 1e-4)
	{
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
	if (dx < 0)
	{
		lrx = -1;
	}
	else
	{
		lrx = 1;
	}
	
	delx = fabs(dely*dx/dy) - lrx*(barWidth/2-x);
	delz = fabs(dely*dz/dy) + z;
	
	rval += delx*delx + delx*delx;
	
	if (dz < 0) printf("dz: %12.04f : Negative and confusing\n",dz);

	
	nbouncesx = delx/barWidth + 1;
	remainderx = delx - barWidth*(nbouncesx-1);
	
	if (nbouncesx % 2 == 1)
	{
		dx = -dx;
		x = lrx*(barWidth/2-remainderx);
	}
	else
	{
		x = lrx*(-barWidth/2+remainderx);
	}
	
	if (x > barWidth/2 - wedgeWidthOff)
	{
		//Hit the 5mm where the wedge is not.  Assume lost
		z = 1337;
		return -1;
	}
	nbouncesz = delz/barDepth + 1;
	remainderz = delz - barDepth*(nbouncesz-1);
	
	if (nbouncesz % 2 == 1)
	{
		dz = -dz;
		z = -remainderz;
	}
	else
	{
		z = remainderz-barDepth;
	}
	
	return sqrt(rval);
}
void DircGopticalSim::warp_wedge(\
	double &x,\
	double &y,\
	double &z,\
	double &dx,\
	double &dy,\
	double &dz)
{
	//No Critical angle checking - shouldn't need it
	//I am abusing the x angle being the same in the bar and wedge here
	//may have to rotate after to the mirror
	//Starts y at 0
	
	bool passed_interface = false;
// 	double x0 = x;
	double dt = 0; //dummy variable representing parametric line
	double n_dot_v = 0; //dummy variable to save dot product
	double n_dot_v0 = 0;
	
	//deal with yz first
	
	//Check for reflection from far wedge plane - max 1 bounce
	if (dz > 0)
	{
		n_dot_v = -(dy*wedgeFarPlaneNy + dz*wedgeFarPlaneNz);
		n_dot_v0 = -(y*wedgeFarPlaneNy + z*wedgeFarPlaneNz);
		
		dt = -(wedgeFarPlaneD+n_dot_v0)/n_dot_v;
//  		printf("lw: %12.04f %12.04f %12.04f :: %12.04f %12.04f %12.04f dt: %12.04f\n",x,y,z,dx,dy,dz,dt);
		//pretending y=0 at bar top for speed
		if (dt*dy < wedgeHeight)
		{
			//reflect it off bottom wedge
			//Does not pass through optical interface
			
			//Will always be true - just use for propagation (short circuit var?)
			
			x_wedge_coerce_check(x,y,z,dx,dy,dz,dt);
			
			
			dy += 2*n_dot_v*wedgeFarPlaneNy;
			dz += 2*n_dot_v*wedgeFarPlaneNz;
		}
		else
		{
// 			printf("IM HERE REFLECTING AND SUCH\n");
// 			n_dot_v = -(dy*upperWedgeFarPlaneNy + dz*upperWedgeFarPlaneNz);
// 			n_dot_v0 = -(y*upperWedgeFarPlaneNy + z*upperWedgeFarPlaneNz);
// 		
// 			dt = -(upperWedgeFarPlaneD+n_dot_v0)/n_dot_v;
// 			//never reached??? probably not  can possibly remove if statement
// 			//No bottom wedge, try top
// 			if (dt*dy < upperWedgeTop)//Should always be true... I hope (remove later?)
// 			{
// 				//Does pass through optical interface
// // 				if (dt*dy
// 				
// 				//Following statement performs the propagation if it does not fail
// 				if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt)))
// 				{
// 					//Goes out the window, fail and return
// 					z = 1337;
// 					return;
// 				}
// 				
// 				//Reflect off top wedge			
// 				dy += 2*n_dot_v*upperWedgeFarPlaneNy;
// 				dz += 2*n_dot_v*upperWedgeFarPlaneNz;
// 			}
		}
	}
	//Now dz < 0 or we have a new starting vector.  Either way, we intersect with the "close" wedge now
	n_dot_v = -(dy*wedgeClosePlaneNy + dz*wedgeClosePlaneNz);
	n_dot_v0 = -(y*wedgeClosePlaneNy + z*wedgeClosePlaneNz);
	
	dt = -(wedgeClosePlaneD+n_dot_v0)/n_dot_v;
	
// 	printf("uw: %12.04f %12.04f %12.04f :: %12.04f %12.04f %12.04f dt: %12.04f\n",x,y,z,dx,dy,dz,dt);

	//Assume close enough to determine which wedge it hit
	double ty = dt*dy + y - barLength/2;
	
	if ((ty < wedgeHeight))
	{
		//reflect it off bottom wedge
		
		//Again, always true
		if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt)))
		{
			//Goes out the window, fail and return
			z = 1337;
			return;
		}
		
		dy += 2*n_dot_v*wedgeClosePlaneNy;
		dz += 2*n_dot_v*wedgeClosePlaneNz;
		
	}
	else if (ty > upperWedgeBottom && ty < upperWedgeTop)
	{
		//passed interface before reflection
		//get dt to propagate to interface
		dt = (quartzLiquidY + barLength/2 - y)/dy;
		if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt)))
		{
			//Goes out the window, fail and return
			z = 1337;
			return;
		}
		//correct ordering below - interface is xz plane, so dx and dz go first
		optical_interface_z(quartzIndex,liquidIndex,dx,dz,dy);
		passed_interface = true;
		
		//NOW intersect with upper wedge
		
		//Signs are weird here, it works, but I think they are twisted up with the initialization too
		n_dot_v = -(dy*upperWedgeClosePlaneNy + dz*upperWedgeClosePlaneNz + dx*upperWedgeClosePlaneNx);
		n_dot_v0 = -(y*upperWedgeClosePlaneNy + z*upperWedgeClosePlaneNz + x*upperWedgeClosePlaneNx);
	
		dt = -(upperWedgeClosePlaneD+n_dot_v0)/n_dot_v;
		
		if (dt*dy + y < upperWedgeTop + barLength/2)
		{
			if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt)))
			{
				//Goes out the window, fail and return
				z = 1337;
				return;
			}
			
			dx += 2*n_dot_v*upperWedgeClosePlaneNx;
			dy += 2*n_dot_v*upperWedgeClosePlaneNy;
			dz += 2*n_dot_v*upperWedgeClosePlaneNz;
			
		}
		else
		{
			//refracted such that it no longer bounces
			dt = (upperWedgeTop + barLength/2 - y)/dy;
			if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt)))
			{
				//Goes out the window, fail and return
				z = 1337;
				return;
			}
		}

	}
	else if (ty < upperWedgeTop)
	{
		//out the window
		z = 1337;
		return;
	}

// 	printf("aw: %12.04f %12.04f %12.04f :: %12.04f %12.04f %12.04f dt: %12.04f\n",x,y,z,dx,dy,dz,dt);
	
	//Have now performed the maximum number of reflections (besides x direction)
	//Finish by taking it to the top
	if (!passed_interface == true)
	{
		//Will go through the interface now before it hits the top
		dt = (quartzLiquidY + barLength/2 - y)/dy;
		if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt)))
		{
			//Goes out the window, fail and return
			z = 1337;
			return;
		}
		//correct ordering below - interface is xz plane, so dx and dz go first
		optical_interface_z(quartzIndex,liquidIndex,dx,dz,dy);
		passed_interface = true;
	}
	
	//Now we just finish
	
	dt = (upperWedgeTop + barLength/2 - y)/dy;
	if (!(x_wedge_coerce_check(x,y,z,dx,dy,dz,dt)))
	{
		//I'm not sure if it can go out the window at this point, but just in case
		//Goes out the window, fail and return
		z = 1337;
		return;
	}
// 	printf("%12.04f %12.04f %12.04f :: %12.04f %12.04f %12.04f dt: %12.04f\n",x,y,z,dx,dy,dz,dt);

	//and we're done.  Everything should be correct now with the ray at the top of the wedge
}
//Possibly inline these or something for speed, but right now, leave them for sanity
bool DircGopticalSim::optical_interface_z(\
	double n1,\
	double n2,\
	double &dx,\
	double &dy,\
	double &dz)
{
	//n1 is starting index of refraction, n2 is ending
	//Assume that it's going through the z plane
	//(this will change the ordering when called in our current geometry)
  
	double before_ang = acos(dz);
	double after_ang = 0;
	
	
	double n12rat = n1/n2;
	double n12rat2 = n12rat*n12rat;
	
	//takes care of trig
	dz = sqrt(1-n12rat2*(1-dz*dz));
	
	if (dz != dz) return false;//total internal reflection
	
	//simpler expression than I expected
	dx *= n12rat;
	dy *= n12rat;
	
// 	printf("in optical interface store_refraction: %s\n",store_refraction ? "true" : "false");
	if (store_refraction == true)
	{
// 	  printf("storing refraction\n");
	  after_ang = acos(dz);
	  refraction_before.push_back(before_ang);
	  refraction_after.push_back(after_ang);
	}
	
	return true;
}
bool DircGopticalSim::x_wedge_coerce_check(\
	double &x,\
	double &y,\
	double &z,\
	double &dx,\
	double &dy,\
	double &dz,\
	double dt)
{
	//already have dt, be sure to add it in when calling this
  
	//assumes x starts between the sides of the wedge
	//assumes dy/dx constant over the whole distance
	//also performs the index plane crossing
	
	//changes x and dx, correctly propagating them
	//propagates y and z as well
	//returns false if it goes out the window on the x step
	//Will not check for y and z directions (assumes dt corresponds to an intersection/reflection point)
	
// 	double cx = x;
// 	double cdx = dx;
	double cdt = 0;
	double tdt = 0;
	double ty = y - barLength/2;
	
	
	while (true)
	{
		//Get time of this bounce based on x direction
		if (dx < 0)
		{
			tdt = x + barWidth/2;
			tdt /= -dx;
		}
		else
		{
			tdt = -x + barWidth/2 - wedgeWidthOff;
			tdt /= dx;
		}
// 		printf("tdt: %12.04f dx: %12.04f\n",tdt,dx);
		if (tdt + cdt < dt)
		{
			//bounced
			//add to total time taken
			cdt += tdt;
			ty += dy*tdt;
			if (ty < upperWedgeBottom && ty > wedgeHeight)
			{
				//out the window on an edge
				return false;
			}
			
			//not out the window, propagate x for the next iteration
			x += tdt*dx;
			if (ty < upperWedgeBottom)
			{
				//bounce off lower wedge
				dx = -dx;
			}
			else
			{
				//if in upper wedge, don't bounce and propagate like normal
				tdt = dt - cdt;
				break;
			}
		}
		else
		{
			//does not bounce - record remaining dt and break;
			tdt = dt - cdt;
			break;
		}
// 		printf("x: %12.04f\n",x);
	}
	//Finish last leg of trip:
	x += dx*tdt;
	y = barLength/2 + ty + dy*tdt;//Note that y starts at zero (bar middle, but ty starts at bar top
	z += dt*dz;//Implemented for convenience - z direction irrelevant to window calculation
	
	return true;//not out the window for the whole trip
}
bool DircGopticalSim::absorbtion_mc(double dx, double dy)
{
	//True if the particle makes it
	//expects vector pointing after bounce
	//Magic number corresponding to minimal distancel traveled
	double min_dist = 650+56;
	
	//approximating here - assume yz distance is independent of where on mirror it hit
	//measuring side plots with a ruler seems to mean this is good to ~5-10%
	double approx_dist = min_dist*sqrt(dx*dx+dy*dy)/dy;
	
	if (store_traveled == true)
	{
		dist_traveled.push_back(approx_dist);
	}
	
	double prob_transmitted = exp(-liquidAbsorbtion*approx_dist);
	if (rand_gen->Uniform(0,1) < prob_transmitted)
	{
		return true;
	}
	else
	{
		return false;
	}
}
double DircGopticalSim::warp_box(\
	double &x,\
	double &y,\
	double &z,\
	double &dx,\
	double &dy,\
	double &dz)
{
  //propagates light from the top of the wedge until just after bouncing off of the focusing mirrors
  //Consider branch profiling or some such other optimizations
  //Also, Fast Math, but don't tell the people in the penthouse
 
  //does not currently include x bouncing off the side - need to add that later
  //possibly in the last "warp to plane" bit
	double rval = 0; //distance traveled

//first reflect off the back of the bar
	double dt;
	
	if (dz > 0)
	{//Condition may not be needed - try without for speed later
		dt = -z/dz;
		
		if (y + dy*dt < focMirrorBottom)
		{
			//reflects off of back
			x += dx*dt;
			y += dy*dt;
			z += dz*dt;
			
			//Abuse the fact that this vector is normalized and hope it stays that way
			dz = -dz;
			rval += dt;
		}
	}
	
	
	double trval  = 0;
	if (three_seg_mirror == true)
	{
// 		printf("preplane x: %8.04f   y: %8.04f   z: %8.04f\n",x,y,z);
		trval = three_seg_reflect(\
		  x,\
		  y,\
		  z,\
		  dx,\
		  dy,\
		  dz);
	}
	else
	{
		 trval = cylindrical_reflect(\
		  x,\
		  y,\
		  z,\
		  dx,\
		  dy,\
		  dz);
	}
	
	if (trval < 0)
	{
		  z = 1337;
		  return -1;
	}

	//short propagate in front of the mirrors:
	//TODO remove later for speed
	x += .1*dx;
	y += .1*dy;
	z += .1*dz;
	
// 	printf("x: %8.04f   y: %8.04f   z: %8.04f   dx: %8.04f    dy: %8.04f    dz: %8.04f\n",x,y,z,dx,dy,dz);
	
	
	rval += trval;
	
	return rval;
}
void DircGopticalSim::plane_reflect(\
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
	double &dt)
{
	//Propagate to the intercept and reflects off a plane
	//Could use this in the Wedge, but would loose some speed

	double n_dot_v = (Nx*dx + Ny*dy + Nz*dz);
	double n_dot_v0 = (Nx*x + Ny*y + Nz*z);
  
	dt = (D - n_dot_v0)/n_dot_v;
	
	x += dt*dx;
	y += dt*dy;
	z += dt*dz;
	
	dx -= 2*n_dot_v*Nx;
	dy -= 2*n_dot_v*Ny;
	dz -= 2*n_dot_v*Nz;
}
double DircGopticalSim::get_z_intercept(\
	double Nx,\
	double Ny,\
	double Nz,\
	double D,\
	double x,\
	double y,\
	double z,\
	double dx,\
	double dy,\
	double dz)
{
	//finds z intercept (without modifying x y and z)
	//Could use this in the Wedge, but would loose some speed
	
	double n_dot_v = (Nx*dx + Ny*dy + Nz*dz);
	double n_dot_v0 = (Nx*x + Ny*y + Nz*z);
	double dt;//optimized by math compression?
	dt = (D - n_dot_v0)/n_dot_v;
	
// 	printf("x: %8.04f   y: %8.04f   z: %8.04f   dx: %8.04f    dy: %8.04f    dz: %8.04f   n_dot_v0: %8.04f tz: %8.04f\n",x,y,z,dx,dy,dz,n_dot_v0,z+dz*dt);
	
	return z + dt*dz;
}
void DircGopticalSim::get_intercept_plane(\
	double Nx,\
	double Ny,\
	double Nz,\
	double D,\
	double &x,\
	double &y,\
	double &z,\
	double dx,\
	double dy,\
	double dz)
{
	//Just intersects a plane (modifying x,y,z)
	//Could use this in the Wedge, but would loose some speed
	
	double n_dot_v = (Nx*dx + Ny*dy + Nz*dz);
	double n_dot_v0 = (Nx*x + Ny*y + Nz*z);
	double dt;//optimized by math compression?
	dt = (D - n_dot_v0)/n_dot_v;
	
// 	printf("x: %8.04f   y: %8.04f   z: %8.04f   dx: %8.04f    dy: %8.04f    dz: %8.04f   n_dot_v0: %8.04f tz: %8.04f\n",x,y,z,dx,dy,dz,n_dot_v0,z+dz*dt);
	
	x += dx*dt;
	y += dy*dt;
	z += dz*dt;
}
double DircGopticalSim::three_seg_reflect(\
		  double &x,\
		  double &y,\
		  double &z,\
		  double &dx,\
		  double &dy,\
		  double &dz)
{
	//Intersects and reflects the three segment plane
	double rval = 0;
	//definitely losing time here - combuting the n_dot_v and n_dot_v0 and dt twice for the chosen plane
	//TODO fix this once the code is correct
	
	//I hope there's a fast way to do these reflections
	
	double tz = 0;
	
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
// 	printf("tz1: %12.04f\n",tz);
	if (tz > threeSeg2Z && tz < 0)
	{//reflect off mirror closest to box
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
			rval);
		return rval;
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
	if (tz > threeSeg3Z && tz < 0)
	{//reflect off middle mirror
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
			rval);
		return rval;
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
	if (tz > -focMirrorZDim && tz < 0)
	{//reflect off mirror closest to box
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
			rval);
		return rval;
	}
	//reflects off of nothing :(
	z = 1337;
	return -1;
}
double DircGopticalSim::sgn(double val) 
{//stolen from http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
    return (0 < val) - (val < 0);
}
double DircGopticalSim::cylindrical_reflect(\
		  double &x,\
		  double &y,\
		  double &z,\
		  double &dx,\
		  double &dy,\
		  double &dz)
{
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
	
	
	
	double norm_loc = sqrt(localNy*localNy + localNz*localNz);//there's gotta be a better way to normalize this
	
// 	printf("new R: %12.04f dy_norm: %12.04f dz_norm: %12.04f\n",norm_loc,dy_norm,dz_norm);
	
	localNy /= norm_loc;
	localNz /= norm_loc;
	
// 	printf("y: %12.04f z: %12.04f lny: %12.04f lnz: %12.04f\n",y,z,localNy,localNz);
	
	//remove for speed in ideal case
// 	rotate_2d(localNx,localNz,cos(foc_yrot/57.3),sin(foc_yrot/57.3));
	
	double n_dot_v = -(dx*localNx + dy*localNy + dz*localNz);
	
// 	printf("dx: %8.04f dy: %8.04f dz: %8.04f\n",dx,dy,dz);
	
	dx += 2*n_dot_v*localNx;
	dy += 2*n_dot_v*localNy;
	dz += 2*n_dot_v*localNz;
	
// 	printf("dx: %8.04f dy: %8.04f dz: %8.04f\n",dx,dy,dz);
	
	return rval;
}