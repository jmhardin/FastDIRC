/*
 * 
 *      This file is part of the Goptical library.
 *  
 *      The Goptical library is free software; you can redistribute it
 *      and/or modify it under the terms of the GNU General Public
 *      License as published by the Free Software Foundation; either
 *      version 3 of the License, or (at your option) any later version.
 *  
 *      The Goptical library is distributed in the hope that it will be
 *      useful, but WITHOUT ANY WARRANTY; without even the implied
 *      warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *      See the GNU General Public License for more details.
 *  
 *      You should have received a copy of the GNU General Public
 *      License along with the Goptical library; if not, write to the
 *      Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 *      Boston, MA 02111-1307 USA
 *  
 *      Copyright (C) 2010-2011 Free Software Foundation, Inc
 *      Author: Alexandre Becoulet
 * 
 */

/* -*- indent-tabs-mode: nil -*- */

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <Goptical/Analysis/Spot>


#include <Goptical/Math/Vector>
#include <Goptical/Math/VectorPair>

#include <Goptical/Material/Base>
#include <Goptical/Material/Sellmeier>
#include <Goptical/Material/Metal>
#include <Goptical/Material/Mirror>
#include <Goptical/Material/Vacuum>
#include <Goptical/Material/Dielectric>



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


#include <Goptical/Curve/Sphere>
#include <Goptical/Curve/Base>
#include <Goptical/Curve/Flat>
#include <Goptical/Shape/Ellipse>
#include <Goptical/Shape/Rectangle>
#include <Goptical/Shape/Disk>

#include <Goptical/Shape/Disk>
#include <Goptical/Shape/EllipticalRing>

#include <Goptical/Trace/Tracer>
#include <Goptical/Trace/Result>
#include <Goptical/Trace/Distribution>
#include <Goptical/Trace/Sequence>
#include <Goptical/Trace/Params>
#include <Goptical/Trace/Ray>
#include "Goptical/common.hh"

#include <Goptical/Io/RendererSvg>
#include <Goptical/Io/Rgb>
#include <Goptical/Io/RendererViewport>

#include "/usr/include/gsl/gsl_multimin.h"
#include <stdlib.h>
#include <math.h>

#include <TRandom3.h>

using namespace Goptical;
class qwartzCrystal : public Material::Base
{
public:
	Io::Rgb get_color() const
	{
		return Io::Rgb(0,.8,.8);
	}
	double get_internal_transmittance(double wavelength,double thickness) const
	{
		return 1;//perfect for now
	}
	double get_refractive_index(double wavelength) const
	{
		//return 1.473;
		//Edge Kludge for now
		//	return 1.3579;
	}
	bool is_opaque() const
	{
		return false;
	}
	bool is_reflecting() const
	{
		return false;
	}
};
class mineralOil : public Material::Base
{
public:
	Io::Rgb get_color() const
	{
		return Io::Rgb(.8,.8,0);
	}
	double get_internal_transmittance(double wavelength,double thickness) const
	{
		return 1;//perfect for now
	}
	double get_refractive_index(double wavelength) const
	{
		//return 1.470;
		//Edge Kludge for now
		//	return 1.3579;
	}
	bool is_opaque() const
	{
		return false;
	}
	bool is_reflecting() const
	{
		return false;
	}
};
class offSphere : public Curve::Base
{
private:
	double r,x0,y0;
public:
	offSphere(double ir, double ix0, double iy0)
	{
		r = ir;
		x0 = ix0;
		y0 = iy0;
	}
	void derivative(const Math::Vector2 &xy, Math::Vector2 &dxdy) const
	{
		Curve::Sphere onSphere(r);
		Math::Vector<2> offset(x0,y0);
		(dynamic_cast<Curve::Base*>(&onSphere))->derivative(xy-offset,dxdy);
	}
	bool intersect(Math::Vector3 &point, const Math::VectorPair3 &ray) const
	{
		Curve::Sphere onSphere(r);
		Math::Vector<2> offset(x0,y0);
		Math::Vector<3>  offset3(x0,y0,0);
		bool rval;
		Math::Vector3 tpoint;
		rval = onSphere.intersect(*(&tpoint),Math::VectorPair<3>(ray[0]-offset3,ray[1]));
		point = tpoint+offset3;
		return rval;
	}
	void normal(Math::Vector3 &normal, const Math::Vector3 &point) const
	{
		Curve::Sphere onSphere(r);
		Math::Vector<2> offset(x0,y0);
		Math::Vector<3>  offset3(x0,y0,0);
		onSphere.normal(normal,point-offset3);
	} 
	double sagitta(const Math::Vector2 &xy) const
	{
		Curve::Sphere onSphere(r);
		Math::Vector<2> offset(x0,y0);
		return (dynamic_cast<Curve::Base*>(&onSphere))->sagitta(xy-offset);
	}
};

class cylinderY : public Curve::Base
{//Slow implementation - faster with symbolic detrivatives and such
private:
	double r,r2;
	double s;
public:
	cylinderY(double ir)
	{
		r = ir;
		r2 = r*r;
		s=1;
		if (r < 0) s = -1;
	}
	double sagitta(const Math::Vector2 &xy) const
	{
		double y = xy.y();
		double y2 = y*y;
		if (y2>r2){
			return 0;
		}
		else{
			return s*sqrt(r2-y2)-r;
		}
	}
	void derivative(const Math::Vector2 &xy, Math::Vector2 &dxdy) const
	{/*
	double y = xy.y();
	double dy = 0;
	double y2 = y*y;
	if (y2<r2){
		dy= s*y/sqrt(r2-y2);
	}
	dxdy[0]=0;
	dxdy[1]=dy;
	*/
	Base::derivative(xy,dxdy);
	}
	void normal(Math::Vector3 &normal, const Math::Vector3 &point) const
	{
		Base::normal(normal,point);
	}
	bool intersect(Math::Vector3 &point, const Math::VectorPair3 &ray) const
	{
		return Base::intersect(point,ray);
	}
	
};
//Shamelessly copied this implementation of the complex hull form wikibooks:
//http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
//TODO find citation 4/27/13
// Implementation of Andrew's monotone chain 2D convex hull algorithm.
// Asymptotic complexity: O(n log n).
// Practical performance: 0.5-1.0 seconds for n=1000000 on a 1GHz machine.

typedef double coord_t;         // coordinate type
typedef double coord2_t;  // must be big enough to hold 2*max(|coordinate|)^2

struct Point {
	coord_t x, y;
	
	bool operator <(const Point &p) const {
		return x < p.x || (x == p.x && y < p.y);
	}
};

// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
coord2_t cross(const Point &O, const Point &A, const Point &B)
{
	return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last point in the returned list is the same as the first one.
std::vector<Point> convex_hull(std::vector<Point> P)
{
	int n = P.size(), k = 0;
	std::vector<Point> H(2*n);
	
	// Sort points lexicographically
	std::sort(P.begin(), P.end());
	
	// Build lower hull
	for (int i = 0; i < n; i++) {
		while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}
	
	// Build upper hull
	for (int i = n-2, t = k+1; i >= 0; i--) {
		while (k >= t && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}
	
	H.resize(k);
	return H;
}
//End Shameless copying
double hullArea(std::vector<Point> H)
{
	//Integrat ydx for internal area
	double integral = 0;
	int n = H.size();
	double dx;
	double avgy;
	Point pi, pip1;
	for(int i = 0; i < n-1; i++)
	{
		pi = H[i];
		//    pip1 = H[(i+1)%n];//takes car of final segment;
		pip1 = H[(i+1)];//takes car of final segment;
		dx = pip1.x-pi.x;
		avgy = (pip1.y+pi.y)/2;
		integral+=dx*avgy;
	}
	return -integral;
}

double convexHullArea(_Goptical::Trace::rays_queue_t rays)
{
	std::vector<Point> points (rays.size());
	for (unsigned int i = 0; i < rays.size();i++)
	{
		points[i].x = ((Math::Vector3)rays[i]->get_intercept_point()).x();
		points[i].y = ((Math::Vector3)rays[i]->get_intercept_point()).y();
	} 
	std::vector<Point> hull = convex_hull(points); 
	return hullArea(hull);
}	
double rectangleArea(_Goptical::Trace::rays_queue_t rays, Math::Vector<3> centroid,double h, double w, double iweight=100)
{
	bool flatWidth=false;
	double cx = centroid.x();
	double cy = centroid.y();
	double x,y,dx,dy;
	double func = 1;
	double maxx,maxy,minx,miny;
	maxx=0;
	minx=0;
	maxy=0;
	miny=0;
	for (unsigned int i = 0; i < rays.size();i++)
	{
		x = ((Math::Vector3)rays[i]->get_intercept_point()).x();
		y = ((Math::Vector3)rays[i]->get_intercept_point()).y();
		dx = fabs(cx-x);
		dy = fabs(cy-y);
		
		if (dx > h/2)
		{
			func += (dx-h/2)*(dx-h/2)/100;
		}
		if (dy > w/2)
		{
			func += (dy-w/2)*(dy-w/2)/100;
		}
		
		
		maxx=std::max(maxx,x);
		minx=std::min(minx,x);
		maxy=std::max(maxy,y);
		miny=std::min(miny,y);
	}
	double tfunc = 0;
	dx = (maxx-minx);
	dy = (maxy-miny);
	if (dx > h)
	{
		tfunc += (dx-h)*10;
	}
	if (dy > w)
	{
		tfunc += (dy-w)*10;
	}
	//  std::cout << dx << " " << dy << " "  <<tfunc << std::endl;
	
	if (flatWidth==true)
	{
		func = tfunc;
	}
	
	return func*iweight;
}	
//int nvars = 12;
//int nobjs = 1;
double evaluate(const gsl_vector* vars, void* params)
{
	TRandom3 *rand_gen = new TRandom3();
	//**********************************************************************
	// Optical system definition
	//  std::cout << "START EVAL" << std::endl;
	
	//TODO:make these params
	double size = 2200;//overal System scale 13.1 deg = 10+3.1 at 4m
	double testSize = 1800;
	//  double sourceDist = ((double*)params)[1];
	
	//  double barLength=2450;
	double barLength=1550;
	//  double barLength=250;
	double barWidth=50;
	double barDepth=17;
	double sourceDist = barDepth*.99;;
	
	double focR = gsl_vector_get(vars,0);
	double focMirrorSize = gsl_vector_get(vars,1);
	double focRot = gsl_vector_get(vars,2)/100;
	double focY = gsl_vector_get(vars,3);
	double focZ = gsl_vector_get(vars,4);
	double sensSize = gsl_vector_get(vars,5);
	double sensRot = gsl_vector_get(vars,6)/100;
	double sensY = gsl_vector_get(vars,7);
	double sensZ = gsl_vector_get(vars,8);
	
	
	double testy = barLength/2+1000;//location of image in y dircetion
	
	double sourceAngle = 0;
	double sDepth = .95*barDepth;
	double sHeight = sDepth*tan(sourceAngle/57.3);
	double emitAngle = 47;
	double particleTheta = 0;
	double particlePhi = 0;
	
	
	Sys::System             sys;
	
	
	ref<qwartzCrystal> qwartz = ref<qwartzCrystal>::create();
	ref<mineralOil> oil = ref<mineralOil>::create();
	
	int numPhots = 200;
	
	double sourcez = -sDepth;
	double sourcey = -barDepth*tan(particleTheta/57.3);
	double sourcex = 0;
	double tsy = sourcey;
	double tsx = 0;
	sourcey = tsy*cos(particlePhi/57.3)-tsx*sin(particlePhi/57.3);
	sourcex = tsy*sin(particlePhi/57.3)+tsx*cos(particlePhi/57.3);
	/*
	 *  std::mt19937_64 mtGenerator;
	 *  std::uniform_real_distribution<double> uniPhi(0,2*3.14159265);
	 *  std::uniform_real_distribution<double> uniDepth(0,barDepth);
	 */
	Sys::SourceRays srcrays(Math::Vector<3>(sourcex,sourcey,sourcez));
	srcrays.set_material(qwartz);
	
	double sourceOff,randPhi;
	
	for (int i = 0; i < numPhots; i++)
	{ 
		sourceOff = rand_gen->Uniform(0,2*3.14159265);
		randPhi = rand_gen->Uniform(0,barDepth);
		std::cout << sourceOff << std::endl;
		srcrays.add_rays(Math::VectorPair<3>(0,0,sourceOff,sin(emitAngle/57.3)*cos(randPhi),sin(emitAngle/57.3)*sin(randPhi),cos(emitAngle/57.3)),&srcrays);
	}
	srcrays.rotate(particleTheta,0,0);
	srcrays.rotate(0,0,particlePhi);
	
// 	sys.add(srcrays);
	
	
	
	ref<Curve::Flat> barCurve = ref<Curve::Flat>::create();
	ref<Shape::Rectangle> barXYShape = ref<Shape::Rectangle>::create(barWidth,barLength);
	ref<Shape::Rectangle> barYZShape = ref<Shape::Rectangle>::create(barDepth,barLength);
	ref<Shape::Rectangle> barXZShape = ref<Shape::Rectangle>::create(barWidth,barDepth);
	
	
	Sys::OpticalSurface qwartzBarFarPlane(Math::VectorPair<3>(0,0,0,0,0,0),barCurve,barXYShape,qwartz,_Goptical::Material::vacuum);
	sys.add(qwartzBarFarPlane);
	
	Sys::OpticalSurface qwartzBarClosePlane(Math::VectorPair<3>(0,0,-barDepth,0,0,0),barCurve,barXYShape,qwartz,_Goptical::Material::vacuum);
	qwartzBarClosePlane.rotate(0,180,0);
	sys.add(qwartzBarClosePlane);
	
	Sys::OpticalSurface qwartzBarYZ1Plane(Math::VectorPair<3>(-barWidth/2,0,-barDepth/2,0,0,0),barCurve,barYZShape,qwartz,_Goptical::Material::vacuum);
	qwartzBarYZ1Plane.rotate(0,90,0);
	sys.add(qwartzBarYZ1Plane);
	
	Sys::OpticalSurface qwartzBarYZ2Plane(Math::VectorPair<3>(barWidth/2,0,-barDepth/2,0,0,0),barCurve,barYZShape,qwartz,_Goptical::Material::vacuum);
	qwartzBarYZ2Plane.rotate(0,-90,0);
	sys.add(qwartzBarYZ2Plane);
	
	Sys::OpticalSurface qwartzBarXZ1Plane(Math::VectorPair<3>(0,barLength/2+0,-barDepth/2,0,0,0),barCurve,barXZShape,qwartz,oil);
	qwartzBarXZ1Plane.rotate(90,0,0);
	//  sys.add(qwartzBarXZ1Plane);
	
	Sys::OpticalSurface qwartzBarXZ2Plane(Math::VectorPair<3>(0,-barLength/2,-barDepth/2,0,0,0),barCurve,barXZShape,qwartz,_Goptical::Material::vacuum);
	qwartzBarXZ2Plane.rotate(-90,0,0);
	sys.add(qwartzBarXZ2Plane);
	
	
	//Wedge
	
	double wedgeWidthOff = 5;
	double wedgeDepthOff = 10;
	double wedgeFarAngle = .006*57.3;
	double wedgeCloseAngle = 30;
	double wedgeWidth=barWidth - wedgeWidthOff;
	double wedgeDepthLow = barDepth+wedgeDepthOff;
	double wedgeDepthHigh = 79;
	double wedgeHeight = 91;
	
	ref<Curve::Flat> wedgeCurve = ref<Curve::Flat>::create();
	ref<Shape::Rectangle> wedgeFarShape = ref<Shape::Rectangle>::create(wedgeWidth,wedgeHeight/cos(wedgeFarAngle/57.3));
	ref<Shape::Rectangle> wedgeCloseShape = ref<Shape::Rectangle>::create(wedgeWidth,wedgeHeight/cos(wedgeCloseAngle/57.3));
	ref<Shape::Rectangle> wedgeSideShape = ref<Shape::Rectangle>::create(wedgeDepthHigh,wedgeHeight);
	// ref<Shape::Rectangle> wedgeBottomShape = ref<Shape::Rectangle>::create(wedgeDepthLow,wedgeWidth);
	// ref<Shape::Rectangle> wedgeTopShape = ref<Shape::Rectangle>::create(wedgeWidth,wedgeDepthHigh);
	ref<Shape::Rectangle> wedgeSmallShape = ref<Shape::Rectangle>::create(wedgeWidthOff,barDepth);
	
	Sys::OpticalSurface wedgeFarPlane(Math::VectorPair<3>((-barWidth+wedgeWidth)/2,barLength/2+wedgeHeight/2,wedgeHeight/(2/tan(wedgeFarAngle/57.3))),wedgeCurve,wedgeFarShape,qwartz,_Goptical::Material::vacuum);
	wedgeFarPlane.rotate(wedgeFarAngle,0,0);
	sys.add(wedgeFarPlane);
	
	
	Sys::OpticalSurface wedgeClosePlane(Math::VectorPair<3>((-barWidth+wedgeWidth)/2,barLength/2+wedgeHeight/2,-barDepth-wedgeDepthOff-wedgeHeight/(2/tan(wedgeCloseAngle/57.3))),wedgeCurve,wedgeCloseShape,qwartz,_Goptical::Material::vacuum);
	wedgeClosePlane.rotate(0,180,0);
	wedgeClosePlane.rotate(wedgeCloseAngle,0,0);
	sys.add(wedgeClosePlane);
	
	Sys::OpticalSurface wedgeLeftPlane(Math::VectorPair<3>(-barWidth/2,barLength/2+wedgeHeight/2,-wedgeDepthHigh/2),wedgeCurve,wedgeSideShape,qwartz,_Goptical::Material::vacuum);
	wedgeLeftPlane.rotate(0,90,0);
	sys.add(wedgeLeftPlane);
	
	Sys::OpticalSurface wedgeRightPlane(Math::VectorPair<3>(barWidth/2-wedgeWidthOff,barLength/2+wedgeHeight/2,-wedgeDepthHigh/2),wedgeCurve,wedgeSideShape,qwartz,_Goptical::Material::vacuum);
	wedgeRightPlane.rotate(0,-90,0);
	sys.add(wedgeRightPlane);
	
	Sys::OpticalSurface wedgeSmallPlane(Math::VectorPair<3>(barWidth/2-wedgeWidthOff/2,barLength/2,0),wedgeCurve,wedgeSmallShape,qwartz,_Goptical::Material::vacuum);
	wedgeSmallPlane.rotate(90,0,0);
	sys.add(wedgeSmallPlane);
	
	/*  
	 *  Sys::OpticalSurface wedgeTopPlane(Math::VectorPair<3>((-barWidth+wedgeWidth)/2,barLength/2+wedgeHeight,-wedgeDepthHigh/2),wedgeCurve,wedgeTopShape,qwartz,oil);
	 *  wedgeTopPlane.rotate(90,0,0);
	 *  sys.add(wedgeTopPlane);
	 */
	
	//Window implemented as a gap
	//TODO - put in double x,y,z for readability
	
	double upperWedgeDepthHigh = 130;
	double upperWedgeTop = 178.6;
	double upperWedgeHeight = 78;
	
	ref<Curve::Flat> upperWedgeCurve = ref<Curve::Flat>::create();
	ref<Shape::Rectangle> upperWedgeFarShape = ref<Shape::Rectangle>::create(wedgeWidth,upperWedgeHeight);
	ref<Shape::Rectangle> upperWedgeCloseShape = ref<Shape::Rectangle>::create(wedgeWidth,upperWedgeHeight/cos(wedgeCloseAngle/57.3));
	ref<Shape::Rectangle> upperWedgeSideShape = ref<Shape::Rectangle>::create(upperWedgeDepthHigh,upperWedgeHeight);
	// ref<Shape::Rectangle> wedgeBottomShape = ref<Shape::Rectangle>::create(wedgeDepthLow,wedgeWidth);
	ref<Shape::Rectangle> upperWedgeTopShape = ref<Shape::Rectangle>::create(wedgeWidth,upperWedgeDepthHigh);
	
	double uwx = (-barWidth+wedgeWidth)/2;
	double uwy = barLength/2+upperWedgeTop-upperWedgeHeight/2;
	double uwz = 0;
	
	Sys::OpticalSurface upperWedgeFarPlane(Math::VectorPair<3>(uwx,uwy,uwz),upperWedgeCurve,upperWedgeFarShape,qwartz,_Goptical::Material::vacuum);
	wedgeFarPlane.rotate(0,0,0);
	sys.add(upperWedgeFarPlane);
	
	uwx = (-barWidth+wedgeWidth)/2;
	uwy = uwy;
	uwz = -barDepth-wedgeDepthOff-(upperWedgeTop-upperWedgeHeight/2)*tan(wedgeCloseAngle/57.3);
	
	Sys::OpticalSurface upperWedgeClosePlane(Math::VectorPair<3>(uwx,uwy,uwz),upperWedgeCurve,upperWedgeCloseShape,qwartz,_Goptical::Material::vacuum);
	upperWedgeClosePlane.rotate(0,180,0);
	upperWedgeClosePlane.rotate(wedgeCloseAngle,0,0);
	sys.add(upperWedgeClosePlane);
	
	uwx = barWidth/2-wedgeWidthOff;
	uwy = uwy;
	uwz = -upperWedgeDepthHigh/2;
	
	Sys::OpticalSurface upperWedgeRightPlane(Math::VectorPair<3>(uwx,uwy,uwz),upperWedgeCurve,upperWedgeSideShape,qwartz,_Goptical::Material::vacuum);
	upperWedgeRightPlane.rotate(0,-90,0);
	sys.add(upperWedgeRightPlane);
	
	uwx = -barWidth/2;
	uwy = uwy;
	uwz = -upperWedgeDepthHigh/2;
	
	Sys::OpticalSurface upperWedgeLeftPlane(Math::VectorPair<3>(uwx,uwy,uwz),upperWedgeCurve,upperWedgeSideShape,qwartz,_Goptical::Material::vacuum);
	upperWedgeLeftPlane.rotate(0,90,0);
	sys.add(upperWedgeLeftPlane);
	
	uwx = (-barWidth+wedgeWidth)/2;
	uwy = barLength/2+upperWedgeTop;
	uwz = -upperWedgeDepthHigh/2;
	
	Sys::OpticalSurface upperWedgeTopPlane(Math::VectorPair<3>(uwx,uwy,uwz),upperWedgeCurve,upperWedgeTopShape,qwartz,oil);
	upperWedgeTopPlane.rotate(90,0,0);
	sys.add(upperWedgeTopPlane);
	//Box
	//Focus Mirror
	double boxSize = barWidth*35;
	double mirrorDepth = 288;
	focRot = 74.11;
	focR = -1200;
	focMirrorSize = mirrorDepth;
	
	double focSag = focR-sqrt(focR*focR-focMirrorSize*focMirrorSize/4);
	double focAng = atan(focSag*2/focMirrorSize) - focRot/57.3;
	
	//  double focz = -focSag*cos(focRot/57.3)+focMirrorSize*sin(focRot/57.3)/2;
	//  double focy = focSag*sin(focRot/57.3)+focMirrorSize*cos(focRot/57.3)/2;
	
	double focz = -focMirrorSize*sin(focRot/57.3)/2;
	double focy = focMirrorSize*cos(focRot/57.3)/2;
	
	double focYoff = 139;
	
	double fx,fy,fz;
	fx = 0;
	fy = focy + barLength/2+upperWedgeTop+focYoff;
	fz = focz;
	
	ref<Shape::Rectangle> focMirrorShape = ref<Shape::Rectangle>::create(boxSize,focMirrorSize);
	cylinderY focCylinder(-focR);
	ref<Sys::Mirror> focMirror = ref<Sys::Mirror>::create(Math::Vector<3>(fx,fy,fz),focCylinder,focMirrorShape,true,Material::mirror,oil);
	focMirror->rotate(focRot,0,0);
	sys.add(focMirror);
	
	
	//Image Plane
	
	sensSize = 312;
	//sensSize = 600;
	sensRot = 47.87;
	double boxCloseZ = -614;
	
	double reflOff = 9;
	double sx,sy,sz;
	sx = 0;
	sy = -sensSize*cos(sensRot/57.3)/2-reflOff+barLength/2;
	sz = boxCloseZ + sensSize*sin(sensRot/57.3)/2;
	
	ref<Shape::Rectangle> boxImageShape = ref<Shape::Rectangle>::create(boxSize,sensSize);
	ref<Curve::Flat> boxImageCurve = ref<Curve::Flat>::create();
	Sys::Image boxImage(Math::Vector<3>(sx,sy,sz),boxImageCurve,boxImageShape);
	boxImage.rotate(sensRot,0,0);
	//  Sys::Image boxImage(Math::Vector<3>(0,imageendy/2+barLength/2+wedgeHeight,gsl_vector_get(vars,0)),boxImageCurve,boxImageShape);
	//  boxImage.rotate(0,0,0);
	sys.add(boxImage);
	
	
	ref<Shape::Rectangle> blockShape = ref<Shape::Rectangle>::create(boxSize,upperWedgeTop);
	ref<Sys::Mirror> opaqueMirror=ref<Sys::Mirror>::create(Math::Vector<3>(0,barLength/2+upperWedgeTop,0),barCurve,blockShape,true,Material::mirror,oil);
	sys.add(opaqueMirror);
	
	
	
	ref<Curve::Flat> imgCurve = ref<Curve::Flat>::create();
	ref<Shape::Rectangle> imgShape = ref<Shape::Rectangle>::create(size/0.5,testSize*2);
	Sys::Image testImage(Math::VectorPair<3>(0,testy,0,0,0,0),imgCurve,imgShape);
	testImage.rotate(90,0,0);
	sys.add(testImage);
	
	//**********************************************************************
	// Render optical layout with rays in svg file
	
	printf("before add\n");
	sys.add(srcrays);

	std::cout << "Before Trace" << std::endl;
	Trace::Tracer         tracer(sys);
	
	tracer.get_params().set_max_bounce(1000);
	tracer.get_trace_result().set_generated_save_state(srcrays);
	
	tracer.trace();
	
	/* anchor layout */
	if (((double*) params)[0] < 0)
	{
		
		Io::RendererSvg       svg_renderer("layout.svg", 1000, 700);
		Io::RendererViewport  &renderer = svg_renderer;
		// 3d system layout
		
		std::cout << renderer.get_feature_size() << std::endl;
		renderer.set_feature_size(20);
		
		renderer.set_perspective();
		renderer.set_fov(45);	
		
		
		sys.draw_3d_fit(renderer, 0);
		renderer.set_camera_transform(renderer.get_camera_transform().linear_rotation(Math::Vector3(0,0,0)));
		sys.draw_3d(renderer);
		
		tracer.get_trace_result().draw_3d(renderer);
		Analysis::Spot spot(sys);
		spot.get_tracer().get_params().set_max_bounce(5000);
		Io::RendererSvg renderSpot("spot.svg", 1000, 1000, Io::rgb_black);
		
		spot.draw_diagram(renderSpot);
	} 
	
	/* anchor layout */
	
	//After each convergence, output
	return -1;
	//  }
	
	}
	void evaluatedf(const gsl_vector* vars, void* params, gsl_vector* g)
	{
		double stepsize = .0001;
		double maxgrad = .1;
		double Fi,Ff;
		double dxi;
		double xi;
		double dfdxi;
		gsl_vector* vars_f = gsl_vector_alloc(vars->size);
		Fi = evaluate(vars,params);
		for (unsigned int i = 0; i < vars->size; i++)
		{
			gsl_vector_memcpy(vars_f,vars);
			xi = gsl_vector_get(vars,i);
			dxi = stepsize*xi;
			gsl_vector_set(vars_f,i,xi+dxi);
			Ff = evaluate(vars_f,params);
			dfdxi = (Ff-Fi)/dxi;
			dfdxi = std::min(fabs(xi)*maxgrad,dfdxi);
			dfdxi = std::max(-fabs(xi)*maxgrad,dfdxi);
			if (Ff==Fi) dfdxi = stepsize*xi;
			
			gsl_vector_set(g,i,dfdxi);
		}
	}
	void evaluatefdf(const gsl_vector* vars, void* params, double* f, gsl_vector* g)
	{
		double stepsize = .0001;
		double maxgrad = .1;
		double Fi,Ff;
		double dxi;
		double xi;
		double dfdxi;
		gsl_vector* vars_f = gsl_vector_alloc(vars->size);
		Fi = evaluate(vars,params);
		for (unsigned int i = 0; i < vars->size; i++)
		{
			gsl_vector_memcpy(vars_f,vars);
			xi = gsl_vector_get(vars,i);
			dxi = stepsize*xi;
			gsl_vector_set(vars_f,i,xi+dxi);
			Ff = evaluate(vars_f,params);
			dfdxi = (Ff-Fi)/dxi;
			dfdxi = std::min(fabs(xi)*maxgrad,dfdxi);
			dfdxi = std::max(-fabs(xi)*maxgrad,dfdxi);
			gsl_vector_set(g,i,dfdxi);
		}
		*f = Fi;
	}
	int gsl_min(double* vals, int nvals, double sdist, double tsize, double* outval) {
		//Modified GSL example
		int nvars = nvals;
		int m_iter = 2;
		int disallow_max = 4000;
		int disallow_count = 0;
		double parnp[5] = {1.0,sdist,tsize,0,0};
		gsl_vector *ss, *x;//ss not always used, but that's ok
		
		const gsl_multimin_fminimizer_type *T = 
		//         gsl_multimin_fminimizer_nmsimplex2rand;
		gsl_multimin_fminimizer_nmsimplex2;
		gsl_multimin_fminimizer *s = NULL;
		gsl_multimin_function minex_func;
		
		/*
		 *       const gsl_multimin_fdfminimizer_type *T;
		 *       gsl_multimin_fdfminimizer *s;      
		 *       T = gsl_multimin_fdfminimizer_conjugate_fr;
		 *       s = gsl_multimin_fdfminimizer_alloc (T, nvars);
		 * 
		 *       gsl_multimin_function_fdf minex_func;
		 *       minex_func.df = evaluatedf;
		 *       minex_func.fdf = evaluatefdf;
		 */
		
		/* Initialize method and iterate */
		minex_func.n = nvars;
		minex_func.f = evaluate;
		minex_func.params = parnp;
		
		
		s = gsl_multimin_fminimizer_alloc (T, nvars);
		
		size_t iter = 0;
		int status;
		double size;
		
		/* Starting point */
		x = gsl_vector_alloc (nvars);
		
		for(int i = 0; i < nvals; i++)
		{
			gsl_vector_set(x,i,vals[i]);
		}     
		
		/* Set initial step sizes to 1 */
		ss = gsl_vector_alloc (nvars);
		gsl_vector_set_all (ss, 550);
		
		gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
		
		do
		{
			iter++;
			status = gsl_multimin_fminimizer_iterate(s);
			
			if (status) 
				break;
			
			size = gsl_multimin_fminimizer_size (s);
			status = gsl_multimin_test_size (size, 1e-2);
			
			if (status == GSL_SUCCESS)
			{
				printf ("converged to minimum at\n");
			}
			
			printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f sdist = %5.1f\n", 
					(int) iter,
					gsl_vector_get (s->x, 0), 
					gsl_vector_get (s->x, 1), 
					s->fval, size,sdist);
			*outval=s->fval;
			
			if (*outval > tsize*tsize)
			{
				disallow_count++;
				if (disallow_count >= disallow_max)
				{
					break;
				}
			}
			else
			{
				disallow_count = 0;
			}
		}
		while (status == GSL_CONTINUE && iter < (unsigned int) m_iter);
						
						/*
						 *       gsl_multimin_fdfminimizer_set (s, &minex_func, x, 0.001, 1e-6);
						 *     
						 *       do
						 *         {
						 *           iter++;
						 *           status = gsl_multimin_fdfminimizer_iterate (s);
						 *     
						 *           if (status)
						 *             break;
						 *     
						 *           status = gsl_multimin_test_gradient (s->gradient, 1e-3);
						 *     
						 *           if (status == GSL_SUCCESS)
						 *             printf ("Minimum found at:\n");
						 *     
						 *           printf ("%5d %.5f %.5f %10.5f\n", iter,
						 *                   gsl_vector_get (s->x, 0), 
						 *                   gsl_vector_get (s->x, 8), 
						 *                   s->f);
						 *     
						 }
						 while (status == GSL_CONTINUE && iter < m_iter);
						 */  
						
						for(int i = 0; i < nvals; i++)
						{
							std::cout << gsl_vector_get(s->x,i) << std::endl;
							vals[i] = gsl_vector_get(s->x,i);
						}
						
						
						gsl_multimin_fminimizer_free (s);
						//       gsl_multimin_fdfminimizer_free (s);
						
						gsl_vector_free(x);
						gsl_vector_free(ss);
						
						return status;
						
						 }
						 int main(int nargs, char* argv[])
						 {
							 int gsl_iter = 10000;
							 gsl_iter = 1;
							 int nvals = 9;
							 double sdist = 100;
							 //	double targetDist = 4000;
							 double targetDist =sdist;
							 double tsize = 3000;
							 double last_stable = sdist;
							 double* vals = new double[nvals];
							 double*outval = new double[1];
							 
							 for(int i = 1; i <= nvals; i++)
							 {
								 vals[i-1] = atof(argv[i]);
								 std::cout << vals[i-1] << std::endl;
							 }
							 
							 srand(20);
							 int status;
							 for(int i = 0; i < gsl_iter; i++)
							 {
								 std::cout << "./newton ";
								 for(int i = 0; i < nvals; i++)
								 {
									 std::cout << vals[i] << " ";
								 }
								 std::cout << std::endl;
								 
								 gsl_vector* x;
								 x = gsl_vector_alloc (nvals);
								 
								 for(int i = 0; i < nvals; i++)
								 {
									 gsl_vector_set(x,i,vals[i]);
								 }     
								 double parp[5] = {-1,sdist,tsize,0,0};
								 std::cout << "eval: " << evaluate(x,parp) << std::endl;
								 
								 status = gsl_min(vals,nvals,sdist,tsize,outval);
								 
								 if (outval[0] < tsize*tsize)
								 {//move in a bit
								 last_stable = std::min(last_stable,sdist);
								 double dsd = (sdist - targetDist)/50;
								 dsd = std::max(dsd,targetDist/200);
								 dsd = std::min(dsd,sdist-targetDist);
								 sdist -= dsd;
								 
							 }
							 else
							 {
								 double dsd = (last_stable-sdist)/2+1;
								 sdist += dsd;
								 for (int j = 0; j < nvals; j++)
								 {
									 //				double perturb = (rand()/(INT_MAX*1.))*.04+.98;
									 //				vals[i]*= perturb;
								 }			
							 }
							 std::cout << "sdist: " << sdist << std::endl;
						 }
						 return status;
	}
	