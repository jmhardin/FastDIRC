#include "dirc_point.h"
#include "../include/dirc_digitizer.h"
#include <vector>
DircDigitizer::DircDigitizer(\
		double iminx,\
		double imaxx,\
		double iresx,\
		double iminy,\
		double imaxy,\
		double iresy)
{
	minx = iminx;
	maxx = imaxx;
	resx = iresx;
	miny = iminy;
	maxy = imaxy;
	resy = iresy;
}
void DircDigitizer::digitize_point(dirc_point &pt)
{
		//overflow/underflow?
	double x = pt.x;
	double y = pt.y;
	int xdig = (x - minx)/resx;
	int ydig = (y - miny)/resy;
	
	double xout = resx*xdig + minx + resx/2;
	double yout = resy*ydig + miny + resy/2;
	
	if (x < minx)
	{
		xout = minx - resx/2;
	}
	else if (x > maxx)
	{
		xout = maxx + resx/2;
	}
	
	if (y < miny)
	{
		yout = miny - resy/2;
	}
	else if (y > maxy)
	{
		yout = maxy + resy/2;
	}
	
	pt.x = xout;
	pt.y = yout;
}
void DircDigitizer::digitize_points(std::vector<dirc_point> &points)
{
	for (unsigned int i = 0; i < points.size(); i++)
	{
		//Wasting time Make it faster later TODO
		digitize_point(points[i]);
	}
}

	