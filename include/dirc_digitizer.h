#include "dirc_point.h"
#include <vector>

//Root is infecting more files :(
#include <TRandom3.h>

#ifndef DIRC_DIGITIZER
#define DIRC_DIGITIZER
class DircDigitizer
{
private:
	double minx,maxx,miny,maxy;
	double resx, resy;
	double t_unc;
	
	TRandom3 *dig_rand;
	
public:
	DircDigitizer(\
		double iminx,\
		double imaxx,\
		double iresx,\
		double iminy,\
		double imaxy,\
		double iresy,\
		double it_unc);

	void digitize_point(dirc_point &pt);
	void digitize_points(std::vector<dirc_point> &points);
};
#endif
