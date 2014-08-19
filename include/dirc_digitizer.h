#include "dirc_point.h"
#include <vector>

#ifndef DIRC_DIGITIZER
#define DIRC_DIGITIZER
class DircDigitizer
{
private:
	double minx,maxx,miny,maxy;
	double resx, resy;
	
public:
	DircDigitizer(\
		double iminx,\
		double imaxx,\
		double iresx,\
		double iminy,\
		double imaxy,\
		double iresy);

	void digitize_point(dirc_point &pt);
	void digitize_points(std::vector<dirc_point> &points);
};
#endif
