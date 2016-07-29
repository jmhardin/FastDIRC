#include "dirc_point.h"
#include "dirc_digitizer.h"
#include <vector>

//Root is infecting more files :(
#include <TRandom3.h>

#ifndef DIRC_BABAR_DIGITIZER
#define DIRC_BABAR_DIGITIZER
class DircBaBarDigitizer : public DircDigitizer
{
private:
	double minx,maxx,miny,maxy;
	double pmt_r;
	double pmt_sep_x,pmt_sep_y;
	double second_grid_minx,second_grid_miny;
	double resx, resy;
	double t_unc;
	double t_bin_size;
	
	TRandom3 *dig_rand;
	
public:
	DircBaBarDigitizer(\
		double iminx,\
		double imaxx,\
		double iminy,\
		double imaxy,\
		double ipmt_r,\
		double it_unc,\
		double it_bin_size);

	void digitize_point(dirc_point &pt);
	//return negative time if point is lost
	void digitize_points(std::vector<dirc_point> &points);
};
#endif
