#include "dirc_point.h"
#include "dirc_lut_enum.h"

#ifndef DIRC_GLUEX_LUT_ENUM
#define DIRC_GLUEX_LUT_ENUM
class DircGluexLUTEnum: public DircLUTEnum
{
private:
	double minx,maxx,miny,maxy,resx,resy;
	int x_n, y_n;
public:
	DircGluexLUTEnum(double iminx,\
                double imaxx,\
                double iresx,\
                double iminy,\
                double imaxy,\
                double iresy);

	int return_enum(dirc_point &pt);
};
#endif
