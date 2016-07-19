#include "dirc_point.h"
#include <vector>

#ifndef DIRC_DIGITIZER
#define DIRC_DIGITIZER
class DircDigitizer
{
public:
	DircDigitizer() {};

	void digitize_point(dirc_point &pt){};
	void digitize_points(std::vector<dirc_point> &points){};
};
#endif
