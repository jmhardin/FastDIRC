#include "dirc_point.h"

#ifndef DIRC_LUT_ENUM
#define DIRC_LUT_ENUM
class DircLUTEnum
{
protected:
	int max_n;
	
public:
	DircLUTEnum();

	virtual int return_enum(dirc_point &pt) {return -1;};
	int get_max_enum();
};
#endif
