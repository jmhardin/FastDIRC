
#ifndef DIRC_POINT
#define DIRC_POINT
struct dirc_point
{
	double x;
	double y;
	double t;
	int updown;
	int last_wall_x;
	int wedge_before_interface;
	double weight;
};
#endif
