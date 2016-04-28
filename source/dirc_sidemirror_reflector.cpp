#include "../include/dirc_sidemirror_reflector"
#include <vector>
#include "../include/dirc_point.h"

DircSideMirrorReflector::DircSideMirrorReflector(double ixl, double ixr)
{
  xl = ixl;
  xr = ixr;
}
DircSideMirrorReflector::reflect_points(std::vector<dirc_point> &points)
{
  double tmpx = 0;
  double dx;
  for (unsigned int i = 0; i < points.size(); i++)
  {
    tmpx = points[i].x;
    while (tmpx < xl || tmpx > xr)
    {
      if (tmpx < xl)
      {
	tmpx = 2*xl - tmpx;
      }
      if (tmpx > xr)
      {
	tmpx = 2*xr - tmpx;
      }
    }
    points[i].x = tmpx;
  }
}
  