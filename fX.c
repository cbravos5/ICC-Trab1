#include math.h

double fX(double x, double y)
{
	return 4*M_PI*(sin(2*M_PI*x)*sinh(M_PI*y) + sin(2*M_PI*(M_PI-x))*sinh(M_PI*(M_PI-y)));
}