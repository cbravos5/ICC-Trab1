#include <stdio.h>
#include <math.h>
#include "fronteira.h"

double u_fronteira_sul(double x){
	return(sin((2*M_PI*M_PI) - (2*M_PI*x))*sinh(M_PI*M_PI));

}

double u_fronteira_norte(double x){
	return(sin(2*M_PI*x)*sinh(M_PI*M_PI));
}


double fX(double x, double y)
{
	return 4*M_PI*(sin(2*M_PI*x)*sinh(M_PI*y) + sin(2*M_PI*(M_PI-x))*sinh(M_PI*(M_PI-y)));
}
