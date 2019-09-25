#include <stdio.h>
#include <math.h>

#define PI 3.14159265

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

void contorno(int nx, int ny, struct matPent *A)
{
	int max = (nx-2)*(nx-2);
	int p = -1;
	double hx2, hy2;
	hx2 = hx2*hx2;
	hy2 = hy2*hy2;
	for (int j = 0; j < max; j++)
	{
		for (int i = 0; i < max; ++i)
		{
			p = p+1;
			A->B[p] = 2*hx2*hy2*fX(j*hx,i*hy)
			//ZEROS DA DIAGONAL INFERIOR
			if (i == 0 && j > 1)
				A->Di[p] = 0; 
			//ZEROS DA DIAGONAL SUPERIOR
			if (i == nx-2 && j < ny-2)
				A->Ds[p] = 0;
			//CONTORNO SUL
			if(j == 1)
				A->B[p] += (hy-2)*hx2*u_fronteira_sul(j*hx);
			//CONTORNO NORTE
			if(j == ny-2)
				A->B[p] -= (hy-2)*hx2*u_fronteira_norte(j*hx); 
		}
	}
}
