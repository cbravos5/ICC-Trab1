#include "gs.h"
#include "sistLinear.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void main(int argc, char *argv[])
{

	double nx, ny, hx, hy, *timeI;
	int maxIter = atoi(argv[6]);
	sscanf(argv[2], "%lf", &nx);
	sscanf(argv[4], "%lf", &ny);
	struct matPent A;
	timeI = malloc((nx -2)*(ny -2)*sizeof(double));

	hx = M_PI/(nx - 1);
	hy = M_PI/(ny - 1);

	
	A.T = 0.0;
	A.r = malloc(maxIter*sizeof(double));
	A.Dp = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Ds = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Di = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Dsa = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Dia = malloc((nx-2)*(ny-2)*sizeof(double));
	A.X = malloc((nx-2)*(ny-2)*sizeof(double));
	A.B = malloc((nx-2)*(ny-2)*sizeof(double));
	PreencheVetores(nx, ny, hx, hy, &A);
	contorno(nx, ny, hx, hy, &A);
	Gs(&A, nx-4 ,(nx-2)*(ny-2), maxIter);


	if(argc > 7)
		out(A, nx, ny, hx, hy, maxIter, argv[8], 1);
	else
		out(A, nx, ny, hx, hy, maxIter, NULL, 0);
}