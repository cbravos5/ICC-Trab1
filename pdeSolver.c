#include "gs.h"
#include "sistLinear.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void main(int argc, char *argv[])
{

	double nx, ny, hx, hy;
	int ite = atoi(argv[6]);
	sscanf(argv[2], "%lf", &nx);
	sscanf(argv[4], "%lf", &ny);
	struct matPent A;

	hx = M_PI/nx;
	hy = M_PI/ny;

	A.Dp = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Ds = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Di = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Dsa = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Dia = malloc((nx-2)*(ny-2)*sizeof(double));
	A.X = malloc((nx-2)*(ny-2)*sizeof(double));
	A.B = malloc((nx-2)*(ny-2)*sizeof(double));
	PreencheVetores(nx, ny, hx, hy, &A);
	contorno(nx, ny, hx, hy, &A);
	Gs(A, nx-4 ,(nx-2)*(ny-2), ite);
	//printVet(nx, ny, &A);
}