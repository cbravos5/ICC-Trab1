#include <math.h>
#include <stdlib.h>
#include <stdio.h>

struct matPent
{
	double *Dp, *Ds, *Di, *Dsa, *Dia, *X, *B;	
};

void PreencheVetores(int nx, int ny, double hx, double hy, struct matPent *A)
{
	double hx2, hy2;
	hx2 = hx*hx;
	hy2 = hy*hy;
	int max = (nx-2)*(ny-2);

	//DIAGONAL PRINCIPAL
	double CTE = 2*(hx2 + hy2 + 4*M_PI*M_PI*hx2*hy2);
	for(int i = 0; i < max; i++)
	{
		A->Dp[i] = CTE;
	}
	
	//DIAGONAL SUPERIOR
	CTE = (hx - 2)*hy2;
	for (int i = 0; i < max - 1 ; ++i)
	{
		A->Ds[i] = CTE;
	}
	A->Ds[max-1] = 0;

	//DIAGONAL INFERIOR
	CTE = CTE*(-1);
	for (int i = 1; i < max; ++i)
	{
		A->Di[i] = CTE;	
	}
	A->Di[0] = 0;

	//DIAGONAL SUPERIOR AFASTADA
	CTE = (hy - 2)*hx2;
	for (int i = 0; i < max; ++i)
	{
		if(i < max-nx)
			A->Dsa[i] = CTE;
		else
			A->Dsa[i] = 0;
	}

	//DIAGONAL SUPERIOR AFASTADA
	CTE = CTE*(-1);
	for (int i = 0; i < max; ++i)
	{
		if(i >= nx)
			A->Dia[i] = CTE;
		else
			A->Dia[i] = 0;
	}
}


void printVet(int nx, int ny, struct matPent *A)
{
	printf("DIAGONAL PRINCIPAL\n");
	for (int i = 0; i < (nx-2)*(ny-2); ++i)
	{
		printf("%f  ", A->Dp[i]);
	}
	printf("\n\n");

		printf("DIAGONAL SUPERIOR\n");
	for (int i = 0; i < (nx-2)*(ny-2); ++i)
	{
		printf("%f  ", A->Ds[i]);
	}
	printf("\n\n");

		printf("DIAGONAL INFERIOR\n");
	for (int i = 0; i < (nx-2)*(ny-2); ++i)
	{
		printf("%f  ", A->Di[i]);
	}
	printf("\n\n");

		printf("DIAGONAL SUPERIOR AFASTADA\n");
	for (int i = 0; i < (nx-2)*(ny-2); ++i)
	{
		printf("%f  ", A->Dsa[i]);
	}
	printf("\n\n");

		printf("DIAGONAL INFERIOR AFASTADA\n");
	for (int i = 0; i < (nx-2)*(ny-2); ++i)
	{
		printf("%f  ", A->Dia[i]);
	}
	printf("\n\n");
}

void main(int argc[], char **argv[])
{
	double nx, ny, hx, hy;
	nx = ny = 5.0;
	hx = hy = 1.567;
	struct matPent A;
	A.Dp = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Ds = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Di = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Dsa = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Dia = malloc((nx-2)*(ny-2)*sizeof(double));
	PreencheVetores(nx, ny, hx, hy, &A);
	printVet(nx, ny, &A);
}






