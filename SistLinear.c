#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "sistLinear.h"
#include "fronteira.h"

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
	// B e X
	for(int i = 0; i < max; i++){
		A->X[i] = 0;
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

		printf("B\n");
	for(int i = 0; i < (nx-2)*(ny-2); i++){
		printf("%f ", A->B[i]);
	}
	printf("\n\n");

		printf("RESULTADO\n");
	for(int i = 0; i < (nx-2)*(ny-2); i++){
		printf("%1.10f ", A->X[i]);
	}
	printf("\n\n");


}

void contorno(int nx, int ny, double hx, double hy, struct matPent *A)
{
	int max = (nx-2)*(nx-2);
	int p = -1;
	double hx2, hy2;
	hx2 = hx*hx;
	hy2 = hy*hy;
	for (int j = 0; j < ny-2; j++)
	{
		for (int i = 0; i < nx-2; ++i)
		{
			p = p+1;
			A->B[p] = 2*hx2*hy2*fX(j*hx,i*hy);
			//ZEROS DA DIAGONAL INFERIOR
			if(i == 0 && j > 0)
				A->Di[p] = 0; 
			//ZEROS DA DIAGONAL SUPERIOR
			if(i == nx-3 && j < ny-3)
				A->Ds[p] = 0;
			//CONTORNO SUL
			if(j == 0)
				A->B[p] += (hy-2)*hx2*u_fronteira_sul(j*hx);
			//CONTORNO NORTE
			if(j == ny-3)
				A->B[p] -= (hy-2)*hx2*u_fronteira_norte(j*hx); 
		}
	}
}







