#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define MAX 55

struct matPent
{
	double *Dp, *Ds, *Di, *Dsa, *Dia, *X, *B;	
};

double fX(double x, double y)
{
	return 4*M_PI*(sin(2*M_PI*x)*sinh(M_PI*y) + sin(2*M_PI*(M_PI-x))*sinh(M_PI*(M_PI-y)));
}

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
		double div = RAND_MAX/(1000 - 100);
		A->B[i] = rand()/div;
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

void gs(struct matPent m, int afast, double eps){
	double Xold[MAX];
	double e = eps + 1.0;
	unsigned int j = 0;
	for(int i = 0; i < MAX; i++)
		Xold[i] = 0;
	while(j < eps/*e > eps*/){
		m.X[0] = (m.B[0] - m.Ds[0]*m.X[1] - m.Dsa[0]*m.X[afast])/m.Dp[0];
		for(int i = 1; i < afast; i++){
			m.X[i] = (m.B[i] - m.Di[i]*Xold[i - 1] - m.Ds[i]*m.X[i + 1] - m.Dsa[i]*m.X[i + afast])/m.Dp[i];
		}
		for(int i = afast; i < MAX - afast; i++)
			m.X[i] = (m.B[i] - m.Dia[i]*Xold[i - afast] - m.Di[i]*Xold[i - 1] - m.Ds[i]*m.X[i + 1] - m.Dsa[i]*m.X[i + afast])/m.Dp[i];
		for(int i = MAX - afast; i < MAX - 1; i++)
			m.X[i] = (m.B[i] - m.Dia[i]*Xold[i - afast] - m.Di[i]*Xold[i - 1] - m.Ds[i]*m.X[i + 1])/m.Dp[i];
		m.X[MAX - 1] = (m.B[MAX - 1] - m.Dia[MAX - 1]*Xold[MAX -1 - afast] - m.Di[MAX - 1]*Xold[MAX - 2])/m.Dp[MAX - 1];
		for(int i = 0; i < MAX; i++)
			Xold[i] = m.X[i];

		j++;

		if((j % 2) == 0){
			printf("iteração %d\n", j);
			for(int i = 0; i < MAX; i++)
				printf("%1.15f  ", m.X[i]);
			printf("\n\n");
		}
		/*e = X[0] - Xold[0];
		for(int i = 1; i < MAX; i++){
			if(fabs(X[i] - Xold [i]) > e)
				e = fabs(X[i] - Xold[i]);
		}*/
	}
}

void main(int argc[], char **argv[])
{
	double nx, ny, hx, hy;
	nx = 13.0;
	ny = 7.0;
	hx = hy = 1.567;
	struct matPent A;
	A.Dp = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Ds = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Di = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Dsa = malloc((nx-2)*(ny-2)*sizeof(double));
	A.Dia = malloc((nx-2)*(ny-2)*sizeof(double));
	A.X = malloc((nx-2)*(ny-2)*sizeof(double));
	A.B = malloc((nx-2)*(ny-2)*sizeof(double));
	PreencheVetores(nx, ny, hx, hy, &A);
	gs(A,11,10);
	printVet(nx, ny, &A);
}






