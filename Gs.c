#include "gs.h"
#include "sistLinear.h"
#include <math.h>
#include <stdio.h> //RETIRAR DEPOIS

double normaAprox(double *Xold, double *Xnew, int MAX)
{
	double r = 0;

	for (int i = 0; i < MAX; ++i)
	{
		r += (Xnew[i] - Xold[i])*(Xnew[i] - Xold[i]);
	}
	return sqrt(r);
}

double normaL2(double *Xnew, struct matPent *A, int MAX, int afast)
{
	double r = 0;
	//PRIMEIRA ATRIBUICAO
		r += (A->B[0] - A->Dp[0]*Xnew[0] + A->Ds[0]*Xnew[1] + A->Dsa[0]*Xnew[afast + 2])*
			 (A->B[0] - A->Dp[0]*Xnew[0] + A->B[0] + A->Ds[0]*Xnew[1] + A->Dsa[0]*Xnew[afast + 2]); 

		//LOOP EA QUE DIAGONAL INFERIOR AFASTADA AINDA NAO APARECEU
		for(int i = 1; i < afast + 2; i++)
			r += (A->B[i] - A->Dp[i]*Xnew[i] + A->Di[i]*Xnew[i - 1] + A->Ds[i]*Xnew[i + 1] + A->Dsa[i]*Xnew[afast + 2 + i])*
				 (A->B[i] - A->Dp[i]*Xnew[i] + A->Di[i]*Xnew[i - 1] + A->Ds[i]*Xnew[i + 1] + A->Dsa[i]*Xnew[afast + 2 + i]);

		//LOOP EA QUE AABAS DIAGONAIS AFASTADAS ESTAO PRESENTES
		for(int i = afast + 2; i < MAX - afast - 2; i++)
			r += (A->B[i] - A->Dp[i]*Xnew[i] + A->Dia[i]*Xnew[i - afast - 2] + A->Di[i]*Xnew[i - 1] + A->Ds[i]*Xnew[i + 1] + A->Dsa[i]*Xnew[afast + 2 + i])*
				 (A->B[i] - A->Dp[i]*Xnew[i] + A->Dia[i]*Xnew[i - afast - 2] + A->Di[i]*Xnew[i - 1] + A->Ds[i]*Xnew[i + 1] + A->Dsa[i]*Xnew[afast + 2 + i]);
		//LOOP EA QUE DIAGONAL SUPERIROR AFASTADA ACABOU
		for(int i = MAX - afast - 2; i < MAX - 1; i++)
			r += (A->B[i] - A->Dp[i]*Xnew[i] + A->Dia[i]*Xnew[i - afast - 2] + A->Di[i]*Xnew[i - 1] + A->Ds[i]*Xnew[i + 1])*
				 (A->B[i] - A->Dp[i]*Xnew[i] + A->Dia[i]*Xnew[i - afast - 2] + A->Di[i]*Xnew[i - 1] + A->Ds[i]*Xnew[i + 1]);
		
		//ULTIAA ATRIBUICAO
		r += (A->B[MAX - 1] - A->Dp[MAX - 1]*Xnew[MAX - 1] + A->Dia[MAX - 1]*Xnew[MAX - afast - 3] + A->Di[MAX - 1]*Xnew[MAX - 2])*
			 (A->B[MAX - 1] - A->Dp[MAX - 1]*Xnew[MAX - 1] + A->Dia[MAX - 1]*Xnew[MAX - afast - 3] + A->Di[MAX - 1]*Xnew[MAX - 2]);

	return sqrt(r);
}

void Gs(struct matPent m, int afast, int MAX, int ite){
	double Xold[MAX], r;
	unsigned int j = 0;
	
	for(int i = 0; i < MAX; i++)
		Xold[i] = 0;
	while(j < ite){
		
		//PRIMEIRA ATRIBUICAO
		m.X[0] = m.B[0]/m.Dp[0] - m.Ds[0]*Xold[1]/m.Dp[0] - m.Dsa[0]*Xold[afast + 2]/m.Dp[0]; 

		//LOOP EM QUE DIAGONAL INFERIOR AFASTADA AINDA NAO APARECEU
		for(int i = 1; i < afast + 2; i++)
			m.X[i] = m.B[i]/m.Dp[i] - m.Di[i]*m.X[i - 1]/m.Dp[i] - m.Ds[i]*Xold[i + 1/m.Dp[i]] - m.Dsa[i]*Xold[afast + 2 + i]/m.Dp[i];

		//LOOP EM QUE AMBAS DIAGONAIS AFASTADAS ESTAO PRESENTES
		for(int i = afast + 2; i < MAX - afast - 2; i++)
			m.X[i] = m.B[i]/m.Dp[i] - m.Dia[i]*m.X[i - afast - 2]/m.Dp[i] - m.Di[i]*m.X[i - 1]/m.Dp[i] - m.Ds[i]*Xold[i + 1]/m.Dp[i] - m.Dsa[i]*Xold[afast + 2 + i]/m.Dp[i];
		
		//LOOP EM QUE DIAGONAL SUPERIROR AFASTADA ACABOU
		for(int i = MAX - afast - 2; i < MAX - 1; i++)
			m.X[i] = m.B[i]/m.Dp[i] - m.Dia[i]*m.X[i - afast - 2]/m.Dp[i] - m.Di[i]*m.X[i - 1]/m.Dp[i] - m.Ds[i]*Xold[i + 1]/m.Dp[i];
		
		//ULTIMA ATRIBUICAO
		m.X[MAX - 1] = m.B[MAX - 1]/m.Dp[MAX - 1] - m.Dia[MAX - 1]*m.X[MAX - afast - 3]/m.Dp[MAX - 1] - m.Di[MAX - 1]*m.X[MAX - 2]/m.Dp[MAX - 1];
		
		r = normaL2(m.X, &m, MAX, afast);

		printf("%1.15f\n", r);

		//X(ANTERIOR) = X(ATUAL)
		for(int i = 0; i < MAX; i++)
			Xold[i] = m.X[i];
		if(r < 0.0001)
			break;
		j++;
		/*
		if((j % 10) == 0){
			printf("iteração %d\n", j);
			for(int i = 0; i < MAX; i++)
				printf("%1.15f  ", m.X[i]);
			printf("\n\n");

		}*/
	}
}