#include "gs.h"
#include "sistLinear.h"
#include <math.h>
#include <stdio.h> //RETIRAR DEPOIS
#include <sys/time.h>

// Retorna tempo em milisegundos
double timestamp(void)
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}


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
	double aux;
	//PRIMEIRA ATRIBUICAO
		aux = (A->B[0] - (A->Dp[0]*Xnew[0] + A->Ds[0]*Xnew[1] + A->Dsa[0]*Xnew[afast + 2]));
		r += aux*aux; //r = (B[0] - A[0][-]*X)^2 

		//LOOP EM QUE DIAGONAL INFERIOR AFASTADA AINDA NAO APARECEU
		for(int i = 1; i < afast + 2; i++){
			aux = (A->B[i] - (A->Dp[i]*Xnew[i] + A->Di[i]*Xnew[i - 1] + A->Ds[i]*Xnew[i + 1] + A->Dsa[i]*Xnew[afast + 2 + i]));
			r += aux*aux;
		}	
		//LOOP EM QUE AMBAS DIAGONAIS AFASTADAS ESTAO PRESENTES
		for(int i = afast + 2; i < MAX - afast - 2; i++){
			aux = (A->B[i] - (A->Dp[i]*Xnew[i] + A->Dia[i]*Xnew[i - afast - 2] + A->Di[i]*Xnew[i - 1] + A->Ds[i]*Xnew[i + 1] + A->Dsa[i]*Xnew[afast + 2 + i]));
			r += aux*aux;
		}
		//LOOP EM QUE DIAGONAL SUPERIROR AFASTADA ACABOU
		for(int i = MAX - afast - 2; i < MAX - 1; i++){
			aux = (A->B[i] - (A->Dp[i]*Xnew[i] + A->Dia[i]*Xnew[i - afast - 2] + A->Di[i]*Xnew[i - 1] + A->Ds[i]*Xnew[i + 1]));
			r += aux*aux;
		}
		//ULTIMA ATRIBUICAO
		aux = (A->B[MAX - 1] - (A->Dp[MAX - 1]*Xnew[MAX - 1] + A->Dia[MAX - 1]*Xnew[MAX - afast - 3] + A->Di[MAX - 1]*Xnew[MAX - 2]));
		r += aux*aux;
	return sqrt(r);// retorna raiz de (r[0]^2+r[1]^2+...+r[n]^2)
}


//FUNCAO QUE CALCULA O VETOR RESULTADO DA ESTRUTURA ATRAVES DO METODO DE
//GAUSS-SEIDEL, ELIMINANDO DO CALCULO AS OCORRENCIAS DE ZEROS
void Gs(struct matPent *m, int afast, int MAX, int ite){
	double Xold[MAX], Taux;
	unsigned int j = 0;
	
	for(int i = 0; i < MAX; i++)
		Xold[i] = 0;
	while(j < ite){
		
		Taux = timestamp();
		//PRIMEIRA ATRIBUICAO
		m->X[0] = m->B[0]/m->Dp[0] - m->Ds[0]*Xold[1]/m->Dp[0] - m->Dsa[0]*Xold[afast + 2]/m->Dp[0]; 

		//LOOP EM QUE DIAGONAL INFERIOR AFASTADA AINDA NAO APARECEU
		for(int i = 1; i < afast + 2; i++)
			m->X[i] = m->B[i]/m->Dp[i] - m->Di[i]*m->X[i - 1]/m->Dp[i] - m->Ds[i]*Xold[i + 1]/m->Dp[i] - m->Dsa[i]*Xold[afast + 2 + i]/m->Dp[i];

		//LOOP EM QUE AMBAS DIAGONAIS AFASTADAS ESTAO PRESENTES
		for(int i = afast + 2; i < MAX - afast - 2; i++)
			m->X[i] = m->B[i]/m->Dp[i] - m->Dia[i]*m->X[i - afast - 2]/m->Dp[i] - m->Di[i]*m->X[i - 1]/m->Dp[i] - m->Ds[i]*Xold[i + 1]/m->Dp[i] - m->Dsa[i]*Xold[afast + 2 + i]/m->Dp[i];
		
		//LOOP EM QUE DIAGONAL SUPERIROR AFASTADA ACABOU
		for(int i = MAX - afast - 2; i < MAX - 1; i++)
			m->X[i] = m->B[i]/m->Dp[i] - m->Dia[i]*m->X[i - afast - 2]/m->Dp[i] - m->Di[i]*m->X[i - 1]/m->Dp[i] - m->Ds[i]*Xold[i + 1]/m->Dp[i];
		
		//ULTIMA ATRIBUICAO
		m->X[MAX - 1] = m->B[MAX - 1]/m->Dp[MAX - 1] - m->Dia[MAX - 1]*m->X[MAX - afast - 3]/m->Dp[MAX - 1] - m->Di[MAX - 1]*m->X[MAX - 2]/m->Dp[MAX - 1];
		
		
		m->T =+ (timestamp() - Taux)/ite; //tempo/maxIte

		m->r[j] = normaL2(m->X, m, MAX, afast);
		//m->r[j] = normaAprox(Xold, m->X, MAX);
		
		//X(ANTERIOR) = X(ATUAL)
		for(int i = 0; i < MAX; i++)
			Xold[i] = m->X[i];
		j++;
	}
}