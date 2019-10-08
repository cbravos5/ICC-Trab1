#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "sistLinear.h"
#include "fronteira.h"



//FUNCAO QUE PREENCHE OS VETORES ALOCADOS DA ESTRUTURA COM SEUS VALORES CORRESPONDENTES
void PreencheVetores(int nx, int ny, double hx, double hy, struct matPent *A)
{
	double hx2, hy2;
	hx2 = hx*hx;
	hy2 = hy*hy;
	int max = (nx-2)*(ny-2);

	//DIAGONAL PRINCIPAL
	double CTE = 4*(hx2 + hy2 + 2*M_PI*M_PI*hx2*hy2);
	for(int i = 0; i < max; i++)
	{
		A->Dp[i] = CTE;
	}
	
	//DIAGONAL SUPERIOR
	CTE = hx*hy2 - 2*hy2;
	for (int i = 0; i < max - 1 ; ++i)
	{
		A->Ds[i] = CTE;
	}
	A->Ds[max-1] = 0;

	//DIAGONAL INFERIOR
	CTE = -2*hy2 - hx*hy2;
	for (int i = 1; i < max; ++i)
	{
		A->Di[i] = CTE;	
	}
	A->Di[0] = 0;

	//DIAGONAL SUPERIOR AFASTADA
	CTE = hy*hx2 - 2*hx2;
	for (int i = 0; i < max; ++i)
	{
		if(i < max-nx)
			A->Dsa[i] = CTE;
		else
			A->Dsa[i] = 0;
	}

	//DIAGONAL INFERIOR AFASTADA
	CTE = -2*hx2 - hy*hx2;
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


//FUNCAO QUE LIDA COM AS CONDICOES DE CONTORNO
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
			A->B[p] = 2*hx2*hy2*fX((i+1)*hx,(j+1)*hy);
			//ZEROS DA DIAGONAL INFERIOR
			if(i == 0 && j > 0)
				A->Di[p] = 0; 
			//ZEROS DA DIAGONAL SUPERIOR
			if(i == nx-3 && j < ny-3)
				A->Ds[p] = 0;
			//CONTORNO SUL
			if(j == 0)
				A->B[p] += (-2*hx2 - hy*hx2)*u_fronteira_sul(i*hx);
			//CONTORNO NORTE
			if(j == ny-3)
				A->B[p] -= (hy*hx2 - 2*hx2)*u_fronteira_norte(i*hx); 
		}
	}
}

//FUNCAO QUE IMPRIME O VETOR RESULTADO E AS ITERACOES NA SAIDA PADRAO OU DENTRO DE UM ARQUIVO
void out(struct matPent A, int nx, int ny, double hx, double hy, int ite, char *arq, int output){
	FILE *pF;
	if(output)
		pF = fopen(arq, "w");
	else
		pF = stdout;
	int i = 0;
	fputs("##########\n",pF);
	fputs("# Tempo método GS: ",pF);
	fprintf(pF, "%1.15f\n#\n", A.T);
	fputs("# Norma L2 do resíduo\n",pF);
	for(i = 0; i < ite; i++){
		fprintf(pF, "# i = %d: %1.15f\n", i, A.r[i]);
	}
	fputs("##########\n",pF);
	for(int j = 0; j < (ny-2); j++){
		for(i = 0; i < (nx-2); i++){
			double ihx, jhy;
			ihx = (i+1)*hx;
			jhy = (j+1)*hy; 
			fprintf(pF, "%f %f %f\n", ihx, jhy, A.X[(i + j*nx-2)]);
		}
	}
	fputs("##########\n",pF);
	fclose(pF);
}