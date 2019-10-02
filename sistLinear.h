#ifndef SISTLINEAR_H_INCLUDE
#define SISTLINEAR_H_INCLUDE

struct matPent
{
	double *Dp, *Ds, *Di, *Dsa, *Dia, *X, *B, *r, T;
};

void PreencheVetores(int nx, int ny, double hx, double hy, struct matPent *A);

void printVet(int nx, int ny, struct matPent *A);

void contorno(int nx, int ny, double hx, double hy, struct matPent *A);

void out(struct matPent A, int nx, int ny, double hx, double hy, char *arq, int output)

#endif