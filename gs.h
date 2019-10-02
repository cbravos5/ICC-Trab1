#ifndef GS_H_INCLUDE
#define GS_H_INCLUDE

#include "sistLinear.h"

double normaL2(double *Xnew, struct matPent *A, int MAX, int afast);

double normaAprox(double *Xold, double *Xnew, int MAX);

void Gs(struct matPent *m, int afast, int MAX, int ite);

#endif