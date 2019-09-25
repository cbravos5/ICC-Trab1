#include <stdio.h>
#include <math.h>

#define PI 3.14159265

double u_fronteira_sul(double x){
	/*double a,b;
	a = 2*PI*PI - 2*PI*x;
	b = PI*PI;
	a = sin(a);
	b = sinh(b);
	return(a*b);*/
	return(sin((2*M_PI*M_PI) - (2*M_PI*x))*sinh(M_PI*M_PI));

}

double u_fronteira_norte(double x){
	return(sin(2*M_PI*x)*sinh(M_PI*M_PI));
}

x[0] = (b[0] - dsup[0]*x[1])/dprin[0];
for(int i = 1; i < TAM - 1; i++){
	x[i] = (b[i] - dinf[i]*x[i - 1] - dsup[i]*x[i + 1])dprin[i];
}
x[TAM] = (b[TAM] - dinf[TAM]*x[TAM - 1]);

void gs(struct matPent m, int afast, double eps){
	double Xold[MAX];
	double e = eps + 1.0;
	for(int i = 0; i < MAX; i++)
		Xold[i] = 0;
	while(e > eps){
		m.X[0] = (m.B[0] - m.Ds[0]*m.X[1] - m.Dsa[0]*m.X[afast])/m.Dp[0];
		for(int i = 1; i < afast - 1; i++){
			m.X[i] = (m.B[i] - m.Di[i]*m.X[i - 1] - m.Ds[i]*m.X[i + 1] - m.Dsa[i]*m.X[i + afast])/m.Dp[i];
		}
		for(int i = afast; i < MAX - afast; i++)
			m.X[i] = (m.B[i] - m.Dia[i]*m.X[i - afast] - m.Di[i]*m.X[i - 1] - m.Ds[i]*m.X[i + 1] - m.Dsa[i]*m.X[i + afast])/m.Dp[i];
		for(int i = MAX - afast; i < MAX - 1; i++)
			m.X[i] = (m.B[i] - m.Dia[i]*m.X[i - afast] - m.Di[i]*m.X[i - 1] - m.Ds[i]*m.X[i + 1])/m.Dp[i];
		m.X[i] = (m.B[i] - m.Dia[i]*m.X[i - 3] - m.Di[i]*m.X[i - 1])/m.Dp[i];
		e = X[0] - Xold[0];
		for(int i = 1; i < MAX; i++){
			if(fabs(X[i] - Xold [i]) > e)
				e = fabs(X[i] - Xold[i]);
		}
	}
}

void main(){
	printf("%lf %lf\n",u_fronteira_sul(0.00324), u_fronteira_norte(0.00324));
}