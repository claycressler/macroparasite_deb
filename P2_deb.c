#include <R.h>
#include <stdio.h>

static double parms[28];

#define imax parms[0]
#define theta parms[1]
#define eta parms[2]
#define rho parms[3]
#define epsA parms[4]
#define epsG parms[5]
#define epsR parms[6]
#define epsI parms[7]
#define Linf parms[8]
#define alpha parms[9]
#define gmin parms[10]
#define m parms[11]
#define mc parms[12]
#define mi parms[13]
#define k parms[14]
#define b parms[15]
#define ui parms[16]
#define umin parms[17]
#define uacc parms[18]
#define epsAmin parms[19]
#define heps parms[20]
#define epsP parms[21]
#define sigmaC parms[22]
#define hC parms[23]
#define uP parms[24]
#define uC parms[25]
#define uI parms[26]
#define InfDose parms[27]
#define MAX(a,b) (((a)>(b))?(a):(b))

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {
  int N=28;
  odeparms(&N, parms);
}

/* derivatives */
void derivs (int *neq, double *t, double *y, double *ydot) {

  // state variables
  double G = y[0];
  double C = y[1];
  double S = y[2];
  double R = y[3];
  double Ii = y[4];
  double P2 = y[6];

  // condition-dependent ingestion rate
  double In = imax * pow(S, 0.66666667) / (1 + exp(eta * (R/S - theta)));
  //Rprintf("\n condition %lg, ingestion %lg", eta * (R/S - theta), exp(eta * (R/S-theta)));
  // parasite-dependent assimilation efficiency
  double eA = epsA * (1 - epsAmin * P2 / (heps + P2));
  // growth rate
  double Grate = 3 * gmin * (pow(alpha, 0.33333333) * Linf * pow(S, 0.66666667) - S);
  // cost of growth
  double CG = epsG * Grate;
  // constitutive defense
  double Ic = k * (S + R);
  // maintenance rate
  double M = (m + mc * Ic)*(S + R) + mi * Ii;

  // rate equations
  ydot[0] = In - rho * G;
  ydot[1] = rho * (1 - eA) * G - rho * C - sigmaC * P2 * C/(hC + C);
  ydot[2] = Grate;
  ydot[3] = epsR * (rho * eA * G - M - CG) - b * R * P2;
  ydot[4] = epsI * b * R * P2 - ui * Ii;
  ydot[5] = umin + uacc * MAX(theta * S/R - 1, 0);
  ydot[6] = epsP * sigmaC*P2*C/(hC+C) - uP * P2 - uC * Ic * P2 - uI * Ii * P2;
}

/* At some point the parasite invades and the parasite abundance goes from 0 to something positive */
void event(int *n, double *t, double *y) {
  y[6] = InfDose;
}


