#ifndef H_PRIM_TO_CONS
#define H_PRIM_TO_CONS

#include <math.h>
#include "global_vars.h"

void prim_2_cons(int imax, int jmax, fluid* f, conservative* Q)
{
 double ke, te;
 int i, j, k;

 for (j=0; j<jmax; j++)
 {
 for (i=0; i<imax; i++)
 {
  k = i + j*imax;

  Q[k].q1 = f[k].rho;
  Q[k].q2 = f[k].rho*f[k].u;
  Q[k].q3 = f[k].rho*f[k].v;
  
  ke = 0.5*(pow(f[k].u,2.0) + pow(f[k].v,2.0)); 
  te = ke + f[k].e;
  Q[k].q4 = f[k].rho*te; 
 }
 } 
 
}


void cons_2_prim(int imax, int jmax, fluid* f, conservative* Q)
{
 double ke;
 int i, j, k;
 
 for (j=0; j<jmax; j++)
 {
 for (i=0; i<imax; i++)
 {
  k = i + j*imax;

  f[k].rho = Q[k].q1;
  f[k].u = Q[k].q2/f[k].rho;
  f[k].v = Q[k].q3/f[k].rho;
  
  ke = 0.5*(pow(f[k].u,2.0) + pow(f[k].v,2.0)); 
  f[k].e = Q[k].q4/f[k].rho - ke; 
  f[k].get_p_and_T();
 }
 }
}


#endif

