#ifndef H_INITIALIZE
#define H_INITIALIZE

#include <math.h>
#include "global_vars.h"

void init(int imax, int jmax, double dx, double dy, fluid *f)
{
 int i, j, k;
 double x, y, r;

 for (j=0; j<jmax; j++)
 {
  for (i=0; i<imax; i++)
  {
   k = i + j*imax;

   x = ((double)(i)+0.5)*dx;   
   y = ((double)(jmax-1-j)+0.5)*dy;   
   r = sqrt(x*x + y*y);
   
   f[k].u = 0.0;
   f[k].v = 0.0;

   if(r < 0.2)
   {
   // driver
   f[k].rho = 12.0;
   f[k].p = 10.0e5;
   }
   else
   {
   // driven
   f[k].rho = 1.2;
   f[k].p = 1.0e5;
   }
   f[k].eos();
  }
 } 
}


void compute_dt(int imax, int jmax, fluid *f, double cfl, double dx, double dy, double &dt)
{
 dt = 0.0;
 double vel, vel_max, dist;
 int i, j, k;
 
 dist = (dx < dy)? dx : dy;

 vel_max = 0.0;
 for (j=0; j<jmax; j++)
 {
 for (i=0; i<imax; i++)
 {
  k = i + j*imax;
  vel = f[k].c + sqrt(pow(f[k].u,2.0) + pow(f[k].v,2.0));
  vel_max = (vel_max > vel)? vel_max : vel;
 }
 }

 dt = cfl*dist/vel_max;

}


#endif

