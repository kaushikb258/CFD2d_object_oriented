#ifndef H_FLUX_UPDATE
#define H_FLUX_UPDATE

#include "global_vars.h"
#include "riemann.h"

void MP5(double, double, double, double, double, double&);
double minmod(double, double);
double median(double, double, double);


void compute_fluxes(int imax, int jmax, double dx, double dy, fluid *f, conservative *Q, conservative *dflux) 
{

 int i, j, k;
 double left[4], right[4];
 double fmh[4], fph[4], gmh[4], gph[4]; 
 int im3, im2, im1, iz, ip1, ip2, ip3;
 double sim3, sim2, sim1;
 int jm3, jm2, jm1, jz, jp1, jp2, jp3;
 double sjm3, sjm2, sjm1;


 for (j=0; j<jmax; j++)
 {
  for (i=0; i<imax; i++)
  { 
   k = i + j*imax;  
   dflux[k].q1 = 0.0;
   dflux[k].q2 = 0.0;
   dflux[k].q3 = 0.0;
   dflux[k].q4 = 0.0;
  } 
 } 
 
//-----------

 for (j=0; j<jmax; j++)
 {
  for (i=0; i<imax; i++)
  {
   k = i + j*imax;  

//----
// i+1/2 

   im2 = k-2;
   im1 = k-1;
   iz = k;
   ip1 = k+1;
   ip2 = k+2;
   ip3 = k+3;

   sim2 = 1.0;
   sim1 = 1.0;   

   if(i==0) { im2 = k+1; im1 = k; sim2 = -1.0; sim1 = -1.0;}
   if(i==1) { im2 = k-1; sim2 = -1.0;}
   if(i==imax-3) { ip3 = k+2;}
   if(i==imax-2) { ip2 = k+1; ip3 = k+1;}
   if(i==imax-1) { ip1 = k; ip2 = k; ip3 = k;}

   MP5(f[im2].rho,f[im1].rho,f[iz].rho,f[ip1].rho,f[ip2].rho,left[0]);
   MP5(sim2*f[im2].u,sim1*f[im1].u,f[iz].u,f[ip1].u,f[ip2].u,left[1]);
   MP5(f[im2].v,f[im1].v,f[iz].v,f[ip1].v,f[ip2].v,left[2]);
   MP5(f[im2].p,f[im1].p,f[iz].p,f[ip1].p,f[ip2].p,left[3]);
   
   MP5(f[ip3].rho,f[ip2].rho,f[ip1].rho,f[iz].rho,f[im1].rho,right[0]);
   MP5(f[ip3].u,f[ip2].u,f[ip1].u,f[iz].u,sim1*f[im1].u,right[1]);
   MP5(f[ip3].v,f[ip2].v,f[ip1].v,f[iz].v,f[im1].v,right[2]);
   MP5(f[ip3].p,f[ip2].p,f[ip1].p,f[iz].p,f[im1].p,right[3]);


   riemann_solver(1,left,right,fph);

//----
// i-1/2 
  
 
   im3 = k-3;
   im2 = k-2;
   im1 = k-1;
   iz = k;
   ip1 = k+1;
   ip2 = k+2;

   sim3 = 1.0;
   sim2 = 1.0;
   sim1 = 1.0;

   if(i==0) { im3 = k+2; im2 = k+1; im1 = k; sim3 = -1.0; sim2 = -1.0; sim1 = -1.0;}
   if(i==1) { im3 = k; im2 = k-1; sim3 = -1.0; sim2 = -1.0;}
   if(i==2) {im3 = k-2; sim3 = -1.0;}
   if(i==imax-2) { ip2 = k+1;}
   if(i==imax-1) { ip1 = k; ip2 = k;}

   MP5(f[im3].rho,f[im2].rho,f[im1].rho,f[iz].rho,f[ip1].rho,left[0]);
   MP5(sim3*f[im3].u,sim2*f[im2].u,sim1*f[im1].u,f[iz].u,f[ip1].u,left[1]);
   MP5(f[im3].v,f[im2].v,f[im1].v,f[iz].v,f[ip1].v,left[2]);
   MP5(f[im3].p,f[im2].p,f[im1].p,f[iz].p,f[ip1].p,left[3]);

   MP5(f[ip2].rho,f[ip1].rho,f[iz].rho,f[im1].rho,f[im2].rho,right[0]);
   MP5(f[ip2].u,f[ip1].u,f[iz].u,f[im1].u,sim1*f[im2].u,right[1]);
   MP5(f[ip2].v,f[ip1].v,f[iz].v,f[im1].v,f[im2].v,right[2]);
   MP5(f[ip2].p,f[ip1].p,f[iz].p,f[im1].p,f[im2].p,right[3]);


   riemann_solver(1,left,right,fmh);

//----
// j+1/2 

   jm2 = k+2*imax;
   jm1 = k+imax;
   jz = k;
   jp1 = k-imax;
   jp2 = k-2*imax;
   jp3 = k-3*imax;

   sjm1 = 1.0;
   sjm2 = 1.0;  

   if(j==0) { jp1 = k; jp2 = k; jp3 = k;}
   if(j==1) { jp2 = k-imax; jp3 = k-imax;}
   if(j==2) { jp3 = k-2*imax;}
   if(j==jmax-2) { jm2 = k+imax; sjm2 = -1.0;}
   if(j==jmax-1) { jm1 = k; jm2 = k-imax; sjm1 = -1.0; sjm2 = -1.0;}

   MP5(f[jm2].rho,f[jm1].rho,f[jz].rho,f[jp1].rho,f[jp2].rho,left[0]);
   MP5(f[jm2].u,f[jm1].u,f[jz].u,f[jp1].u,f[jp2].u,left[1]);
   MP5(sjm2*f[jm2].v,sjm1*f[jm1].v,f[jz].v,f[jp1].v,f[jp2].v,left[2]);
   MP5(f[jm2].p,f[jm1].p,f[jz].p,f[jp1].p,f[jp2].p,left[3]);
   
   MP5(f[jp3].rho,f[jp2].rho,f[jp1].rho,f[jz].rho,f[jm1].rho,right[0]);
   MP5(f[jp3].u,f[jp2].u,f[jp1].u,f[jz].u,f[jm1].u,right[1]);
   MP5(f[jp3].v,f[jp2].v,f[jp1].v,f[jz].v,sjm1*f[jm1].v,right[2]);
   MP5(f[jp3].p,f[jp2].p,f[jp1].p,f[jz].p,f[jm1].p,right[3]);


   riemann_solver(2,left,right,gph);

//----
// j-1/2 

   jm3 = k+3*imax;
   jm2 = k+2*imax;
   jm1 = k+1*imax;
   jz = k;
   jp1 = k-1*imax;
   jp2 = k-2*imax;

   sjm1 = 1.0;
   sjm2 = 1.0;  
   sjm3 = 1.0;

   if(j==0) { jp1 = k; jp2 = k;}
   if(j==1) { jp2 = k-imax;}
   if(j==jmax-3) {jm3 = k+2*imax; sjm3 = -1.0;}
   if(j==jmax-2) { jm3 = k; jm2 = k+imax; sjm3 = -1.0; sjm2 = -1.0;}
   if(j==jmax-1) { jm3 = k-2*imax; jm2 = k-imax; jm1 = k; sjm3 = -1.0; sjm2 = -1.0; sjm1 = -1.0;}

   MP5(f[jm3].rho,f[jm2].rho,f[jm1].rho,f[jz].rho,f[jp1].rho,left[0]);
   MP5(f[jm3].u,f[jm2].u,f[jm1].u,f[jz].u,f[jp1].u,left[1]);
   MP5(sjm3*f[jm3].v,sjm2*f[jm2].v,sjm1*f[jm1].v,f[jz].v,f[jp1].v,left[2]);
   MP5(f[jm3].p,f[jm2].p,f[jm1].p,f[jz].p,f[jp1].p,left[3]);
   
   MP5(f[jp2].rho,f[jp1].rho,f[jz].rho,f[jm1].rho,f[jm2].rho,right[0]);
   MP5(f[jp2].u,f[jp1].u,f[jz].u,f[jm1].u,f[jm2].u,right[1]);
   MP5(f[jp2].v,f[jp1].v,f[jz].v,sjm1*f[jm1].v,sjm2*f[jm2].v,right[2]);
   MP5(f[jp2].p,f[jp1].p,f[jz].p,f[jm1].p,f[jm2].p,right[3]);


   riemann_solver(2,left,right,gmh);

//----

  dflux[k].q1 += (fph[0] - fmh[0])/dx;
  dflux[k].q2 += (fph[1] - fmh[1])/dx;
  dflux[k].q3 += (fph[2] - fmh[2])/dx;
  dflux[k].q4 += (fph[3] - fmh[3])/dx;

  dflux[k].q1 += (gph[0] - gmh[0])/dy;
  dflux[k].q2 += (gph[1] - gmh[1])/dy;
  dflux[k].q3 += (gph[2] - gmh[2])/dy;
  dflux[k].q4 += (gph[3] - gmh[3])/dy;

  }
 }  

}


void update(int imax, int jmax, conservative *Q, conservative *dflux, double dt, double dx, double dy)
{
 int i, j, k; 
 for (j=0; j<jmax; j++)
 {
  for (i = 0; i<imax; i++)
  {
   k = i + j*imax;  
   Q[k].q1 = Q[k].q1 - dt*dflux[k].q1;
   Q[k].q2 = Q[k].q2 - dt*dflux[k].q2;
   Q[k].q3 = Q[k].q3 - dt*dflux[k].q3;
   Q[k].q4 = Q[k].q4 - dt*dflux[k].q4;
  }
 } 
}


void MP5(double s1, double s2, double s3, double s4, double s5, double &sph)
{
 double sorig = (2.0*s1 -13.0*s2 + 47.0*s3 + 27.0*s4 - 3.0*s5)/60.0; 
 const double alpha = 1.5;
 double sa, sb;
 sa = s4-s3;
 sb = alpha*(s3-s2); 
 double smp = s3 + minmod(sa,sb);
 sph = median(sorig,s3,smp); 
}

double minmod(double a, double b)
{
 if (a*b < 0) return 0.0;
 else if(abs(a) < abs(b)) return a;
 else return b;
}

double median(double a, double b, double c)
{
 if(b<=a && c>=a) return a;
 if(c<=a && b>=a) return a;
 if(a<=b && c>=b) return b;
 if(c<=b && a>=b) return b;
 if(a<=c && b>=c) return c;
 if(b<=c && a>=c) return c;
 cout<<"should not be here \n";
}


#endif

