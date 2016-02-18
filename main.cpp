#include <iostream>
#include <math.h>
#include "global_vars.h"
#include "fluid_class.h"
#include "initialize.h"
#include "input_output.h"
#include "prim_to_cons.h"
#include "flux_update.h"

using namespace std;

void init(int, int, double, double, fluid*);
void read_input_file(int &, double &, double &, int &, int &, double &); 
void prim_2_cons(int, int, fluid*, conservative*);
void compute_fluxes(int, int, double, double, fluid*, conservative*, conservative*); 
void update(int, int, conservative*, conservative*, double, double, double);
void compute_dt(int, int, fluid*, double, double, double, double &);
void output_results(int, int, double, double, fluid*);

int main()
{

 int nmax, imax, jmax, n; 
 double dx, dy, cfl, dt;
 
 // read input file
 read_input_file(nmax, dx, dy, imax, jmax, cfl); 

 fluid *f = new fluid [imax*jmax];
 
 // initialize
 init(imax,jmax,dx,dy,f);

 conservative *Q = new conservative [imax*jmax];

 prim_2_cons(imax,jmax,f,Q);

 conservative *dflux = new conservative [imax*jmax];
 
 // evolve solution
 for (n=0; n<nmax; n++)
 {
  compute_dt(imax,jmax,f,cfl,dx,dy,dt);
  cout<<"n = "<<n<<" dt = "<<dt<<endl;
  compute_fluxes(imax,jmax,dx,dy,f,Q,dflux); 
  update(imax,jmax,Q,dflux,dt,dx,dy);
  cons_2_prim(imax,jmax,f,Q);
 } 

 // output results 
 output_results(imax,jmax,dx,dy,f);


 delete [] f;
 delete [] Q;
 delete [] dflux;


 cout<<"done \n";
 return 0;
}

