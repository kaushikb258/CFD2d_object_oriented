#ifndef H_FLUID_CLASS
#define H_FLUID_CLASS

#include <math.h>
#include "global_vars.h"


class fluid
{
 private:
 
 public: 

  double rho;
  double u;
  double v;
  double e;
  double p;
  double T;
  double c;

  
  // constructor
  fluid()
  {
   rho = 1.2;
   u = 0.0;
   v = 0.0;
   p = 1.0e5;
   eos();
  }  
  
  // destructor
  ~fluid() {}

  // equation of state
  void eos()
  {
   e = p/rho/(gam-1.0);
   T = p/rho/Rgas;
   c = sqrt(gam*Rgas*T);
  }

  void get_p_and_T()
  {
   p = (gam-1.0)*rho*e;
   T = p/rho/Rgas; 
   c = sqrt(gam*p/rho);
  }

};


class conservative
{
 private:

 public:

  double q1;
  double q2;
  double q3; 
  double q4;

  conservative()
  {
   q1 = 0.0;
   q2 = 0.0;
   q3 = 0.0;
   q4 = 0.0;
  }
 
  ~conservative() {}

};



#endif


