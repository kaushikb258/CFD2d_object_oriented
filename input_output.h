#ifndef H_INPUT_OUTPUT
#define H_INPUT_OUTPUT

#include <iostream>
#include <fstream>
#include <string>
#include "global_vars.h"
#include "out_vtk.h"


using namespace std;

void read_input_file(int &nmax, double &dx, double &dy, int &imax, int &jmax, double &cfl)
{

// read input file
string s;
ifstream f;
f.open ("input_file", ios_base::in);
f>>s>>nmax;
f>>s>>dx;
f>>s>>dy;
f>>s>>imax;
f>>s>>jmax;
f>>s>>cfl;
f.close();

cout<<"------------"<<endl;
cout<<"nmax = "<<nmax<<endl;
cout<<"dx = "<<dx<<endl;
cout<<"dy = "<<dy<<endl;
cout<<"imax = "<<imax<<endl;
cout<<"jmax = "<<jmax<<endl;
cout<<"CFL = "<<cfl<<endl;
cout<<"------------"<<endl;

}


void output_results(int imax, int jmax, double dx, double dy, fluid *f)
{

 string filename = "output_paraview.vtk";
 string title = "cfd2d";
 int i, j, k;
 int nprim = 4;
 double coord[imax][jmax][3];
 double q[imax][jmax][4];

 // initialize coordinates
 for (j=0; j<jmax; j++)
 {
  for (i=0; i<imax; i++)
  {
   coord[i][j][0] = ((double) (i) + 0.5)*dx; //x
   coord[i][j][1] = ((double) (jmax-1-j) + 0.5)*dy; //y
   coord[i][j][2] = 0.0; //z = 0 as 2D grid
  }
 }

 // q
 for (j=0; j<jmax; j++)
 {
  for (i=0; i<imax; i++)
  {
   k = i + j*imax;
   q[i][j][0] = f[k].rho;
   q[i][j][1] = f[k].u;
   q[i][j][2] = f[k].v;
   q[i][j][3] = f[k].p;
  }
 }

 // write vtk file 
 vtk_write(filename,title,imax,jmax,&coord[0][0][0],nprim,&q[0][0][0]);

}


#endif
