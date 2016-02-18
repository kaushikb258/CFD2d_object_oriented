#ifndef H_OUT_VTK
#define H_OUT_VTK

#include <iostream>
#include <string>
#include <sstream> 
#include <fstream>
#include <math.h>

using namespace std;

void vtk_write(string filename,string title,int imax, int jmax, double* xyz, int nprim, double* prim)
{
 string cell_size_string, node_num_string; 
 int node;
 int i, j, k, n;
 stringstream s_node_num, s_cells, s_imax, s_jmax, s_kmax; 
 ofstream f;
 int kmax = 1;

 s_node_num << (imax*jmax);
 s_cells << (imax-1)*(jmax-1);
 s_imax << imax;
 s_jmax << jmax;
 s_kmax << kmax;

 f.open (filename.c_str(),ios_base::out);
 f<< "# vtk DataFile Version 2.0\n";
 f<< title<<"\n";
 f<< "ASCII\n";
 f<< "DATASET STRUCTURED_GRID\n";
 f<< "DIMENSIONS "<<"\t"<<s_imax.str()<<"\t\t"<<s_jmax.str()<<"\t\t"<<s_kmax.str()<<"\n"; 
 f<< "POINTS "<<"\t"<<s_node_num.str()<<"\t"<<"double\n";
 
 n = 0;
 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {
   for (k=0; k<3; k++)
   { 
    f<<xyz[n]<<" ";
    n++;
   }
   f<<"\n";
  }
 }  

 f<< "CELL_DATA "<<"\t"<<s_cells.str()<<"\n";
 f<< "POINT_DATA "<<"\t"<<s_node_num.str()<<"\n";

 f<< "SCALARS density double \n";
 f<< "LOOKUP_TABLE default \n";
 n = 0;
 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {
   f<< prim[n]<<" ";
   n += nprim;
  }
 } 
 f<<"\n";

 f<< "SCALARS u-velocity double \n";
 f<< "LOOKUP_TABLE default \n";
 n = 1;
 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {
   f<< prim[n]<<" ";
   n += nprim;
  }
 } 
 f<<"\n";

 f<< "SCALARS v-velocity double \n";
 f<< "LOOKUP_TABLE default \n";
 n = 2;
 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {
   f<< prim[n]<<" ";
   n += nprim;
  }
 } 
 f<<"\n";

 f<< "SCALARS pressure double \n";
 f<< "LOOKUP_TABLE default \n";
 n = 3;
 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {
   // pressure in bar
   f<< prim[n]/(1.0e5)<<" ";
   n += nprim;
  }
 } 

}

#endif
