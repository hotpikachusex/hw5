#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std;

void eulerforw(const int N, double* x ,double* vx, const double dt);
void eulerbackw(const int N, double* y ,double* vy, const double dt);
void output(const int N, double* x, double* vx, double* y, double* vy, double* z, const double dt);
void cosinus(const int N, double* z, const double dt);

       int main(){
 
 	  
 	  const int N = 2000;
 	  const double tmax = 20.0*M_PI, tmin = 0.0; 
 	  const double dt = (tmax-tmin)/(N-1.0);
 	  double* x = new double[N];
 	  double* vx = new double[N];
	  double* y = new double[N];
 	  double* vy = new double[N];
	  double* z = new double[N];
 
 	  x[0] = 1.0;
 	  vx[0] = 0.0;
	  y[0] = 1.0;
	  vy[0] = 0.0; 
 
 
 	  eulerforw(N,x,vx,dt);
  	  eulerbackw(N,y,vy,dt);
 	  cosinus(N,z,dt);
	  output(N,x,vx,y,vy,z,dt);
 
 	  delete[] x;
 	  delete[] vx;
	  delete[] y;
 	  delete[] vy;
	  delete[] z;
 		  return 0;
   
 }
 
 void eulerforw(const int N, double* x, double* vx, const double dt){
         for(int i = 0; i<N; i++){
                 x[i+1] = x[i] + dt*vx[i];
                 vx[i+1] = vx[i]-dt*x[i];
     }
 }


 void eulerbackw(const int N, double* y, double* vy, const double dt){
         
         for(int i =0;i<N; i++){
                 y[i+1] = 1/(dt*dt+1.0) * y[i] + dt*vy[i];
                 vy[i+1] = 1/(dt*dt+1.0) * vy[i] - dt *y[i];
         }
 }
 
 void cosinus(const int N, double* z, const double dt){
         
         for(int i =0;i<N; i++){
                 z[i]=cos(i * dt);
         }
 }


void output(const int N, double* x, double* vx, double* y, double* vy, double* z, const double dt){
        ofstream out("Eulers.txt");

        for(int i = 0; i < N; i++){
                out << i * dt << "\t" << x[i] << "\t" << vx[i] << "\t" <<  y[i] << "\t" << vy[i] << "\t" << z[i] << endl;
        }

        out.close();
	

} 