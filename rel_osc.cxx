#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0, const double t);
void RKstep(double* const yn, const double* const y0, const double t, const double dt, double* k1, double* k2, double* k3, double* k4);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
	const double dt = 0.01;
	double t = 0;
	const double T= 20;
	const int dim = 2;
	
	double y0[2], yn[2];
	for (double p0 = 0.1; p0 < 5; p0 += 0.1){
	  y0[0] = p0;
	  y0[1] = 0;
	  double k1[dim], k2[dim], k3[dim], k4[dim];
	  double b1, b2, b3, b4;
	  t = 0;
	  while(t<=T)
	  {
		  t += dt;
		  RKstep(yn, y0, t, dt, k1,k2,k3,k4);
		  if (y0[1]>0 && yn[1]<0)
		    break;
		for(int i=0; i<dim; i++) y0[i] = yn[i];
	  
  // 	  out << t << "\t" << y0[0] << "\t" << y0[1] << endl;
	  }
	  
	  
	  double thetal = 0;
	  double thetar = 1;
	  double tol = 1e-8;
	  double thetam;
	  while (tol < abs(yn[1])){
	  thetam = (thetal + thetar)/2.0;
	  
	  b1 = thetam - (3. * thetam * thetam)/2.0 + ( 2. * thetam * thetam * thetam ) / 3.;
	  b2 = thetam * thetam - ( 2. * thetam * thetam * thetam ) / 3.;
	  b3 = thetam * thetam - ( 2. * thetam * thetam * thetam ) / 3.;
	  b4 = - (thetam * thetam) / 2. + ( 2. * thetam * thetam * thetam ) / 3.;
	  
	  yn[1] = y0[1] + dt * b1 * k1[1] + b2 * k2[1] + b3 * k3[1] + b4 * k4[1];
	  
	  if (yn[1] > 0)
	    thetal = thetam;
	  else if (yn[1] <= 0)
	    thetar = thetam;
	  
  // 	cout << thetam << "\t" << yn[1] << endl;
	  }
	  
	  out << p0 << "\t" << t - dt + thetam*dt << endl;
	}
	out.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0, const double t, const double dt, double* k1, double* k2, double* k3, double* k4)
{
	const int dim = 2;
	

        for(int i=0; i<dim; i++) k1[i] = y0[i];
	f(k1, t);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dt * k1[i];
         f(k2, t+0.5*dt);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dt * k2[i];
	f(k3, t+0.5*dt);

        for(int i=0;i<dim; i++) k4[i] = y0[i] + dt * k3[i];
	f(k4,  t+dt);

	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + 1./6.*dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}
//-------------------

void f(double* const y0, const double t)
{	
	double y[2] = { y0[0], y0[1]};
	
	y0[0] = y[1];
	y0[1] = -y[0]/(sqrt(1+y[0]*y[0]));
	
}
