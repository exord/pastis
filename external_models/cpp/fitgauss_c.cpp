//fun_c.cpp

#include "fitgauss_c.h"
#include <math.h>

void gaussian_res(double* p, double* xi, double *yi, int n, double* rez)
{   
   double p0   = p[0];
   double p1   = p[1];
   double p2   = p[2];
   double p3_2 = p[3] * p[3];
   for (int i = 0; i < n; i++)
     rez[i] = p0 + p1 * exp(-(xi[i] - p2) * (xi[i] - p2) / p3_2) - yi[i];
}
//                                                                        
void gaussian_res_J(double* p, double* xi, double *yi, int n, double* rez)
{
   double p1   = p[1];
   double p2   = p[2];
   double p3   = p[3];
   double p3_2 = p[3] * p[3];
   
   // new numpy array is row-major!  (http://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html)
   double* rez0 = rez;
   double* rez1 = rez + n;
   double* rez2 = rez + 2 * n;
   double* rez3 = rez + 3 * n;
   
   //Derivatives of Gaussian Function
   //f = p0 + p1 * exp(-(xi - p2)^2 / p3^2) - yi
   //df/dp0 = 1
   //df/dp1 = exp(-(xi - p2)^2 / p3^2)
   //df/dp2 = 2 * p1 * exp(-(xi - p2)^2 / p3^2) * (xi[i] - p2)   / p3^2    
   //df/dp3 = 2 * p1 * exp(-(xi - p2)^2 / p3^2) * (xi[i] - p2)^2 / p3^3
   
   for (int i = 0; i < n; i++)
     {
	rez0[i] = 1;
	rez1[i] = exp(-(xi[i] - p2) * (xi[i] - p2) / p3_2);	
	rez2[i] = 2 * rez1[i] * p1 * (xi[i] - p2) / p3_2;
	rez3[i] = rez2[i] * (xi[i] - p2) / p3;
     }
}
//                                                                            
#include <iostream>
using namespace std;
void asym_gaussian_res(double* p, double* xi, double *yi, int n, double* rez)
{   
   
   double p0   = p[0];
   double p1   = p[1];
   double p2   = p[2];
   double p3   = p[3];
   
//   cout<<p0<<" "<<p1<<" "<<p2<<" "<<p3<<endl;
   double den_less =  2 * p2 * (1 - p3) * p2 * (1 - p3); //2 * (p2 * (1 - p3))^2
   double den_more =  2 * p2 * (1 + p3) * p2 * (1 + p3); //2 * (p2 * (1 + p3))^2 
   for (int i = 0; i < n; i++)
     {
	double den;
	if (xi[i] < p1)
	  {
	     den = den_less;
	  }
	else  // if xi[i] >= p1
	  {
	     den = den_more;
	  }
//	rez[i] = 1 + p0 * exp(-(xi[i] - p1) * (xi[i] - p1) / den) - yi[i];
	rez[i] = 1 - p0/100 * exp(-(xi[i] - p1) * (xi[i] - p1) / den) - yi[i];
     }
}
//                                                                            
void asym_gaussian_res_J(double* p, double* xi, double *yi, int n, double* rez)
{   
   
   double p0   = p[0];
   double p1   = p[1];
   double p2   = p[2];
   double p3   = p[3];

   // new numpy array is row-major!  (http://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html)
   double* rez0 = rez;
   double* rez1 = rez + n;
   double* rez2 = rez + 2 * n;
   double* rez3 = rez + 3 * n;
   
   //Derivatives of left part of Asymmetric gaussian function 
   //f = 1 - p0/100 * exp(-(xi - p1)^2 / 2 /(p2 * (1 - p3))^2)
   //den:= 2 * (p2 * (1 - p3))^2
   //df/dp0 = -1/100 * exp(-(xi - p1)^2 / den)
   //df/dp1 = -1/100 * p0 * exp(-(xi - p1)^2 / den) * 2 * (xi - p1)  / den
   //df/dp2 = -1/100 * p0 * exp(-(xi - p1)^2 / den) * 2 * (xi - p1)^2 / 2 /(1 - p3)^2 / p2^3 
   //df/dp3 = -1/100 * p0 * exp(-(xi - p1)^2 / den) * (-2) * (xi - p1)^2 / 2 / p2^2 / (1-p3)^3

   //Derivatives of right part of Asymmetric gaussian function 
   //f = 1 - p0/100 * exp(-(xi - p1)^2 / 2 /(p2 * (1 + p3))^2)
   //den:= 2 * (p2 * (1 + p3))^2
   //df/dp0 = -1/100 * exp(-(xi - p1)^2 / den)
   //df/dp1 = -1/100 * p0 * exp(-(xi - p1)^2 / den) * 2 * (xi - p1)  / den
   //df/dp2 = -1/100 * p0 * exp(-(xi - p1)^2 / den) * 2 * (xi - p1)^2 / 2 /(1 + p3)^2 / p2^3 
   //df/dp3 = -1/100 * p0 * exp(-(xi - p1)^2 / den) * (-2) * (xi - p1)^2 / 2 / p2^2 / (1+p3)^3

   double den_less = 2 * p2 * (1 - p3) * p2 * (1 - p3); //(p2 * (1 - p3))^2
   double den_more = 2 * p2 * (1 + p3) * p2 * (1 + p3); //(p2 * (1 + p3))^2 
   
   for (int i = 0; i < n; i++)
     {
	double den;
	double den_rez3;
	if (xi[i] < p1)
	  {
	     den      = den_less;
	     den_rez3 = -(1 - p3);  
	  }
	else  // if xi[i] >= p1
	  {
	     den      = den_more;
	     den_rez3 = 1 + p3; 
	  }

	rez0[i] = -0.01 * exp(-(xi[i] - p1) * (xi[i] - p1) / den);
	rez1[i] =  rez0[i] * 2 * p0 * (xi[i] - p1) / den;
	rez2[i] =  rez1[i] * (xi[i] - p1) / p2 ;
	rez3[i] =  rez1[i] * (xi[i] - p1) / den_rez3;
     }
}
//                                                                            
