//fitgauss_c.h

extern "C" {
   
   //function rezi = p0 + p1 * exp(-(xi - p2)^2 / p3^2) - yi;
   void gaussian_res  (double* p, double* xi, double *yi, int n, double* rez);
   
   //jacobian for gaussian_res
   void gaussian_res_J(double* p, double* xi, double *yi, int n, double* rez);
   
   //if (xi < p1)  rezi = 1 - p0/100 * exp(-(xi - p1)^2 / 2 / (p2 * (1 - p3))^2) - yi
   //if (xi >= p1) rezi = 1 - p0/100 * exp(-(xi - p1)^2 / 2 / (p2 * (1 + p3))^2) - yi
   void asym_gaussian_res (double* p, double* xi, double *yi, int n, double* rez);

   //jacobian for asym_gaussian_res
   void asym_gaussian_res_J(double* p, double* xi, double *yi, int n, double* rez);
}

