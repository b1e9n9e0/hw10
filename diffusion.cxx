#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
//---------------------------------------
using namespace std;
//---------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N, const double t);
void initialize(double* const u0, const double dx,const double dt, const double xmin,
                const int N);
void step(double* const u1,  double* const u0,  const double dt,
          const double dx, const double D,const int N);

//---------------------------------------
int main(){

  const double tEnd = 5 ;
  const double D = 1;

  const int N  = 200;
  const double xmin = -20;
  const double xmax = 20;
  const double dx = (xmax-xmin)/(N-1) ;

  double dt = dx;
  double t = 0;
  const int Na = 10;
  const int Nk = int(tEnd/Na/dt);


  double* u0 = new double[N];
  double* u1 = new double[N];
  double* h;
  stringstream strm;

  initialize(u0,dx,dt, xmin,N);

  writeToFile(u0, "u_0", dx, xmin, N,t);

  cout << "Nk = " << Nk << endl;

  for(int i=1; i<=Na; i++)
  {
   for(int j=0; j<Nk; j++){

   	step(u1,u0,dt,dx,D,N);

      h = u0;
      u0 = u1;
      u1 = h;
      t +=dt;

   }
   strm.str("");
   strm << "u_" << i;
   writeToFile(u0, strm.str(), dx, xmin, N,t);
  }

  cout << "t = " << t << endl;

  delete[] u0;
  delete[] u1;
  return 0;
}
//-----------------------------------------------
void step(double* const f1, double* const f0,
          const double dt, const double dx,
          const double D, const int N)
{
	double* d=new double[N];
  double* u=new double[N];
  double* l=new double[N];
  double C = D*dt/(dx*dx);
  double x;

  for(int i=0;i<N;i++) d[i] = -C*(2.0+(dx*dx));
  for(int i=0;i<N;i++) u[i] = C;
  for(int i=0;i<N;i++) l[i] = C;
  	//obere diagonalmatrix erstellen
  for(int i = 1; i < N; i++){
  	x = l[i]/d[i-1];
  	d[i] -= x*u[i-1];
  	//l[i] = 0; auch ohne rechenne ist  das jetzt so...
  	f0[i] -= x*f0[i-1];
  }
  for(int i = N-2; i >= 0; i--){
  	f1[i] = f0[i] + D*dt/(dx*dx)*(f0[(i+1)%N] - 2*f0[i] + f0[(N+i-1)%N]);
  }

  delete[] d;
  delete[] u;
  delete[] l;
}
//-----------------------------------------------
void initialize(double* const u0, const double dx,
                const double dt, const double xmin,  const int N)
{
   double u,ux, uxx;
   for(int i=0; i<N; i++)
   {
     double x = xmin + i*dx;
     u0[i] = 1.0/sqrt(4*M_PI)*exp(-x*x/4.0);

   }
}
//-----------------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N, const double t)
{
   ofstream out(s.c_str());
   for(int i=0; i<N; i++){
     double x = xmin + i * dx;
     double ana = 1.0/sqrt(4*M_PI*(t+1.0)) * exp(-x*x/(4*(t+1)));
     out << x << "\t" << u[i] << "\t" << ana << endl;
   }
   out.close();
}
