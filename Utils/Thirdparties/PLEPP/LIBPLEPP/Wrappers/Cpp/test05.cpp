#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <boost/numeric/odeint.hpp>
/*
  + 
    http://www.boost.org/doc/libs/1_60_0/libs/numeric/odeint/doc/html/index.html

  + Driven harmonic oscillators
    https://en.wikipedia.org/wiki/Harmonic_oscillator
    http://math.colgate.edu/~wweckesser/math308Fall02/handouts/ForcedHarmonicOsc.pdf
*/


typedef std::vector< double > state_type;

class harm_osc 
{
  double _k, _c, _m, _F; 

  public:
  harm_osc( double k, double c, double m, double F=0.0) { _k=k; _c=c; _m=m; _F=F; }

  void operator() ( const state_type &U , state_type &dUdt , const double /* t */ )
  {                                                   //    U = (x,v);  m a = F - k x - c v -> dvdt = F/m - k/m x - c/m v; dxdt = v    
    dUdt[0] =           U[1];                         // dUdt =       v  
    dUdt[1] =  -_c/_m * U[1] - _k/_m * U[0] + _F/_m;  // dUdt =  -k/m x - c/m v + F  
  }
};

const double k = 5.790;    // N/m2 
const double c = 0.325e-3; // g/s2 -> 1kg/1e3g 
const double m = 2.979e-3; // g -> kg 
const double D = 0.160e-2; // cm -> m 


int main(int /* argc */ , char** /* argv */ )
{
  using namespace std;
  using namespace boost::numeric::odeint;

  state_type x(2);
  x[0] = D; 
  x[1] = 0.0; 

  runge_kutta4< state_type > stepper;

  const double dt = 0.25e-3;

  std::ofstream  fout("test05.dat");

  for( double t=0.0, j=0; t<10.0 ; t+= dt, j+=1 )
  {
    stepper.do_step( harm_osc(k,c,m), x , t , dt );
  //stepper.do_step( harm_osc(10.0,2.0,1.0,cos(2*t)), x , t , dt );  // test!! x=-0.5, v=4.0  
  
    fout<< j <<" "<< t <<" ";
    for(int i=0; i<x.size(); i++) fout<< x[i] <<" "; 
    fout<<"\n"; 
  }
  fout.close();


  return 0; 
} 
