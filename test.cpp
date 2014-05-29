#include <HIntLib/integrand.h>
#include <HIntLib/integrator.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/esterr.h>
#include <HIntLib/digitalnet2.h>
#include <HIntLib/sobolmatrix.h>
#include <HIntLib/niederreitermatrix.h>
#include <iostream>

int main()
{

  int s = 5;
  
  class TestIntegrand: public Integrand
  {
    public:
      real operator() (const real [] x)
      {
        
      }
  }

  TestIntegrand f(s);
  
  MonteCarloPointSet<MersenneTwister> ps;
  MCIntegrator integrator_mc (&ps);
  integrator_mc.setMinNumEval (300);
  
  NiederreiterMatrix matrix;
  DigitalNet2PointSet<real> ps (matrix, true, DigitalNet::CENTER);
  QMCIntegrator integrator_nied (&ps);
  
  SobolMatrix matrix;
  DigitalSeq2PointSet<real> ps (matrix, true);
  QMCIntegrator integrator_sobol (&ps);
  
  void calculateIntegral (Integrand& f, Integrator& integrator)
  {
    Hypercube h (f.getDimension());
    cout << "Result: "
    << integrator(f, h, 1000000, 0, 1/1000)
    << "\n";
  }
  
  calculateIntegral(&f, &integrator_mc)
  calculateIntegral(&f, &integrator_sobol)
  calculateIntegral(&f, &integrator_nied)
}
