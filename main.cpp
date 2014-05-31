#include <HIntLib/integrand.h>
#include <HIntLib/integrator.h>
#include <HIntLib/mcintegrator.h>
#include <HIntLib/mersennetwister.h>
#include <HIntLib/qmcintegrator.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/esterr.h>
#include <HIntLib/digitalnet2.h>
#include <HIntLib/sobolmatrix.h>
#include <HIntLib/niederreitermatrix.h>
#include <HIntLib/defaults.h>
#include <HIntLib/mcpointset.h>
#include <iostream>
#include <numeric>

using namespace HIntLib;
using namespace std;

class TestFunction : public Integrand
{
public:
    TestFunction(int s) : Integrand(s) {}
    virtual ~TestFunction() {}
    virtual real operator() (const real *x)
    {
      vector<real> v(x, x + this->getDimension());
      real sum = 0;
      for (real n : v) sum += n;
      return sum;
    }

};

void calculateIntegral (TestFunction& f, Integrator& integrator)
{
  Hypercube h(f.getDimension());
  std::cout << "Result: "
  << integrator(f, h, 1000000, 0, 1/1000)
  << "\n";
}

int main()
{
    int s=5;
    TestFunction f(s);

    MonteCarloPointSet<MersenneTwister> ps_mc;
    MCIntegrator integrator_mc (&ps_mc);
    integrator_mc.setMinNumEval (300);

    NiederreiterMatrix matrix_nied;
    DigitalNet2PointSet<real> ps_nied (matrix_nied, true, DigitalNet::CENTER);
    QMCIntegrator integrator_nied (&ps_nied);

    SobolMatrix matrix_sobol;
    DigitalSeq2PointSet<real> ps_sobol (matrix_sobol, true);
    QMCIntegrator integrator_sobol (&ps_sobol);

    calculateIntegral(f, integrator_mc);
    calculateIntegral(f, integrator_sobol);
    calculateIntegral(f, integrator_nied);

}
