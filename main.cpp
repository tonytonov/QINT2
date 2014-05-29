#include <HIntLib/integrand.h>
#include <HIntLib/integrator.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/esterr.h>
#include <HIntLib/digitalnet2.h>
#include <HIntLib/sobolmatrix.h>
#include <HIntLib/niederreitermatrix.h>
#include <iostream>

void calculateIntegral (Integrand& f, Integrator& integrator)
{
  Hypercube h (f.getDimension());
  cout << "Result: "
  << integrator(f, h, 1000000, 0, 1/1000)
  << "\n";
}

int main()
{
    TestIntegrand f(s);

    MonteCarloPointSet<MersenneTwister> ps_mc;
    MCIntegrator integrator_mc (&ps_mc);
    integrator_mc.setMinNumEval (300);

    NiederreiterMatrix matrix_nied;
    DigitalNet2PointSet<real> ps_nied (matrix_nied, true, DigitalNet::CENTER);
    QMCIntegrator integrator_nied (&ps_nied);

    SobolMatrix matrix_sobol;
    DigitalSeq2PointSet<real> ps_sobol (matrix_sobol, true);
    QMCIntegrator integrator_sobol (&ps_sobol);

    calculateIntegral(&f, &integrator_mc);
    calculateIntegral(&f, &integrator_sobol);
    calculateIntegral(&f, &integrator_nied);

}
