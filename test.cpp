#include <HIntLib/integrand.h>
#include <HIntLib/integrator.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/esterr.h>
#include <HIntLib/digitalnet2.h>
#include <HIntLib/sobolmatrix.h>
#include <HIntLib/niederreitermatrix.h>
#include <iostream>

class TestIntegrand: public Integrand
{
    public real operator() (const real [] x)
    {
      return x;
    }
}
