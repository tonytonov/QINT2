#ifndef RQMCINTEGRATOR_H
#define RQMCINTEGRATOR_H

#include <HIntLib/qmcintegrator.h>
#include <HIntLib/integrand.h>
#include <HIntLib/pointset.h>
#include <HIntLib/defaults.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/digitalnet2.h>

using namespace HIntLib;

class RQMCIntegrator : public Integrator
{
public:
    RQMCIntegrator(Digital2PointSet* _ps, unsigned int _rn = 10, unsigned int _gs = 1)
        : ps(_ps), randCount(_rn), globalSeed(_gs) {}

    virtual
    Status integrate (
          Integrand &, const Hypercube &, Index maxEval,
          real reqAbsError, real reqRelError, EstErr &ee);
private:
    Digital2PointSet* ps;
    unsigned int randCount;
    unsigned int globalSeed;
};

#endif // RQMCINTEGRATOR_H
