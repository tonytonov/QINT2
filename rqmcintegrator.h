#ifndef RQMCINTEGRATOR_H
#define RQMCINTEGRATOR_H

#include <HIntLib/qmcintegrator.h>
#include <HIntLib/integrand.h>
#include <HIntLib/pointset.h>
#include <HIntLib/defaults.h>
#include <HIntLib/hypercube.h>

using namespace HIntLib;

class RQMCIntegrator : public Integrator
{
public:
    RQMCIntegrator(PointSet* _ps, unsigned int _rn = 10) : ps(_ps), randNum(_rn) {}

    virtual
    Status integrate (
          Integrand &, const Hypercube &, Index maxEval,
          real reqAbsError, real reqRelError, EstErr &ee);
private:
    PointSet* ps;
    unsigned int randNum;
};

#endif // RQMCINTEGRATOR_H
