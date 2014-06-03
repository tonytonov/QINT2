#ifndef RQMCINTEGRATOR_H
#define RQMCINTEGRATOR_H

#include <HIntLib/qmcintegrator.h>
#include <HIntLib/integrand.h>
using namespace HIntLib;

class RQMCIntegrator : public QMCIntegrator
{
public:
    RQMCIntegrator(PointSet* _ps) : QMCIntegrator(_ps) {}

    virtual
    Status integrate (
          Integrand &, const Hypercube &, Index maxEval,
          real reqAbsError, real reqRelError, EstErr &ee);
};

#endif // RQMCINTEGRATOR_H
