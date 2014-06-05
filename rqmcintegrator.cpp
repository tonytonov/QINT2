#include "rqmcintegrator.h"

Integrator::Status RQMCIntegrator::integrate(
        Integrand &f,
        const Hypercube &h,
        Index n, real, real,
        EstErr &ee)
{
    checkDimension(h, f);

    if (n == 0)
    {
       ee.set (0.0, 0.0);
       return ERROR;
    }

    ps->setCube (&h);
    ps->randomize (getSeed());
    n = ps->getOptimalNumber (n, h);
    Statistic<> stat;
    Point point (h.getDimension());
    ps->integrate(point, f, n, stat);
    ee.setNoErr (stat.getMean() * h.getVolume());
    return MAX_EVAL_REACHED;
}
