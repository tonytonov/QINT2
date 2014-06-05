#include "rqmcintegrator.h"
#include <numeric>

double sum(std::vector<double> v)
{
    double sum = 0;
    for (double n : v) sum += n;
    return sum;
}

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
    n = ps->getOptimalNumber (n, h);

    if (randNum == 0) return ERROR;
    int m = n / randNum;

    std::vector<Statistic<>> stats;
    stats.reserve(m);
    for (unsigned int i = 0; i < randNum; i++)
    {
        Statistic<> s;
        Point point (h.getDimension());
        ps->randomize (getSeed());
        ps->integrate(point, f, m, s);
        stats.push_back(s);
    }

    std::vector<double> means(m);
    //std::vector<double> stds(m);
    for (auto x : stats) means.push_back(x.getMean());
    ee.setNoErr(sum(means) / randNum * h.getVolume());
    return MAX_EVAL_REACHED;
}
