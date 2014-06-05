#include "rqmcintegrator.h"
#include <numeric>
#include <cmath>

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
       ee.set(0.0, 0.0);
       return ERROR;
    }

    ps->setCube(&h);
    n = ps->getOptimalNumber(n, h);
    ps->enableRandomize();

    if (randCount == 0) return ERROR;
    int m = n / randCount;

    std::vector<Statistic<>> stats(randCount);
    std::default_random_engine e(globalSeed);

    for (unsigned int i = 0; i < randCount; i++)
    {
        Statistic<> s;
        Point point(h.getDimension());
        int seed = e();
        ps->randomize(seed);
        ps->integrate(point, f, m, s);
        stats[i] = s;
    }

    std::vector<double> estimates;
    estimates.reserve(randCount);
    std::vector<double> sqdiffs;
    sqdiffs.reserve(randCount);
    for (auto x : stats) estimates.push_back(x.getMean() * h.getVolume());
    double rqmcEst = sum(estimates) / randCount;
    for (auto x : estimates) sqdiffs.push_back(std::pow(rqmcEst - x, 2));
    double rqmcStdError = std::sqrt(sum(sqdiffs) / randCount / (randCount - 1));
    ee.set(rqmcEst, rqmcStdError);
    return MAX_EVAL_REACHED;
}
