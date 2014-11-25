#include <rqmcintegrator.h>
#include <utils.h>
#include <vector>
#include <cmath>

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
    for (auto x : stats) estimates.push_back(x.getMean() * h.getVolume());
    double rqmcEst = sum(estimates) / randCount;
    double rqmcStdError = std::sqrt(var(estimates) / randCount);
    if (randCount == 1) rqmcStdError = -1;
    ee.set(rqmcEst, rqmcStdError);
    return MAX_EVAL_REACHED;
}
