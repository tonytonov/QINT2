#include "qintintegrator.h"
#include <vector>
#include <cmath>

Integrator::Status QINTIntegrator::integrateAndIntercept(
        InterceptableIntegrand &f,
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
        f.eraseInterceptedPoints();
        f.reserveInterceptedPoints(m);
        ps->integrate(point, f, m, s);
        auto inter = f.getInterceptedPoints();
        stats[i] = s;
    }

    std::vector<double> estimates;
    estimates.reserve(randCount);
    std::vector<double> variances;
    variances.reserve(randCount);
    for (auto x : stats) estimates.push_back(x.getMean() * h.getVolume());
    double qintEst = 0;//sum(estimates) / randCount;
    //variances =
    //double qintStdError = std::sqrt(sum(sqdiffs) / randCount / (randCount - 1));
    ee.set(qintEst, 0);
    return MAX_EVAL_REACHED;
}

Integrator::Status QINTIntegrator::integrate(
        Integrand &,
        const Hypercube &,
        Index, real, real,
        EstErr &)
{
    return ERROR;
}
