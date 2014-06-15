#include "qintintegrator.h"
#include <vector>
#include <cmath>

Integrator::Status QINTIntegrator::integrate(
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

    std::vector<std::vector<std::vector<double>>> interceptedSequence;
    interceptedSequence.reserve(randCount);
    auto iif = dynamic_cast<InterceptableIntegrand*>(&f);
    for (unsigned int i = 0; i < randCount; i++)
    {
        Statistic<> s;
        Point point(h.getDimension());
        int seed = e();
        ps->randomize(seed);
        if (iif)
        {
            iif->eraseIntercepted();
            iif->reserveIntercepted(m);
        }
        ps->integrate(point, f, m, s);
        if (iif)
        {
            interceptedSequence.push_back(iif->getInterceptedPoints());
        }
        stats[i] = s;
    }

    std::vector<double> estimates;
    estimates.reserve(randCount);
    std::vector<double> variances;
    variances.reserve(randCount);
    for (auto x : stats) estimates.push_back(x.getMean() * h.getVolume());
    double qintEst = sum(estimates) / randCount;
    //variances = getQintVariances(interceptedSequence, indexSequence);
    //double qintVar = sum(variances) / randCount / randCount;
    //double qintStdErr = std::sqrt(qintVar / randCount);
    ee.set(qintEst, 0);
    return MAX_EVAL_REACHED;
}
