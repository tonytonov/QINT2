#include "qintintegrator.h"
#include "utils.h"
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

    std::vector<t_sequence> interceptedSequences;
    interceptedSequences.reserve(randCount);
    auto iif = dynamic_cast<InterceptableIntegrand*>(&f);
    repeat(randCount, [&]
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
            interceptedSequences.push_back(iif->getInterceptedPoints());
        }
        stats.push_back(s);
    });

    std::vector<double> estimates;
    estimates.reserve(randCount);
    std::vector<double> variances;
    variances.reserve(randCount);
    for (auto x : stats) estimates.push_back(x.getMean() * h.getVolume());
    double qintEst = sum(estimates) / randCount;
    CubicShapeIndexer indexer(2);
    std::vector<std::vector<int>> indexes;
    indexes.reserve(randCount);
    for (auto x : interceptedSequences)
    {
        indexes.push_back(indexer.CreateIndex(x));
    }
    //variances = getQintVariances(interceptedSequence, indexSequence);
    //double qintVar = sum(variances) / randCount / randCount;
    //double qintStdErr = std::sqrt(qintVar / randCount);
    ee.set(qintEst, 0);
    return MAX_EVAL_REACHED;
}


std::vector<int> CubicShapeIndexer::CreateIndex(t_sequence sequence)
{
    std::vector<int> res;
    res.reserve(sequence.size());
    for (const auto v : sequence)
    {
        int d = v.size();
        std::vector<int> partTimes;
        partTimes.reserve(d);
        std::vector<int> binaryIndex;
        binaryIndex.reserve(d);
        int index = 0;
        auto dv = std::div(sParam, d);
        int a = dv.quot;
        int b = dv.rem;
        for (int i = 0; i < d; i++)
        {
            (i < b) ? partTimes.push_back(a + 1) : partTimes.push_back(a);
            binaryIndex.push_back(std::floor(v[i] * std::pow(2, partTimes[i])));
            index += binaryIndex[i] * std::pow(2, sParam - sum(partTimes));
        }
        res.push_back(index);
    }
    return res;
}
