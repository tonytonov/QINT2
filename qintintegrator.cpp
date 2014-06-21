#include <qintintegrator.h>
#include <utils.h>
#include <vector>
#include <cmath>

Integrator::Status QintIntegrator::integrate(
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

    std::vector<Statistic<>> stats;
    stats.reserve(randCount);
    std::default_random_engine e(globalSeed);

    std::vector<t_sequence> interceptedSequences;
    interceptedSequences.reserve(randCount);
    std::vector<std::vector<real>> interceptedValues;
    interceptedValues.reserve(randCount);
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
            interceptedValues.push_back(iif->getInterceptedValues());
        }
        stats.push_back(s);
    });

    std::vector<double> estimates;
    estimates.reserve(randCount);
    std::vector<double> variances;
    variances.reserve(randCount);
    std::vector<double> mixed_vals;
    mixed_vals.reserve(n);
    for (const auto x : stats) estimates.push_back(x.getMean() * h.getVolume());
    double qintEst = sum(estimates) / randCount;
    CubicShapeIndexer indexer(sParam);
    std::vector<std::vector<int>> indexes;
    indexes.reserve(randCount);
    unsigned i=0;
    for (const auto x : interceptedSequences)
    {
        auto currentIndex = indexer.CreateIndex(x);
        variances.push_back(estimateQintVariance(interceptedValues[i], currentIndex));
        mixed_vals.insert(mixed_vals.end(), interceptedValues[i].begin(), interceptedValues[i].end());
        ++i;
    }
    double rqmcStdError = std::sqrt(var(estimates) / randCount);
    double mcStdError = std::sqrt(var(mixed_vals) / n);
    double qintStdError = std::sqrt(sum(variances) / randCount / randCount);
    ee.set(qintEst, qintStdError);
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

double QintIntegrator::estimateQintVariance(std::vector<double> values, std::vector<int> index)
{
    if (values.size() != index.size()) throw ("Sequence and index sizes do not match!");
    double totalVar = var(values);
    double alphaTerm = 0;
    std::vector<double> alphas(std::pow(2, sParam));
    std::vector<unsigned> alphaCounters(std::pow(2, sParam));
    for (unsigned i=0; i<values.size(); i++)
    {
        alphas[index[i]] += values[i];
        ++alphaCounters[index[i]];
    }
    for (unsigned j=0; j<alphas.size(); j++)
    {
        alphas[j] = alphas[j] / alphaCounters[j] / std::pow(2, sParam);
    }
    for (unsigned i=0; i<alphas.size(); i++)
    {
        for (unsigned j=i; j<alphas.size(); j++)
        {
            alphaTerm += (alphas[i] - alphas[j]) * (alphas[i] - alphas[j]);
        }
    }
    return (totalVar - alphaTerm) / values.size();
}
