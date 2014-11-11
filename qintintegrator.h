#ifndef QINTINTEGRATOR_H
#define QINTINTEGRATOR_H

#include <HIntLib/qmcintegrator.h>
#include <HIntLib/integrand.h>
#include <HIntLib/pointset.h>
#include <HIntLib/defaults.h>
#include <HIntLib/hypercube.h>
#include <HIntLib/digitalnet2.h>
#include <interceptableintegrand.h>
#include <vector>

using namespace HIntLib;
typedef std::vector<std::vector<double>> t_sequence;

class QintIntegrator : public Integrator
{
public:
    QintIntegrator(Digital2PointSet* _ps,
                   unsigned _rn = 10,
                   unsigned _sp = 1,
                   unsigned _gs = 1,
                   unsigned _vo = 1)
        : ps(_ps), randCount(_rn), sParam(_sp), globalSeed(_gs), varOption(_vo) {}

    virtual
    Status integrate(
            Integrand &, const Hypercube &, Index maxEval,
            real reqAbsError, real reqRelError, EstErr &ee);
private:
    Digital2PointSet* ps;
    unsigned randCount;
    unsigned sParam;
    unsigned globalSeed;
    unsigned varOption;
    double estimateQintVariance(std::vector<double> &values, std::vector<int> &index);
};

class QINTPartitionIndexer
{
public:
    QINTPartitionIndexer(unsigned _p)
        : sParam(_p), partitionCount(1 << _p) {}
    virtual
    std::vector<int> CreateIndex(const t_sequence &sequence) = 0;
protected:
    unsigned sParam;
    unsigned partitionCount;
};

class CubicShapeIndexer : public QINTPartitionIndexer
{
public:
    CubicShapeIndexer(unsigned _p) : QINTPartitionIndexer(_p) {}
    virtual
    std::vector<int> CreateIndex(const t_sequence &sequence);
};

#endif // QINTINTEGRATOR_H
