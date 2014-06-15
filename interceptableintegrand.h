#ifndef INTERCEPTABLEINTEGRAND_H
#define INTERCEPTABLEINTEGRAND_H

#include <HIntLib/defaults.h>
#include <HIntLib/integrand.h>
#include <vector>

using namespace HIntLib;

struct InterceptedSet {
    real value;
    std::vector<real> vec;
};

class InterceptableIntegrand : public Integrand
{
public:
    InterceptableIntegrand (int s): Integrand(s) {}
    virtual ~InterceptableIntegrand() {}
    virtual real operator() (const real* x) {
        auto res = intercept(x);
        interceptedPoints.push_back(res.vec);
        interceptedValues.push_back(res.value);
        return res.value;
    }
    virtual InterceptedSet intercept (const real []) {
        //FIXME
    }

private:
    std::vector<std::vector<real>> interceptedPoints;
    std::vector<real> interceptedValues;
public:
    std::vector<std::vector<real>> getInterceptedPoints() { return interceptedPoints; }
    std::vector<std::vector<real>> getInterceptedValues() { return interceptedValues; }
    void eraseIntercepted()
    {
        interceptedPoints.clear();
        interceptedValues.clear();
    }
    void reserveIntercepted(int n)
    {
        interceptedPoints.reserve(n);
        interceptedValues.reserve(n);
    }
};

#endif // INTERCEPTABLEINTEGRAND_H
