#ifndef INTERCEPTABLEINTEGRAND_H
#define INTERCEPTABLEINTEGRAND_H

#include <HIntLib/defaults.h>
#include <HIntLib/integrand.h>
#include <vector>
#include <limits>
#include <cmath>

using namespace HIntLib;

struct InterceptedSet {
    real value;
    std::vector<real> vec;
};

class InterceptableIntegrand : public Integrand
{
public:
    InterceptableIntegrand (int s, real exactValue = std::numeric_limits<double>::quiet_NaN()): Integrand(s), exactValue(exactValue) {}
    virtual ~InterceptableIntegrand() {}
    virtual real operator() (const real* x) {
        auto res = intercept(x);
        interceptedPoints.push_back(res.vec);
        interceptedValues.push_back(res.value);
        return res.value;
    }
    virtual InterceptedSet intercept (const real []) = 0;

private:
    std::vector<std::vector<real>> interceptedPoints;
    std::vector<real> interceptedValues;
    real exactValue;
public:
    std::vector<std::vector<real>> getInterceptedPoints() { return interceptedPoints; }
    std::vector<real> getInterceptedValues() { return interceptedValues; }
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
    real getExactValue() { return exactValue; }
};

#endif // INTERCEPTABLEINTEGRAND_H
