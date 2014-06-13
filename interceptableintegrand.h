#ifndef INTERCEPTABLEINTEGRAND_H
#define INTERCEPTABLEINTEGRAND_H

#include <HIntLib/defaults.h>
#include <HIntLib/integrand.h>
#include <vector>

using namespace HIntLib;

struct InterceptedValue {
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
        return res.value;
    }
    virtual InterceptedValue intercept (const real []) = 0;

private:
    std::vector<std::vector<real>> interceptedPoints;
public:
    std::vector<std::vector<real>> getInterceptedPoints() { return interceptedPoints; }
    void eraseInterceptedPoints() { interceptedPoints.clear(); }
    void reserveInterceptedPoints(int n) { interceptedPoints.reserve(n); }
};

#endif // INTERCEPTABLEINTEGRAND_H
