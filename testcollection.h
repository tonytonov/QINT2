#ifndef TESTCOLLECTION_H
#define TESTCOLLECTION_H

#include <vector>
#include <HIntLib/qmcintegrator.h>
#include <HIntLib/integrand.h>
#include <HIntLib/pointset.h>
#include <HIntLib/defaults.h>
#include <interceptableintegrand.h>

class TestFunction : public Integrand
{
public:
    TestFunction(int s) : Integrand(s) {}
    virtual ~TestFunction() {}
    virtual real operator() (const real *);
};

class TestFunctionInterceptable : public InterceptableIntegrand
{
public:
    TestFunctionInterceptable(int s) : InterceptableIntegrand(s) {}
    virtual ~TestFunctionInterceptable() {}
    virtual InterceptedSet intercept(const real *);
};

#endif // TESTCOLLECTION_H
