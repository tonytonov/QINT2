#ifndef TESTCOLLECTION_H
#define TESTCOLLECTION_H

#include <vector>
#include <HIntLib/qmcintegrator.h>
#include <HIntLib/integrand.h>
#include <HIntLib/pointset.h>
#include <HIntLib/defaults.h>
#include <interceptableintegrand.h>

class F00_simpleSum : public Integrand
{
public:
    F00_simpleSum(int s) : Integrand(s) {}
    virtual real operator() (const real *);
};

class FI00_simpleSum : public InterceptableIntegrand
{
public:
    FI00_simpleSum(int s) : InterceptableIntegrand(s, (double) s / 2)
    { name = "Simple sum"; }
    virtual InterceptedSet intercept(const real *);
};

class FI01_MorCaf1 : public InterceptableIntegrand
{
public:
    FI01_MorCaf1(int s) : InterceptableIntegrand(s, 1.0)
    { name = "Morokoff-Caflisch function #1";}
    virtual InterceptedSet intercept(const real *);
};

class FI03_PieceLin : public InterceptableIntegrand
{
public:
    FI03_PieceLin(int s) : InterceptableIntegrand(s, s)
    { name = "Piecewise linear function";}
    virtual InterceptedSet intercept(const real *);
};

#endif // TESTCOLLECTION_H
