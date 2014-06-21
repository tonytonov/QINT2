#include <vector>
#include <utils.h>
#include <testcollection.h>

real F00_simpleSum::operator()(const real *x)
{
    std::vector<real> v(x, x + this->getDimension());
    return sum(v);
}

InterceptedSet FI00_simpleSum::intercept(const real *x)
{
    std::vector<real> v(x, x + this->getDimension());
    InterceptedSet res {sum(v), v};
    return res;
}

InterceptedSet FI01_fMorCaf::intercept(const real *x)
{
    int d = this->getDimension();
    std::vector<real> v(x, x + d);
    real f = 1;
    for (const auto x : v) f *= (1 + 1.0 / d) * std::pow(x, 1.0 / d);
    InterceptedSet res {f, v};
    return res;
}
