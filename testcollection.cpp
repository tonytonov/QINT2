#include <vector>
#include <testcollection.h>

real TestFunction::operator()(const real *x)
{
    std::vector<real> v(x, x + this->getDimension());
    real sum = 0;
    for (auto n : v) sum += n;
    return sum;
}

InterceptedSet TestFunctionInterceptable::intercept(const real *x)
{
    std::vector<real> v(x, x + this->getDimension());
    real sum = 0;
    for (auto n : v) sum += n;
    InterceptedSet res {sum, v};
    return res;
}
