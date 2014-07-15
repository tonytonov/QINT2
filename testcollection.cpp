#include <vector>
#include <random>
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

InterceptedSet FI01_MorCaf1::intercept(const real *x)
{
    int d = this->getDimension();
    std::vector<real> v(x, x + d);
    real f = 1;
    for (const auto x : v) f *= (1 + 1.0 / d) * std::pow(x, 1.0 / d);
    InterceptedSet res {f, v};
    return res;
}

InterceptedSet FI02_NormalDensity::intercept(const real *x)
{
    real mean = 0.5;
    real mult = 0.3829249;
    int d = this->getDimension();
    std::vector<real> v(x, x + d);
    real f = 1;
    for (const auto x : v) f *= dnorm(x - mean) / mult;
    InterceptedSet res {f, v};
    return res;
}

InterceptedSet FI03_PieceLin::intercept(const real *x)
{
    int d = this->getDimension();
    std::vector<real> v(x, x + d);
    real value = 0.0;
    real total = 0.0;
    for (int j = 0; j < d; j++)
    {
        if (v[j] <= 0.5 - double(j) / 2 / (j + 10)) {total = 0.0;}
        else if (v[j] >= 0.5 + double(j) / 2 / (j + 10)) {total = 1.0;}
        else {total = v[j] * (10 + (double)j) / j - 5.0 / j;}
        value += total * 2;
    }
    InterceptedSet res {value, v};
    return res;
}

InterceptedSet FI04_PieceLinEx::intercept(const real *x)
{
    int d = this->getDimension();
    std::vector<real> v(x, x + d);
    real value = 0.0;
    real total = 1.0;
    for (int j = 0; j < d; j++)
    {
        if (v[j] <= 0.5 - double(j) / 2 / (j + 10)) {total = 0.0;}
        else if (v[j] >= 0.5 + double(j) / 2 / (j + 10)) {total = 1.0;}
        else {total = v[j] * (10 + (double)j) / j - 5.0 / j;}
        value *= total * 2;
    }
    InterceptedSet res {value, v};
    return res;
}
