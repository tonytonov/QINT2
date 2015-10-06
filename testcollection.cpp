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

InterceptedSet FI05_CubicPolynomial::intercept(const real *x)
{
    int d = this->getDimension();
    std::vector<real> v(x, x + d);
    real f = 1;
    for (const auto x : v) f *= (x * x * x + 0.75);
    InterceptedSet res {f, v};
    return res;
}

InterceptedSet FI06_Oscillatory::intercept(const real *x)
{
    int d = this->getDimension();
    std::vector<real> v(x, x + d);
    real f = 1;
    for (int j = 0; j < d; j++)
    {
        f *= (j + 1) * (cos((j + 1) * v[j]) / sin(j + 1));
    }
    InterceptedSet res {f, v};
    return res;
}

InterceptedSet FI07_Singular::intercept(const real *x)
{
    int d = this->getDimension();
    std::vector<real> v(x, x + d);
    real f = 1;
    for (int j = 0; j < d; j++)
    {
        f *= 1.0 / (j + 1) * (pow(v[j], 1.0 / (j + 1) - 1));
    }
    InterceptedSet res {f, v};
    return res;
}

InterceptedSet FI08_InfVariation::intercept(const real *x)
{
    int d = this->getDimension();
    std::vector<real> v(x, x + d);
    real rad = pow(tgamma((double) d / 2 + 1) / (2 * pow(M_PI, (double) d / 2)), 1.0 / d);

    real r = 0;
    for (int j = 0; j < d; j++)
    {
        r += (v[j] - 0.5) * (v[j] - 0.5);
    }
    r = sqrt(r);

    real f = 0;
    if (r <= rad) f = 2; else f = 0;
    InterceptedSet res {f, v};
    return res;
}
