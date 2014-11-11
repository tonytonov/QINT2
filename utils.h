#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <cmath>

template<typename T>
std::string S(T t)
{
    return std::to_string(t);
}

template<typename T>
T sum(const std::vector<T> v)
{
    double sum = 0;
    for (auto n : v) sum += n;
    return sum;
}

template<typename T>
T mean(const std::vector<T> v)
{
    return sum(v) / v.size();
}

template<typename T>
std::vector<T> square(const std::vector<T> v)
{
    std::vector<T> res;
    res.reserve(v.size());
    for (const auto x : v) res.push_back(x*x);
    return res;
}

template<typename T>
T sumsq(const std::vector<T> v)
{
    return (sum(square(v)));
}

template<typename T>
T var(const std::vector<T> v)
{
    T m = mean(v);
    return (sum(square(v)) - v.size() * m * m) / (v.size() - 1);
}

template<typename F>
void repeat(unsigned n, F f) {
    while (n--) f();
}

inline double pnorm(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = std::fabs(x)/std::sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*std::exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

inline double dnorm(double x, double m = 0, double s = 1)
{
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}
#endif // UTILS_H
