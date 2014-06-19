#ifndef UTILS_H
#define UTILS_H

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
T var(const std::vector<T> v)
{
    T m = mean(v);
    return (sum(square(v)) - v.size() * m * m) / (v.size() - 1);
}

template<typename F>
void repeat(unsigned n, F f) {
    while (n--) f();
}

#endif // UTILS_H
