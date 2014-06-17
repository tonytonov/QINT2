#ifndef UTILS_H
#define UTILS_H

template<typename T>
T sum(std::vector<T> v)
{
    double sum = 0;
    for (auto n : v) sum += n;
    return sum;
}

template<typename F>
void repeat(unsigned n, F f) {
    while (n--) f();
}

#endif // UTILS_H
