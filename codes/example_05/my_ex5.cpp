#include <iostream>

#include <nanobind/nanobind.h>

namespace nb = nanobind;

template <typename T>
T process(T t) { return 2 * t; };

NB_MODULE(my_ex5, m) {
    m.def("process", &process<int>);
    m.def("process", &process<double>);
}