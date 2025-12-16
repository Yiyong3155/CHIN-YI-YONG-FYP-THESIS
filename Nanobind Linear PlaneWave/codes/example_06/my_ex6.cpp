#include "my_ex6.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;

template <typename T>
void print(T n) {
    for (int i = 0; i < 5; ++i) {
        std::cout << i * n << "\n";
    }
}

template <typename T>
Compute<T>::Compute(T a_) : a(a_) { };

template <typename T>
void Compute<T>::run() {
    print<T>(a);
};

NB_MODULE(my_ex6, m) {
    nb::class_<Compute<int>>(m, "Compute")
        .def(nb::init<int>())
        .def("run", &Compute<int>::run);
    nb::class_<Compute<double>>(m, "Compute")
        .def(nb::init<double>())
        .def("run", &Compute<double>::run);
}
