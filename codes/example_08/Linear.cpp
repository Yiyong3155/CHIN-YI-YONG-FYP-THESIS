#include "Linear.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;

template <typename T>
LinearSpectral2D<T>::LinearSpectral2D(basix::FiniteElement<T> element) {};

NB_MODULE(Linear, m) {
    nb::class_<LinearSpectral2D<float>>(m, "LinearSpectral2D")
        .def(nb::init<basix::FiniteElement<float>>());
}