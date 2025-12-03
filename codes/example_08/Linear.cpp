#include "Linear.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;

template <typename T, int P>
LinearSpectral2D<T, P>::LinearSpectral2D(
    basix::FiniteElement<T> element,
    mesh::Mesh<T> mesh) {};

NB_MODULE(Linear, m) {
    nb::class_<LinearSpectral2D<float, 4>>(m, "LinearSpectral2D")
        .def(nb::init<
            basix::FiniteElement<float>,
            mesh::Mesh<float>>());
}