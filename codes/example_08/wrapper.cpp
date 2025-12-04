#include "fusx/Linear.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>

namespace nb = nanobind;

NB_MODULE(Linear, m) {
    nb::class_<LinearSpectral2D<double, 4>>(m, "LinearSpectral2D")
        .def(nb::init<
            basix::FiniteElement<double>,
            std::shared_ptr<mesh::Mesh<double>>,
            std::shared_ptr<mesh::MeshTags<std::int32_t>>,
            std::shared_ptr<fem::Function<double>>,
            std::shared_ptr<fem::Function<double>>,
            const double&,
            const double&,
            const double&>());
}