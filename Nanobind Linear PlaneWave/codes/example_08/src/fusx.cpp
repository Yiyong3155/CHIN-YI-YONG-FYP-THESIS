#include "Linear.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>

namespace nb = nanobind;

NB_MODULE(_fusx, m) {
    nb::class_<LinearSpectral2D<double, 4>>(m, "LinearSpectral2D")
        .def(nb::init<
            basix::FiniteElement<double>,
            std::shared_ptr<mesh::Mesh<double>>,
            std::shared_ptr<mesh::MeshTags<std::int32_t>>,
            std::shared_ptr<fem::Function<double>>,
            std::shared_ptr<fem::Function<double>>,
            const double&,
            const double&,
            const double&>(),
            nb::call_guard<nb::gil_scoped_release>()  // GIL released during constructor
        )
        .def("init",
             &LinearSpectral2D<double, 4>::init,
             nb::call_guard<nb::gil_scoped_release>()     // release GIL
        )
        .def("rk4",
             &LinearSpectral2D<double, 4>::rk4,
             nb::call_guard<nb::gil_scoped_release>()     // release GIL
        )
        .def("u_sol",
             &LinearSpectral2D<double, 4>::u_sol,
             nb::call_guard<nb::gil_scoped_release>()     // release GIL
        );
}
