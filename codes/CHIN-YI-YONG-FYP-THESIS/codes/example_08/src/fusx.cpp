#include "Linear.hpp"

#include <cstddef>
#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/array.h>

using namespace dolfinx;

namespace nb = nanobind;

NB_MODULE(_fusx, m) {
    // ===== 2D Linear Spectral Solver =====
    nb::class_<LinearSpectral2D<double, 4>>(m, "LinearSpectral2D")
        .def(nb::init<  // ✅ ADD < HERE!
            basix::FiniteElement<double>,
            std::shared_ptr<mesh::Mesh<double>>,
            std::shared_ptr<mesh::MeshTags<std::int32_t>>,
            std::shared_ptr<fem::Function<double>>,
            std::shared_ptr<fem::Function<double>>,
            const double&,
            const double&,
            const double&,
            const std::array<double, 2>&,
            const std::array<double, 2>&,
            std::size_t,
            std::size_t>(),
            nb::arg("element"),
            nb::arg("mesh"),
            nb::arg("facet_tags"),
            nb::arg("speed_of_sound"),
            nb::arg("density"),
            nb::arg("source_frequency"),
            nb::arg("source_amplitude"),
            nb::arg("source_speed"),
            nb::arg("region_min") = std::array<double, 2>{-2.25, -2.25},
            nb::arg("region_max") = std::array<double, 2>{2.25, 2.25},
            nb::arg("grid_nx") = std::size_t(251),
            nb::arg("grid_ny") = std::size_t(251),
            nb::call_guard<nb::gil_scoped_release>()
        )
        .def("init",
             &LinearSpectral2D<double, 4>::init,
             nb::call_guard<nb::gil_scoped_release>())
        .def("rk4",
             &LinearSpectral2D<double, 4>::rk4,
             nb::call_guard<nb::gil_scoped_release>())
        .def("u_sol",
             &LinearSpectral2D<double, 4>::u_sol,
             nb::call_guard<nb::gil_scoped_release>());

    // ===== 3D Linear Spectral Solver =====
    nb::class_<LinearSpectral3D<double, 4>>(m, "LinearSpectral3D")
        .def(nb::init<  // ✅ ADD < HERE!
            basix::FiniteElement<double>,
            std::shared_ptr<mesh::Mesh<double>>,
            std::shared_ptr<mesh::MeshTags<std::int32_t>>,
            std::shared_ptr<fem::Function<double>>,
            std::shared_ptr<fem::Function<double>>,
            const double&,
            const double&,
            const double&,
            const std::array<double, 3>&,
            const std::array<double, 3>&,
            std::size_t,
            std::size_t,
            std::size_t>(),
            nb::arg("element"),
            nb::arg("mesh"),
            nb::arg("facet_tags"),
            nb::arg("speed_of_sound"),
            nb::arg("density"),
            nb::arg("source_frequency"),
            nb::arg("source_amplitude"),
            nb::arg("source_speed"),
            nb::arg("region_min") = std::array<double, 3>{-2.25, -2.25, -2.25},
            nb::arg("region_max") = std::array<double, 3>{2.25, 2.25, 2.25},
            nb::arg("grid_nx") = std::size_t(101),
            nb::arg("grid_ny") = std::size_t(101),
            nb::arg("grid_nz") = std::size_t(101),
            nb::call_guard<nb::gil_scoped_release>()
        )
        .def("init",
             &LinearSpectral3D<double, 4>::init,
             nb::call_guard<nb::gil_scoped_release>())
        .def("rk4",
             &LinearSpectral3D<double, 4>::rk4,
             nb::call_guard<nb::gil_scoped_release>())
        .def("u_sol",
             &LinearSpectral3D<double, 4>::u_sol,
             nb::call_guard<nb::gil_scoped_release>())
        .def("number_of_dofs",
             &LinearSpectral3D<double, 4>::number_of_dofs);
}