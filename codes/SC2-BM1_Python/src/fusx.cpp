#include "Linear.hpp"

#include <cstddef>
#include <string>
#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>

using namespace dolfinx;
namespace nb = nanobind;

NB_MODULE(_fusx, m) {

    nb::class_<LinearSpectral<double, 4>>(m, "LinearSpectral")
        .def("__init__", [](LinearSpectral<double, 4>* self,
                            basix::FiniteElement<double> element,
                            std::shared_ptr<mesh::Mesh<double>> mesh,
                            std::shared_ptr<mesh::MeshTags<std::int32_t>> facet_tags,
                            std::shared_ptr<fem::Function<double>> speed_of_sound,
                            std::shared_ptr<fem::Function<double>> density,
                            const double& source_frequency,
                            const double& source_amplitude,
                            const double& source_speed) {
            new (self) LinearSpectral<double, 4>(
                element, mesh, facet_tags, speed_of_sound, density,
                source_frequency, source_amplitude, source_speed);
        },
        nb::arg("element"), nb::arg("mesh"), nb::arg("facet_tags"),
        nb::arg("speed_of_sound"), nb::arg("density"),
        nb::arg("source_frequency"), nb::arg("source_amplitude"),
        nb::arg("source_speed"),
        nb::call_guard<nb::gil_scoped_release>())
        .def("init", &LinearSpectral<double, 4>::init,
             nb::call_guard<nb::gil_scoped_release>())
        .def("rk4", &LinearSpectral<double, 4>::rk4,
             nb::arg("start_time"), nb::arg("final_time"), nb::arg("time_step"),
             nb::call_guard<nb::gil_scoped_release>())
        .def("u_sol", &LinearSpectral<double, 4>::u_sol,
             nb::call_guard<nb::gil_scoped_release>())
        .def("number_of_dofs", &LinearSpectral<double, 4>::number_of_dofs);
}
