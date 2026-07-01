#include "Linear.hpp"
 
#include <cstddef>
#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
 
using namespace dolfinx;
namespace nb = nanobind;
 
NB_MODULE(_fusx, m) {
 
    // ===== 2D Linear Spectral Solver =====
    nb::class_<LinearSpectral2D<double, 4>>(m, "LinearSpectral2D")
        .def("__init__", [](LinearSpectral2D<double, 4>* self,
                            basix::FiniteElement<double> element,
                            std::shared_ptr<mesh::Mesh<double>> mesh,
                            std::shared_ptr<mesh::MeshTags<std::int32_t>> facet_tags,
                            std::shared_ptr<fem::Function<double>> speed_of_sound,
                            std::shared_ptr<fem::Function<double>> density,
                            const double& source_frequency,
                            const double& source_amplitude,
                            const double& source_speed,
                            const std::array<double, 2>& region_min,
                            const std::array<double, 2>& region_max,
                            std::size_t grid_nx,
                            std::size_t grid_ny) {
            new (self) LinearSpectral2D<double, 4>(
                element, mesh, facet_tags, speed_of_sound, density,
                source_frequency, source_amplitude, source_speed,
                region_min, region_max, grid_nx, grid_ny);
        },
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
        nb::call_guard<nb::gil_scoped_release>())
        .def("init", &LinearSpectral2D<double, 4>::init,
             nb::call_guard<nb::gil_scoped_release>())
        .def("rk4", &LinearSpectral2D<double, 4>::rk4,
             nb::arg("start_time"),
             nb::arg("final_time"),
             nb::arg("time_step"),
             nb::arg("output_interval") = 0,
             nb::arg("vtx_filename") = "region_output",
             nb::call_guard<nb::gil_scoped_release>(),
             "RK4 time integration with optional VTX output on region submesh")
        .def("u_sol", &LinearSpectral2D<double, 4>::u_sol,
             nb::call_guard<nb::gil_scoped_release>());
 
    // ===== 3D Linear Spectral Solver =====
    nb::class_<LinearSpectral3D<double, 4>>(m, "LinearSpectral3D")
        .def("__init__", [](LinearSpectral3D<double, 4>* self,
                            basix::FiniteElement<double> element,
                            std::shared_ptr<mesh::Mesh<double>> mesh,
                            std::shared_ptr<mesh::MeshTags<std::int32_t>> facet_tags,
                            std::shared_ptr<fem::Function<double>> speed_of_sound,
                            std::shared_ptr<fem::Function<double>> density,
                            const double& source_frequency,
                            const double& source_amplitude,
                            const double& source_speed,
                            const std::array<double, 3>& region_min,
                            const std::array<double, 3>& region_max,
                            std::size_t grid_nx,
                            std::size_t grid_ny,
                            std::size_t grid_nz) {
            new (self) LinearSpectral3D<double, 4>(
                element, mesh, facet_tags, speed_of_sound, density,
                source_frequency, source_amplitude, source_speed,
                region_min, region_max, grid_nx, grid_ny, grid_nz);
        },
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
        nb::call_guard<nb::gil_scoped_release>())
        .def("init", &LinearSpectral3D<double, 4>::init,
             nb::call_guard<nb::gil_scoped_release>())
        .def("rk4", &LinearSpectral3D<double, 4>::rk4,
             nb::arg("start_time"),
             nb::arg("final_time"),
             nb::arg("time_step"),
             nb::arg("output_interval") = 0,
             nb::arg("vtx_filename") = "region_output",
             nb::call_guard<nb::gil_scoped_release>(),
             "RK4 time integration with optional VTX output on region submesh")
        .def("u_sol", &LinearSpectral3D<double, 4>::u_sol,
             nb::call_guard<nb::gil_scoped_release>())
        .def("number_of_dofs", &LinearSpectral3D<double, 4>::number_of_dofs);
}
 