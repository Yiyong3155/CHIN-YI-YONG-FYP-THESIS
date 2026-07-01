#include "Linear.hpp"
#include "forms.h"
#include <cmath>
#include <dolfinx.h>
#include <dolfinx/fem/Constant.h>
#include <dolfinx/io/XDMFFile.h>
#include <iomanip>
#include <iostream>

#define T_MPI MPI_DOUBLE
using T = double;

int main(int argc, char* argv[]) {
  dolfinx::init_logging(argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  {
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    const T sourceFrequency = 0.5e6;
    const T sourceAmplitude = 60000;
    const T period = 1 / sourceFrequency;
    const T speedOfSound = 1500;
    const T density = 1000;
    const T domainLength = 0.12;
    const int degreeOfBasis = 4;

    auto element = fem::CoordinateElement<T>(mesh::CellType::hexahedron, 1);
    io::XDMFFile fmesh(MPI_COMM_WORLD,
        "/home/shared/Nanobind_Linear_PlaneWave/codes/example_14/mesh.xdmf", "r");
    auto mesh = std::make_shared<mesh::Mesh<T>>(
        fmesh.read_mesh(element, mesh::GhostMode::none, "planar_3d_0"));
    mesh->topology()->create_connectivity(2, 3);
    auto mt_cell = std::make_shared<mesh::MeshTags<std::int32_t>>(
        fmesh.read_meshtags(*mesh, "planar_3d_0_cells"));
    auto mt_facet = std::make_shared<mesh::MeshTags<std::int32_t>>(
        fmesh.read_meshtags(*mesh, "planar_3d_0_facets"));

    const int tdim = mesh->topology()->dim();
    const int num_cell = mesh->topology()->index_map(tdim)->size_local();
    std::vector<int> num_cell_range(num_cell);
    std::iota(num_cell_range.begin(), num_cell_range.end(), 0);
    std::vector<T> mesh_size_local = mesh::h(*mesh, num_cell_range, tdim);
    auto min_it = std::min_element(mesh_size_local.begin(), mesh_size_local.end());
    T meshSizeMinLocal = *min_it;
    T meshSizeMinGlobal;
    MPI_Reduce(&meshSizeMinLocal, &meshSizeMinGlobal, 1, T_MPI, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Bcast(&meshSizeMinGlobal, 1, T_MPI, 0, MPI_COMM_WORLD);

    // DG0 space for material parameters
    auto V_DG = std::make_shared<fem::FunctionSpace<T>>(
        fem::create_functionspace(
            mesh,
            basix::create_element<T>(
                basix::element::family::P,
                mesh::cell_type_to_basix_type(mesh::CellType::hexahedron),
                0,
                basix::element::lagrange_variant::unset,
                basix::element::dpc_variant::unset,
                true)));
    auto c0   = std::make_shared<fem::Function<T>>(V_DG);
    auto rho0 = std::make_shared<fem::Function<T>>(V_DG);

    auto cells_1 = mt_cell->find(1);
    std::span<T> c0_arr = c0->x()->mutable_array();
    std::for_each(cells_1.begin(), cells_1.end(),
                  [&](std::int32_t& i) { c0_arr[i] = speedOfSound; });
    c0->x()->scatter_fwd();
    std::span<T> rho0_arr = rho0->x()->mutable_array();
    std::for_each(cells_1.begin(), cells_1.end(),
                  [&](std::int32_t& i) { rho0_arr[i] = density; });
    rho0->x()->scatter_fwd();

    const T CFL = 0.65;
    T timeStepSize = CFL * meshSizeMinGlobal / (speedOfSound * degreeOfBasis * degreeOfBasis);
    const int stepPerPeriod = period / timeStepSize + 1;
    timeStepSize = period / stepPerPeriod;
    const T startTime  = 0.0;
    const T finalTime  = domainLength / speedOfSound + 8.0 / sourceFrequency;
    const int numberOfStep = (finalTime - startTime) / timeStepSize + 1;

    // Spectral element for solution space
    auto sol_element = basix::create_element<T>(
        basix::element::family::P,
        mesh::cell_type_to_basix_type(mesh::CellType::hexahedron),
        degreeOfBasis,
        basix::element::lagrange_variant::gll_warped,
        basix::element::dpc_variant::unset,
        false);

    auto model = LinearSpectral<T, 4>(
        sol_element, mesh, mt_facet, c0, rho0,
        sourceFrequency, sourceAmplitude, speedOfSound);

    auto nDofs = model.number_of_dofs();

    if (mpi_rank == 0) {
      std::cout << "Benchmark: 1\nSource: 2\n";
      std::cout << "Polynomial basis degree: " << degreeOfBasis << "\n";
      std::cout << "Minimum mesh size: " << std::setprecision(2) << meshSizeMinGlobal << "\n";
      std::cout << "Degrees of freedom: " << nDofs << "\n";
      std::cout << "CFL number: " << CFL << "\n";
      std::cout << "Time step size: " << timeStepSize << "\n";
      std::cout << "Number of steps per period: " << stepPerPeriod << "\n";
      std::cout << "Total number of steps: " << numberOfStep << "\n";
    }

    common::Timer tsolve("Solve time");
    model.init();
    tsolve.start();
    model.rk4(startTime, finalTime, timeStepSize);
    tsolve.stop();

    // Write final solution to VTX
    io::VTXWriter<T> vtx(MPI_COMM_WORLD, "output.bp", {model.u_sol()}, "BP4");
    vtx.write(finalTime);

    if (mpi_rank == 0) {
      std::cout << "Solve time: " << tsolve.elapsed()[0] << std::endl;
      std::cout << "Time per step: " << tsolve.elapsed()[0] / numberOfStep << std::endl;
    }
  }
  PetscFinalize();
  return 0;
}
