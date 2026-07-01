// Copyright (C) 2024 Adeeb Arif Kor
// SPDX-License-Identifier:    MIT
 
#pragma once
 
#include "forms.h"
#include "spectral_op.hpp"
#include "Linear.hpp"
 
#include <fstream>
#include <memory>
#include <string>
#include <optional>
#include <filesystem>
 
#include <dolfinx.h>
#include <dolfinx/geometry/utils.h>
#include <dolfinx/la/Vector.h>
#include <dolfinx/io/ADIOS2Writers.h>
 
using namespace dolfinx;
 
template <typename T, int P>
class LossySpectral2D {
public:
  LossySpectral2D(basix::FiniteElement<T> element,
                  std::shared_ptr<mesh::Mesh<T>> Mesh,
                  std::shared_ptr<mesh::MeshTags<std::int32_t>> FacetTags,
                  std::shared_ptr<fem::Function<T>> speedOfSound,
                  std::shared_ptr<fem::Function<T>> density,
                  std::shared_ptr<fem::Function<T>> diffusivityOfSound,
                  const T& sourceFrequency,
                  const T& sourceAmplitude,
                  const T& sourceSpeed,
                  const std::array<T, 2>& eval_region_min = {-2.25, -2.25},
                  const std::array<T, 2>& eval_region_max = {2.25, 2.25},
                  std::size_t eval_grid_nx = 251,
                  std::size_t eval_grid_ny = 251) {
 
    element_ = std::optional<basix::FiniteElement<T>>(element);
    region_min = eval_region_min;
    region_max = eval_region_max;
    grid_nx = eval_grid_nx;
    grid_ny = eval_grid_ny;
 
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
 
    c0 = speedOfSound; rho0 = density; delta0 = diffusivityOfSound;
    freq = sourceFrequency;
    w0 = 2 * M_PI * sourceFrequency;
    p0 = sourceAmplitude; s0 = sourceSpeed;
    period = 1.0 / sourceFrequency;
    window_length = 4.0;
    mesh = Mesh; ft = FacetTags;
 
    V = std::make_shared<fem::FunctionSpace<T>>(
        fem::create_functionspace(mesh, element));
    index_map = V->dofmap()->index_map;
    bs = V->dofmap()->index_map_bs();
 
    u = std::make_shared<fem::Function<T>>(V);
    u_n = std::make_shared<fem::Function<T>>(V);
    v_n = std::make_shared<fem::Function<T>>(V);
    g = std::make_shared<fem::Function<T>>(V);
    g_ = g->x()->mutable_array();
    dg = std::make_shared<fem::Function<T>>(V);
    dg_ = dg->x()->mutable_array();
 
    std::span<T> u_ = u->x()->mutable_array();
    std::fill(u_.begin(), u_.end(), 1.0);
 
    std::vector<std::int32_t> ft_unique(ft->values().size());
    std::copy(ft->values().begin(), ft->values().end(), ft_unique.begin());
    std::sort(ft_unique.begin(), ft_unique.end());
    auto it = std::unique(ft_unique.begin(), ft_unique.end());
    ft_unique.erase(it, ft_unique.end());
 
    std::map<fem::IntegralType, std::vector<std::pair<std::int32_t, std::vector<std::int32_t>>>> fd;
    std::map<fem::IntegralType, std::vector<std::pair<std::int32_t, std::span<const std::int32_t>>>> fd_view;
 
    std::vector<std::int32_t> facet_domains;
    for (auto& tag : ft_unique) {
      facet_domains = fem::compute_integration_domains(
        fem::IntegralType::exterior_facet, *V->mesh()->topology_mutable(),
        ft->find(tag), mesh->topology()->dim()-1);
      fd[fem::IntegralType::exterior_facet].push_back({tag, facet_domains});
    }
    for (auto const& [key, val] : fd)
      for (auto const& [tag, vec] : val)
        fd_view[key].push_back({tag, std::span(vec.data(), vec.size())});
 
    a = std::make_shared<fem::Form<T>>(fem::create_form<T>(
        *form_forms_a, {V}, {{"u", u}, {"c0", c0}, {"rho0", rho0}, {"delta0", delta0}}, {},
        fd_view, {}));
    m = std::make_shared<la::Vector<T>>(index_map, bs);
    m_ = m->mutable_array();
    std::fill(m_.begin(), m_.end(), 0.0);
    fem::assemble_vector(m_, *a);
    m->scatter_rev(std::plus<T>());
 
    L = std::make_shared<fem::Form<T>>(fem::create_form<T>(
        *form_forms_L, {V},
        {{"g", g}, {"dg", dg}, {"v_n", v_n}, {"c0", c0}, {"rho0", rho0}, {"delta0", delta0}}, {},
        fd_view, {}, {}));
    b = std::make_shared<la::Vector<T>>(index_map, bs);
    b_ = b->mutable_array();
 
    lin_op = std::make_shared<StiffnessSpectral2D<T, P>>(V);
    att_op = std::make_shared<StiffnessSpectral2D<T, P>>(V);
 
    std::span<const T> c0_ = c0->x()->array();
    std::span<const T> rho0_ = rho0->x()->array();
    std::span<const T> delta0_ = delta0->x()->array();
 
    lin_coeff = std::make_shared<fem::Function<T>>(rho0->function_space());
    lin_coeff_ = lin_coeff->x()->mutable_array();
    att_coeff = std::make_shared<fem::Function<T>>(rho0->function_space());
    att_coeff_ = att_coeff->x()->mutable_array();
 
    for (std::size_t i = 0; i < rho0_.size(); ++i) {
      lin_coeff_[i] = -1.0 / rho0_[i];
      att_coeff_[i] = -delta0_[i] / rho0_[i] / c0_[i] / c0_[i];
    }
    lin_coeff->x()->scatter_fwd();
    att_coeff->x()->scatter_fwd();
  }
 
  void init() { u_n->x()->set(0.0); v_n->x()->set(0.0); }
 
  void f0(T& t, std::shared_ptr<la::Vector<T>> u,
          std::shared_ptr<la::Vector<T>> v, std::shared_ptr<la::Vector<T>> result) {
    kernels::copy<T>(*v, *result);
  }
 
  void f1(T& t, std::shared_ptr<la::Vector<T>> u,
          std::shared_ptr<la::Vector<T>> v, std::shared_ptr<la::Vector<T>> result) {
    if (t < period * window_length) {
      window = 0.5 * (1.0 - cos(freq * M_PI * t / window_length));
      dwindow = 0.5 * M_PI * freq / window_length * sin(freq * M_PI * t / window_length);
    } else {
      window = 1.0; dwindow = 0.0;
    }
 
    std::fill(g_.begin(), g_.end(), window * 2.0 * p0 * w0 / s0 * cos(w0 * t));
    std::fill(dg_.begin(), dg_.end(),
              dwindow * 2.0 * p0 * w0 / s0 * cos(w0 * t)
                  - window * 2.0 * p0 * w0 * w0 / s0 * sin(w0 * t));
 
    u->scatter_fwd(); kernels::copy<T>(*u, *u_n->x());
    v->scatter_fwd(); kernels::copy<T>(*v, *v_n->x());
 
    std::fill(b_.begin(), b_.end(), 0.0);
    lin_op->operator()(*u_n->x(), lin_coeff_, *b);
    att_op->operator()(*v_n->x(), att_coeff_, *b);
    fem::assemble_vector(b_, *L);
    b->scatter_rev(std::plus<T>());
 
    out = result->mutable_array(); _b = b->array(); _m = m->array();
    std::transform(_b.begin(), _b.end(), _m.begin(), out.begin(),
                   [](const T& bi, const T& mi) { return bi / mi; });
  }
 
  void rk4(const T& startTime, const T& finalTime, const T& timeStep,
           int output_interval = 0,
           const std::string& vtx_filename = "lossy_region") {
 
    // ── Build submesh from custom region ─────────────────────────────────
    const std::size_t Nx = grid_nx, Ny = grid_ny;
    std::vector<double> point_coordinates(3 * Nx * Ny);
    for (std::size_t i = 0; i < Nx; ++i)
      for (std::size_t j = 0; j < Ny; ++j) {
        std::size_t idx = i * Ny + j;
        point_coordinates[3*idx]     = region_min[0] + i * (region_max[0] - region_min[0]) / (Nx - 1);
        point_coordinates[3*idx + 1] = region_min[1] + j * (region_max[1] - region_min[1]) / (Ny - 1);
        point_coordinates[3*idx + 2] = 0.0;
      }
 
    const int tdim = mesh->topology()->dim();
    mesh->topology()->create_entities(tdim);
    auto map = mesh->topology()->index_map(tdim);
    const std::int32_t num_entities = map->size_local() + map->num_ghosts();
    std::vector<std::int32_t> entities(num_entities);
    std::iota(entities.begin(), entities.end(), 0);
 
    auto bb_tree = geometry::BoundingBoxTree(*mesh, tdim, entities);
    auto cell_candidates = compute_collisions<double>(bb_tree, point_coordinates);
    auto colliding_cells = geometry::compute_colliding_cells<double>(
        *mesh, cell_candidates, point_coordinates);
 
    std::vector<std::int32_t> cells;
    std::vector<double> points_on_proc;
    for (std::size_t i = 0; i < Nx * Ny; ++i) {
      auto link = colliding_cells.links(i);
      if (link.size() > 0) {
        points_on_proc.push_back(point_coordinates[3*i]);
        points_on_proc.push_back(point_coordinates[3*i + 1]);
        points_on_proc.push_back(point_coordinates[3*i + 2]);
        cells.push_back(link[0]);
      }
    }
 
    std::vector<std::int32_t> submesh_cells = cells;
    std::sort(submesh_cells.begin(), submesh_cells.end());
    submesh_cells.erase(std::unique(submesh_cells.begin(), submesh_cells.end()), submesh_cells.end());
 
    auto [sub_mesh, sub_cell_map, sub_vertex_map, sub_geom_map] =
        mesh::create_submesh(*mesh, tdim,
            std::span<const std::int32_t>(submesh_cells.data(), submesh_cells.size()));
    auto submesh_ptr = std::make_shared<mesh::Mesh<T>>(std::move(sub_mesh));
    auto V_sub = std::make_shared<fem::FunctionSpace<T>>(
        fem::create_functionspace(submesh_ptr, element_.value()));
    auto u_sub = std::make_shared<fem::Function<T>>(V_sub);
    u_sub->x()->set(0.0);
 
    std::unique_ptr<io::VTXWriter<T>> region_writer;
    if (output_interval > 0) {
      std::filesystem::create_directories("vtx_output");
      region_writer = std::make_unique<io::VTXWriter<T>>(
          MPI_COMM_WORLD,
          std::filesystem::path("vtx_output") / (vtx_filename + ".bp"),
          dolfinx::io::adios2_writer::U<T>{u_sub}, "BP4");
    }
 
    if (mpi_rank == 0) {
      std::cout << "\n📊 Lossy Submesh VTX output setup (2D):\n";
      std::cout << "   Grid: " << Nx << " x " << Ny << " = " << (Nx*Ny) << " points\n";
      std::cout << "   Submesh cells: " << submesh_cells.size() << "\n";
      std::cout << "   Region: [" << region_min[0] << ", " << region_min[1]
                << "] to [" << region_max[0] << ", " << region_max[1] << "]\n";
      if (output_interval > 0)
        std::cout << "   Output file: vtx_output/" << vtx_filename << ".bp\n";
    }
 
    // ── Time-stepping ─────────────────────────────────────────────────────
    T t = startTime, tf = finalTime, dt = timeStep;
    int totalStep = static_cast<int>((finalTime - startTime) / timeStep) + 1;
    int step = 0;
 
    auto u_ = std::make_shared<la::Vector<T>>(index_map, bs);
    auto v_ = std::make_shared<la::Vector<T>>(index_map, bs);
    kernels::copy<T>(*u_n->x(), *u_); kernels::copy<T>(*v_n->x(), *v_);
 
    auto un = std::make_shared<la::Vector<T>>(index_map, bs);
    auto vn = std::make_shared<la::Vector<T>>(index_map, bs);
    auto u0 = std::make_shared<la::Vector<T>>(index_map, bs);
    auto v0 = std::make_shared<la::Vector<T>>(index_map, bs);
    auto ku = std::make_shared<la::Vector<T>>(index_map, bs);
    auto kv = std::make_shared<la::Vector<T>>(index_map, bs);
    kernels::copy<T>(*u_, *ku); kernels::copy<T>(*v_, *kv);
 
    std::array<T, 4> a_runge = {0.0, 0.5, 0.5, 1.0};
    std::array<T, 4> b_runge = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
    std::array<T, 4> c_runge = {0.0, 0.5, 0.5, 1.0};
    T tn;
 
    if (region_writer) {
      kernels::copy(*u_, *u_n->x()); u_n->x()->scatter_fwd();
      {
        std::vector<std::int32_t> sub_cells_idx(submesh_cells.size());
        std::iota(sub_cells_idx.begin(), sub_cells_idx.end(), 0);
        u_sub->interpolate(*u_n,
            std::span<const std::int32_t>(sub_cell_map.data(), sub_cell_map.size()),
            std::span<const std::int32_t>(sub_cells_idx.data(), sub_cells_idx.size()));
      }
      region_writer->write(startTime);
    }
 
    while (t < tf) {
      dt = std::min(dt, tf - t);
      kernels::copy<T>(*u_, *u0); kernels::copy<T>(*v_, *v0);
      for (int i = 0; i < 4; ++i) {
        kernels::copy<T>(*u0, *un); kernels::copy<T>(*v0, *vn);
        kernels::axpy<T>(*un, dt * a_runge[i], *ku, *un);
        kernels::axpy<T>(*vn, dt * a_runge[i], *kv, *vn);
        tn = t + c_runge[i] * dt;
        f0(tn, un, vn, ku); f1(tn, un, vn, kv);
        kernels::axpy<T>(*u_, dt * b_runge[i], *ku, *u_);
        kernels::axpy<T>(*v_, dt * b_runge[i], *kv, *v_);
      }
      t += dt; step += 1;
 
      if (region_writer && step % output_interval == 0) {
        kernels::copy(*u_, *u_n->x()); u_n->x()->scatter_fwd();
        {
          std::vector<std::int32_t> sub_cells_idx(submesh_cells.size());
          std::iota(sub_cells_idx.begin(), sub_cells_idx.end(), 0);
          u_sub->interpolate(*u_n,
              std::span<const std::int32_t>(sub_cell_map.data(), sub_cell_map.size()),
              std::span<const std::int32_t>(sub_cells_idx.data(), sub_cells_idx.size()));
        }
        region_writer->write(t);
      }
 
      if (step % 100 == 0 && mpi_rank == 0)
        std::cout << "t: " << t << ",\t Steps: " << step << "/" << totalStep << "\n";
    }
 
    kernels::copy<T>(*u_, *u_n->x()); kernels::copy<T>(*v_, *v_n->x());
    u_n->x()->scatter_fwd(); v_n->x()->scatter_fwd();
 
    if (mpi_rank == 0) {
      std::cout << "\n========================================\n";
      std::cout << "✅ Lossy simulation complete!\n";
      if (output_interval > 0)
        std::cout << "   VTX submesh: vtx_output/" << vtx_filename << ".bp\n";
      std::cout << "   Region: [" << region_min[0] << ", " << region_min[1]
                << "] to [" << region_max[0] << ", " << region_max[1] << "]\n";
      std::cout << "   Submesh cells: " << submesh_cells.size() << "\n";
      std::cout << "========================================\n";
    }
  }
 
  std::shared_ptr<fem::Function<T>> u_sol() const { return u_n; }
  std::int64_t number_of_dofs() const { return V->dofmap()->index_map->size_global(); }
 
private:
  int mpi_rank, mpi_size, bs;
  T freq, p0, w0, s0, period, window_length, window, dwindow;
 
  std::optional<basix::FiniteElement<T>> element_;
  std::array<T, 2> region_min, region_max;
  std::size_t grid_nx, grid_ny;
 
  std::shared_ptr<mesh::Mesh<T>> mesh;
  std::shared_ptr<mesh::MeshTags<std::int32_t>> ft;
  std::shared_ptr<const common::IndexMap> index_map;
  std::shared_ptr<fem::FunctionSpace<T>> V;
  std::shared_ptr<fem::Function<T>> u, u_n, v_n, g, dg, c0, rho0, delta0;
  std::shared_ptr<fem::Form<T>> a, L;
  std::shared_ptr<la::Vector<T>> m, b;
  std::span<T> g_, dg_, m_, b_, out;
  std::span<const T> _m, _b;
  std::shared_ptr<StiffnessSpectral2D<T, P>> lin_op, att_op;
  std::shared_ptr<fem::Function<T>> lin_coeff, att_coeff;
  std::span<T> lin_coeff_, att_coeff_;
};
 
template <typename T, int P>
class LossySpectral3D {
public:
  LossySpectral3D(basix::FiniteElement<T> element,
                  std::shared_ptr<mesh::Mesh<T>> Mesh,
                  std::shared_ptr<mesh::MeshTags<std::int32_t>> FacetTags,
                  std::shared_ptr<fem::Function<T>> speedOfSound,
                  std::shared_ptr<fem::Function<T>> density,
                  std::shared_ptr<fem::Function<T>> diffusivityOfSound,
                  const T& sourceFrequency,
                  const T& sourceAmplitude,
                  const T& sourceSpeed,
                  const std::array<T, 3>& eval_region_min = {-2.25, -2.25, -2.25},
                  const std::array<T, 3>& eval_region_max = {2.25, 2.25, 2.25},
                  std::size_t eval_grid_nx = 101,
                  std::size_t eval_grid_ny = 101,
                  std::size_t eval_grid_nz = 101) {
 
    element_ = std::optional<basix::FiniteElement<T>>(element);
    region_min = eval_region_min;
    region_max = eval_region_max;
    grid_nx = eval_grid_nx; grid_ny = eval_grid_ny; grid_nz = eval_grid_nz;
 
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
 
    c0 = speedOfSound; rho0 = density; delta0 = diffusivityOfSound;
    freq = sourceFrequency;
    w0 = 2 * M_PI * sourceFrequency;
    p0 = sourceAmplitude; s0 = sourceSpeed;
    period = 1.0 / sourceFrequency;
    window_length = 4.0;
    mesh = Mesh; ft = FacetTags;
 
    V = std::make_shared<fem::FunctionSpace<T>>(
        fem::create_functionspace(mesh, element));
    index_map = V->dofmap()->index_map;
    bs = V->dofmap()->index_map_bs();
 
    u = std::make_shared<fem::Function<T>>(V);
    u_n = std::make_shared<fem::Function<T>>(V);
    v_n = std::make_shared<fem::Function<T>>(V);
    g = std::make_shared<fem::Function<T>>(V);
    g_ = g->x()->mutable_array();
    dg = std::make_shared<fem::Function<T>>(V);
    dg_ = dg->x()->mutable_array();
 
    std::span<T> u_ = u->x()->mutable_array();
    std::fill(u_.begin(), u_.end(), 1.0);
 
    std::vector<std::int32_t> ft_unique(ft->values().size());
    std::copy(ft->values().begin(), ft->values().end(), ft_unique.begin());
    std::sort(ft_unique.begin(), ft_unique.end());
    auto it = std::unique(ft_unique.begin(), ft_unique.end());
    ft_unique.erase(it, ft_unique.end());
 
    std::map<fem::IntegralType, std::vector<std::pair<std::int32_t, std::vector<std::int32_t>>>> fd;
    std::map<fem::IntegralType, std::vector<std::pair<std::int32_t, std::span<const std::int32_t>>>> fd_view;
 
    std::vector<std::int32_t> facet_domains;
    for (auto& tag : ft_unique) {
      facet_domains = fem::compute_integration_domains(
        fem::IntegralType::exterior_facet, *V->mesh()->topology_mutable(),
        ft->find(tag), mesh->topology()->dim()-1);
      fd[fem::IntegralType::exterior_facet].push_back({tag, facet_domains});
    }
    for (auto const& [key, val] : fd)
      for (auto const& [tag, vec] : val)
        fd_view[key].push_back({tag, std::span(vec.data(), vec.size())});
 
    a = std::make_shared<fem::Form<T>>(fem::create_form<T>(
        *form_forms_a, {V}, {{"u", u}, {"c0", c0}, {"rho0", rho0}, {"delta0", delta0}}, {},
        fd_view, {}));
    m = std::make_shared<la::Vector<T>>(index_map, bs);
    m_ = m->mutable_array();
    std::fill(m_.begin(), m_.end(), 0.0);
    fem::assemble_vector(m_, *a);
    m->scatter_rev(std::plus<T>());
 
    L = std::make_shared<fem::Form<T>>(fem::create_form<T>(
        *form_forms_L, {V},
        {{"g", g}, {"dg", dg}, {"v_n", v_n}, {"c0", c0}, {"rho0", rho0}, {"delta0", delta0}}, {},
        fd_view, {}, {}));
    b = std::make_shared<la::Vector<T>>(index_map, bs);
    b_ = b->mutable_array();
 
    lin_op = std::make_shared<StiffnessSpectral3D<T, P>>(V);
    att_op = std::make_shared<StiffnessSpectral3D<T, P>>(V);
 
    std::span<const T> c0_ = c0->x()->array();
    std::span<const T> rho0_ = rho0->x()->array();
    std::span<const T> delta0_ = delta0->x()->array();
 
    lin_coeff = std::make_shared<fem::Function<T>>(rho0->function_space());
    lin_coeff_ = lin_coeff->x()->mutable_array();
    att_coeff = std::make_shared<fem::Function<T>>(rho0->function_space());
    att_coeff_ = att_coeff->x()->mutable_array();
 
    for (std::size_t i = 0; i < rho0_.size(); ++i) {
      lin_coeff_[i] = -1.0 / rho0_[i];
      att_coeff_[i] = -delta0_[i] / rho0_[i] / c0_[i] / c0_[i];
    }
    lin_coeff->x()->scatter_fwd();
    att_coeff->x()->scatter_fwd();
  }
 
  void init() { u_n->x()->set(0.0); v_n->x()->set(0.0); }
 
  void f0(T& t, std::shared_ptr<la::Vector<T>> u,
          std::shared_ptr<la::Vector<T>> v, std::shared_ptr<la::Vector<T>> result) {
    kernels::copy<T>(*v, *result);
  }
 
  void f1(T& t, std::shared_ptr<la::Vector<T>> u,
          std::shared_ptr<la::Vector<T>> v, std::shared_ptr<la::Vector<T>> result) {
    if (t < period * window_length) {
      window = 0.5 * (1.0 - cos(freq * M_PI * t / window_length));
      dwindow = 0.5 * M_PI * freq / window_length * sin(freq * M_PI * t / window_length);
    } else {
      window = 1.0; dwindow = 0.0;
    }
 
    std::fill(g_.begin(), g_.end(), window * 2.0 * p0 * w0 / s0 * cos(w0 * t));
    std::fill(dg_.begin(), dg_.end(),
              dwindow * 2.0 * p0 * w0 / s0 * cos(w0 * t)
                  - window * 2.0 * p0 * w0 * w0 / s0 * sin(w0 * t));
 
    u->scatter_fwd(); kernels::copy<T>(*u, *u_n->x());
    v->scatter_fwd(); kernels::copy<T>(*v, *v_n->x());
 
    std::fill(b_.begin(), b_.end(), 0.0);
    lin_op->operator()(*u_n->x(), lin_coeff_, *b);
    att_op->operator()(*v_n->x(), att_coeff_, *b);
    fem::assemble_vector(b_, *L);
    b->scatter_rev(std::plus<T>());
 
    out = result->mutable_array(); _b = b->array(); _m = m->array();
    std::transform(_b.begin(), _b.end(), _m.begin(), out.begin(),
                   [](const T& bi, const T& mi) { return bi / mi; });
  }
 
  void rk4(const T& startTime, const T& finalTime, const T& timeStep,
           int output_interval = 0,
           const std::string& vtx_filename = "lossy_region_3d") {
 
    // ── Build submesh from custom region ─────────────────────────────────
    const std::size_t Nx = grid_nx, Ny = grid_ny, Nz = grid_nz;
    std::vector<double> point_coordinates(3 * Nx * Ny * Nz);
    for (std::size_t i = 0; i < Nx; ++i)
      for (std::size_t j = 0; j < Ny; ++j)
        for (std::size_t k = 0; k < Nz; ++k) {
          std::size_t idx = i * Ny * Nz + j * Nz + k;
          point_coordinates[3*idx]     = region_min[0] + i * (region_max[0] - region_min[0]) / (Nx - 1);
          point_coordinates[3*idx + 1] = region_min[1] + j * (region_max[1] - region_min[1]) / (Ny - 1);
          point_coordinates[3*idx + 2] = region_min[2] + k * (region_max[2] - region_min[2]) / (Nz - 1);
        }
 
    const int tdim = mesh->topology()->dim();
    mesh->topology()->create_entities(tdim);
    auto map = mesh->topology()->index_map(tdim);
    const std::int32_t num_entities = map->size_local() + map->num_ghosts();
    std::vector<std::int32_t> entities(num_entities);
    std::iota(entities.begin(), entities.end(), 0);
 
    auto bb_tree = geometry::BoundingBoxTree(*mesh, tdim, entities);
    auto cell_candidates = compute_collisions<double>(bb_tree, point_coordinates);
    auto colliding_cells = geometry::compute_colliding_cells<double>(
        *mesh, cell_candidates, point_coordinates);
 
    std::vector<std::int32_t> cells;
    std::vector<double> points_on_proc;
    for (std::size_t i = 0; i < Nx * Ny * Nz; ++i) {
      auto link = colliding_cells.links(i);
      if (link.size() > 0) {
        points_on_proc.push_back(point_coordinates[3*i]);
        points_on_proc.push_back(point_coordinates[3*i + 1]);
        points_on_proc.push_back(point_coordinates[3*i + 2]);
        cells.push_back(link[0]);
      }
    }
 
    std::vector<std::int32_t> submesh_cells = cells;
    std::sort(submesh_cells.begin(), submesh_cells.end());
    submesh_cells.erase(std::unique(submesh_cells.begin(), submesh_cells.end()), submesh_cells.end());
 
    auto [sub_mesh, sub_cell_map, sub_vertex_map, sub_geom_map] =
        mesh::create_submesh(*mesh, tdim,
            std::span<const std::int32_t>(submesh_cells.data(), submesh_cells.size()));
    auto submesh_ptr = std::make_shared<mesh::Mesh<T>>(std::move(sub_mesh));
    auto V_sub = std::make_shared<fem::FunctionSpace<T>>(
        fem::create_functionspace(submesh_ptr, element_.value()));
    auto u_sub = std::make_shared<fem::Function<T>>(V_sub);
    u_sub->x()->set(0.0);
 
    std::unique_ptr<io::VTXWriter<T>> region_writer;
    if (output_interval > 0) {
      std::filesystem::create_directories("vtx_output");
      region_writer = std::make_unique<io::VTXWriter<T>>(
          MPI_COMM_WORLD,
          std::filesystem::path("vtx_output") / (vtx_filename + ".bp"),
          dolfinx::io::adios2_writer::U<T>{u_sub}, "BP4");
    }
 
    if (mpi_rank == 0) {
      std::cout << "\n📊 Lossy Submesh VTX output setup (3D):\n";
      std::cout << "   Grid: " << Nx << " x " << Ny << " x " << Nz
                << " = " << (Nx*Ny*Nz) << " points\n";
      std::cout << "   Submesh cells: " << submesh_cells.size() << "\n";
      std::cout << "   Region: [" << region_min[0] << ", " << region_min[1]
                << ", " << region_min[2] << "] to [" << region_max[0]
                << ", " << region_max[1] << ", " << region_max[2] << "]\n";
      if (output_interval > 0)
        std::cout << "   Output file: vtx_output/" << vtx_filename << ".bp\n";
    }
 
    // ── Time-stepping ─────────────────────────────────────────────────────
    T t = startTime, tf = finalTime, dt = timeStep;
    int totalStep = static_cast<int>((finalTime - startTime) / timeStep) + 1;
    int step = 0;
 
    auto u_ = std::make_shared<la::Vector<T>>(index_map, bs);
    auto v_ = std::make_shared<la::Vector<T>>(index_map, bs);
    kernels::copy<T>(*u_n->x(), *u_); kernels::copy<T>(*v_n->x(), *v_);
 
    auto un = std::make_shared<la::Vector<T>>(index_map, bs);
    auto vn = std::make_shared<la::Vector<T>>(index_map, bs);
    auto u0 = std::make_shared<la::Vector<T>>(index_map, bs);
    auto v0 = std::make_shared<la::Vector<T>>(index_map, bs);
    auto ku = std::make_shared<la::Vector<T>>(index_map, bs);
    auto kv = std::make_shared<la::Vector<T>>(index_map, bs);
    kernels::copy<T>(*u_, *ku); kernels::copy<T>(*v_, *kv);
 
    std::array<T, 4> a_runge = {0.0, 0.5, 0.5, 1.0};
    std::array<T, 4> b_runge = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
    std::array<T, 4> c_runge = {0.0, 0.5, 0.5, 1.0};
    T tn;
 
    if (region_writer) {
      kernels::copy(*u_, *u_n->x()); u_n->x()->scatter_fwd();
      {
        std::vector<std::int32_t> sub_cells_idx(submesh_cells.size());
        std::iota(sub_cells_idx.begin(), sub_cells_idx.end(), 0);
        u_sub->interpolate(*u_n,
            std::span<const std::int32_t>(sub_cell_map.data(), sub_cell_map.size()),
            std::span<const std::int32_t>(sub_cells_idx.data(), sub_cells_idx.size()));
      }
      region_writer->write(startTime);
    }
 
    while (t < tf) {
      dt = std::min(dt, tf - t);
      kernels::copy<T>(*u_, *u0); kernels::copy<T>(*v_, *v0);
      for (int i = 0; i < 4; ++i) {
        kernels::copy<T>(*u0, *un); kernels::copy<T>(*v0, *vn);
        kernels::axpy<T>(*un, dt * a_runge[i], *ku, *un);
        kernels::axpy<T>(*vn, dt * a_runge[i], *kv, *vn);
        tn = t + c_runge[i] * dt;
        f0(tn, un, vn, ku); f1(tn, un, vn, kv);
        kernels::axpy<T>(*u_, dt * b_runge[i], *ku, *u_);
        kernels::axpy<T>(*v_, dt * b_runge[i], *kv, *v_);
      }
      t += dt; step += 1;
 
      if (region_writer && step % output_interval == 0) {
        kernels::copy(*u_, *u_n->x()); u_n->x()->scatter_fwd();
        {
          std::vector<std::int32_t> sub_cells_idx(submesh_cells.size());
          std::iota(sub_cells_idx.begin(), sub_cells_idx.end(), 0);
          u_sub->interpolate(*u_n,
              std::span<const std::int32_t>(sub_cell_map.data(), sub_cell_map.size()),
              std::span<const std::int32_t>(sub_cells_idx.data(), sub_cells_idx.size()));
        }
        region_writer->write(t);
      }
 
      if (step % 100 == 0 && mpi_rank == 0)
        std::cout << "t: " << t << ",\t Steps: " << step << "/" << totalStep << "\n";
    }
 
    kernels::copy<T>(*u_, *u_n->x()); kernels::copy<T>(*v_, *v_n->x());
    u_n->x()->scatter_fwd(); v_n->x()->scatter_fwd();
 
    if (mpi_rank == 0) {
      std::cout << "\n========================================\n";
      std::cout << "✅ Lossy 3D simulation complete!\n";
      if (output_interval > 0)
        std::cout << "   VTX submesh: vtx_output/" << vtx_filename << ".bp\n";
      std::cout << "   Region: [" << region_min[0] << ", " << region_min[1]
                << ", " << region_min[2] << "] to [" << region_max[0]
                << ", " << region_max[1] << ", " << region_max[2] << "]\n";
      std::cout << "   Submesh cells: " << submesh_cells.size() << "\n";
      std::cout << "========================================\n";
    }
  }
 
  std::shared_ptr<fem::Function<T>> u_sol() const { return u_n; }
  std::int64_t number_of_dofs() const { return V->dofmap()->index_map->size_global(); }
 
private:
  int mpi_rank, mpi_size, bs;
  T freq, p0, w0, s0, period, window_length, window, dwindow;
 
  std::optional<basix::FiniteElement<T>> element_;
  std::array<T, 3> region_min, region_max;
  std::size_t grid_nx, grid_ny, grid_nz;
 
  std::shared_ptr<mesh::Mesh<T>> mesh;
  std::shared_ptr<mesh::MeshTags<std::int32_t>> ft;
  std::shared_ptr<const common::IndexMap> index_map;
  std::shared_ptr<fem::FunctionSpace<T>> V;
  std::shared_ptr<fem::Function<T>> u, u_n, v_n, g, dg, c0, rho0, delta0;
  std::shared_ptr<fem::Form<T>> a, L;
  std::shared_ptr<la::Vector<T>> m, b;
  std::span<T> g_, dg_, m_, b_, out;
  std::span<const T> _m, _b;
  std::shared_ptr<StiffnessSpectral3D<T, P>> lin_op, att_op;
  std::shared_ptr<fem::Function<T>> lin_coeff, att_coeff;
  std::span<T> lin_coeff_, att_coeff_;
};
 
template <typename T>
const T compute_diffusivity_of_sound(const T w0, const T c0, const T alpha) {
  return 2 * alpha * c0 * c0 * c0 / w0 / w0;
}