#pragma once

#include "forms.h"
#include "spectral_op.hpp"

#include <fstream>
#include <memory>
#include <string>
#include <chrono>

#include <dolfinx.h>
#include <dolfinx/geometry/utils.h>
#include <dolfinx/la/Vector.h>

using namespace dolfinx;

namespace kernels {
template <typename T>
void copy(const la::Vector<T>& in, la::Vector<T>& out) {
  std::span<const T> _in = in.array();
  std::span<T> _out = out.mutable_array();
  std::copy(_in.begin(), _in.end(), _out.begin());
}

template <typename T>
void axpy(la::Vector<T>& r, T alpha, const la::Vector<T>& x, const la::Vector<T>& y) {
  std::transform(x.array().begin(), x.array().begin() + x.index_map()->size_local(), y.array().begin(),
                 r.mutable_array().begin(),
                 [&alpha](const T& vx, const T& vy) { return vx * alpha + vy; });
}
} // namespace kernels

template <typename T, int P>
class LinearSpectral2D {
public:
    LinearSpectral2D(
        basix::FiniteElement<T> element,
        std::shared_ptr<mesh::Mesh<T>> Mesh,
        std::shared_ptr<mesh::MeshTags<std::int32_t>> FacetTags,
        std::shared_ptr<fem::Function<T>> speedOfSound,
        std::shared_ptr<fem::Function<T>> density,
        const T& sourceFrequency,
        const T& sourceAmplitude,
        const T& sourceSpeed
    ) {
        // MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        // Initialize timing
        time_window = 0.0;
        time_bc_update = 0.0;
        time_scatter = 0.0;
        time_copy = 0.0;
        time_stiffness = 0.0;
        time_assemble = 0.0;
        time_scatter_rev = 0.0;
        time_solve = 0.0;
        time_f1_total = 0.0;

        // Physical parameters
        c0 = speedOfSound;
        rho0 = density;
        freq = sourceFrequency;
        w0 = 2 * M_PI * sourceFrequency;
        p0 = sourceAmplitude;
        s0 = sourceSpeed;
        period = 1.0 / sourceFrequency;
        window_length = 4.0;

        // Mesh data
        mesh = Mesh;
        ft = FacetTags;

        // Define function space
        V = std::make_shared<fem::FunctionSpace<T>>(
            fem::create_functionspace(mesh, element));

        // Define field functions
        index_map = V->dofmap()->index_map;
        bs = V->dofmap()->index_map_bs();

        u = std::make_shared<fem::Function<T>>(V);
        u_n = std::make_shared<fem::Function<T>>(V);
        v_n = std::make_shared<fem::Function<T>>(V);

        // Define source function
        g = std::make_shared<fem::Function<T>>(V);
        g_ = g->x()->mutable_array();

        // Define forms
        std::span<T> u_ = u->x()->mutable_array();
        std::fill(u_.begin(), u_.end(), 1.0);

        // Compute exterior facets
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
        fd[fem::IntegralType::exterior_facet].push_back(
            {tag, facet_domains});
        }

        for (auto const& [key, val] : fd) {
        for (auto const& [tag, vec] : val) {
            fd_view[key].push_back({tag, std::span(vec.data(), vec.size())});
            }
        }

        // Define LHS form
        a = std::make_shared<fem::Form<T>>(
            fem::create_form<T>(*form_forms_a, {V}, {{"u", u}, {"c0", c0}, {"rho0", rho0}}, {}, {}, {}));

        m = std::make_shared<la::Vector<T>>(index_map, bs);
        m_ = m->mutable_array();
        std::fill(m_.begin(), m_.end(), 0.0);
        fem::assemble_vector(m_, *a);
        m->scatter_rev(std::plus<T>());

        // Define RHS form
        L = std::make_shared<fem::Form<T>>(fem::create_form<T>(
            *form_forms_L, {V}, {{"g", g}, {"v_n", v_n}, {"c0", c0}, {"rho0", rho0}}, {},
            fd_view, {}, {}));
        b = std::make_shared<la::Vector<T>>(index_map, bs);
        b_ = b->mutable_array();

        // Define operator
        stiff_op = std::make_shared<StiffnessSpectral2D<T, P>>(V);

        // Define coefficient for the operator
        std::span<const T> c0_ = c0->x()->array();
        std::span<const T> rho0_ = rho0->x()->array();

        coeff = std::make_shared<fem::Function<T>>(rho0->function_space());
        c_ = coeff->x()->mutable_array();

        for (std::size_t i = 0; i < rho0_.size(); ++i)
            c_[i] = -1.0 / rho0_[i];

        coeff->x()->scatter_fwd();
        
        // ✅ ADDED: Cache array views to avoid repeated Python binding overhead
        u_n_array_cache = u_n->x()->mutable_array();
        v_n_array_cache = v_n->x()->mutable_array();
        g_array_cache = g->x()->mutable_array();
    }

    void init() {
        u_n->x()->set(0.0);
        v_n->x()->set(0.0);
    }

    void f0(T& t, 
            std::shared_ptr<la::Vector<T>> u, 
            std::shared_ptr<la::Vector<T>> v,
            std::shared_ptr<la::Vector<T>> result) {
        auto start = std::chrono::high_resolution_clock::now();
        kernels::copy<T>(*v, *result);
        auto end = std::chrono::high_resolution_clock::now();
        time_copy += std::chrono::duration<double>(end - start).count();
    }

    void f1(T& t,
            std::shared_ptr<la::Vector<T>> u,
            std::shared_ptr<la::Vector<T>> v,
            std::shared_ptr<la::Vector<T>> result) {
        
        auto f1_start = std::chrono::high_resolution_clock::now();
        
        auto t_start = std::chrono::high_resolution_clock::now();
        
        // Apply windowing
        if (t < period * window_length) {
            window = 0.5 * (1.0 - cos(freq * M_PI * t / window_length));
        } else {
            window = 1.0;
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        time_window += std::chrono::duration<double>(t1 - t_start).count();

        // ✅ MODIFIED: Update boundary condition using cached array
        std::fill(g_array_cache.begin(), g_array_cache.end(), window * p0 * w0 / s0 * cos(w0 * t));
        auto t2 = std::chrono::high_resolution_clock::now();
        time_bc_update += std::chrono::duration<double>(t2 - t1).count();

        u->scatter_fwd();
        auto t3 = std::chrono::high_resolution_clock::now();
        time_scatter += std::chrono::duration<double>(t3 - t2).count();
        
        // ✅ MODIFIED: Copy using cached array instead of u_n->x()
        std::span<const T> u_span = u->array();
        std::copy(u_span.begin(), u_span.end(), u_n_array_cache.begin());
        auto t4 = std::chrono::high_resolution_clock::now();
        time_copy += std::chrono::duration<double>(t4 - t3).count();

        v->scatter_fwd();
        auto t5 = std::chrono::high_resolution_clock::now();
        time_scatter += std::chrono::duration<double>(t5 - t4).count();
        
        // ✅ MODIFIED: Copy using cached array instead of v_n->x()
        std::span<const T> v_span = v->array();
        std::copy(v_span.begin(), v_span.end(), v_n_array_cache.begin());
        auto t6 = std::chrono::high_resolution_clock::now();
        time_copy += std::chrono::duration<double>(t6 - t5).count();

        // Assemble RHS
        std::fill(b_.begin(), b_.end(), 0.0);
        auto t7 = std::chrono::high_resolution_clock::now();
        
        stiff_op->operator()(*u_n->x(), c_, *b);
        auto t8 = std::chrono::high_resolution_clock::now();
        time_stiffness += std::chrono::duration<double>(t8 - t7).count();
        
        fem::assemble_vector(b_, *L);
        auto t9 = std::chrono::high_resolution_clock::now();
        time_assemble += std::chrono::duration<double>(t9 - t8).count();
        
        b->scatter_rev(std::plus<T>());
        auto t10 = std::chrono::high_resolution_clock::now();
        time_scatter_rev += std::chrono::duration<double>(t10 - t9).count();

        // Solve
        {
            out = result->mutable_array();
            _b = b->array();
            _m = m->array();

            std::transform(_b.begin(), _b.end(), _m.begin(), out.begin(),
                           [](const T& bi, const T& mi) { return bi / mi; });
        }
        auto t11 = std::chrono::high_resolution_clock::now();
        time_solve += std::chrono::duration<double>(t11 - t10).count();
        
        auto f1_end = std::chrono::high_resolution_clock::now();
        time_f1_total += std::chrono::duration<double>(f1_end - f1_start).count();
    }

    void rk4(const T& startTime, const T& finalTime, const T& timeStep) {
        T t = startTime;
        T tf = finalTime;
        T dt = timeStep;
        int totalStep = (finalTime - startTime) / timeStep + 1;
        int step = 0;

        std::shared_ptr<la::Vector<T>> u_, v_, un, vn, u0, v0, ku, kv;

        u_ = std::make_shared<la::Vector<T>>(index_map, bs);
        v_ = std::make_shared<la::Vector<T>>(index_map, bs);

        kernels::copy<T>(*u_n->x(), *u_);
        kernels::copy<T>(*v_n->x(), *v_);

        un = std::make_shared<la::Vector<T>>(index_map, bs);
        vn = std::make_shared<la::Vector<T>>(index_map, bs);
        u0 = std::make_shared<la::Vector<T>>(index_map, bs);
        v0 = std::make_shared<la::Vector<T>>(index_map, bs);
        ku = std::make_shared<la::Vector<T>>(index_map, bs);
        kv = std::make_shared<la::Vector<T>>(index_map, bs);

        kernels::copy<T>(*u_, *ku);
        kernels::copy<T>(*v_, *kv);

        std::array<T, 4> a_runge = {0.0, 0.5, 0.5, 1.0};
        std::array<T, 4> b_runge = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
        std::array<T, 4> c_runge = {0.0, 0.5, 0.5, 1.0};

        T tn;

        while (t < tf) {
            dt = std::min(dt, tf - t);

            kernels::copy<T>(*u_, *u0);
            kernels::copy<T>(*v_, *v0);

            for (int i = 0; i < 4; i++) {
                kernels::copy<T>(*u0, *un);
                kernels::copy<T>(*v0, *vn);

                kernels::axpy<T>(*un, dt * a_runge[i], *ku, *un);
                kernels::axpy<T>(*vn, dt * a_runge[i], *kv, *vn);

                tn = t + c_runge[i] * dt;

                f0(tn, un, vn, ku);
                f1(tn, un, vn, kv);

                kernels::axpy<T>(*u_, dt * b_runge[i], *ku, *u_);
                kernels::axpy<T>(*v_, dt * b_runge[i], *kv, *v_);
            }
            
            t += dt;
            step += 1;

            if (step % 100 == 0) {
                if (mpi_rank == 0) {
                std::cout << "t: " << t << ",\t Steps: " << step << "/" << totalStep << std::endl;
                }
            }
        }

        kernels::copy<T>(*u_, *u_n->x());
        kernels::copy<T>(*v_, *v_n->x());
        u_n->x()->scatter_fwd();
        v_n->x()->scatter_fwd();
        
        // Print timing breakdown
        if (mpi_rank == 0) {
            double sum_subtimers = time_window + time_bc_update + time_scatter + 
                                   time_copy + time_stiffness + time_assemble + 
                                   time_scatter_rev + time_solve;
            
            std::cout << "\n========== DETAILED TIMING BREAKDOWN ==========\n";
            std::cout << "f1() total time:          " << time_f1_total * 1000 << " ms\n";
            std::cout << "  Window computation:     " << time_window * 1000 << " ms\n";
            std::cout << "  BC update (std::fill):  " << time_bc_update * 1000 << " ms\n";
            std::cout << "  Scatter forward:        " << time_scatter * 1000 << " ms\n";
            std::cout << "  Vector copy:            " << time_copy * 1000 << " ms\n";
            std::cout << "  Stiffness operator:     " << time_stiffness * 1000 << " ms\n";
            std::cout << "  Assembly (fem):         " << time_assemble * 1000 << " ms\n";
            std::cout << "  Scatter reverse:        " << time_scatter_rev * 1000 << " ms\n";
            std::cout << "  Solve (division):       " << time_solve * 1000 << " ms\n";
            std::cout << "f1() overhead:            " << (time_f1_total - sum_subtimers) * 1000 << " ms\n";
            std::cout << "===============================================\n";
        }
    }

    std::shared_ptr<fem::Function<T>> u_sol() const { return u_n; }

private:
    int mpi_rank, mpi_size;
    int bs;
    T freq, p0, w0, s0, period, window_length, window;

    std::shared_ptr<mesh::Mesh<T>> mesh;
    std::shared_ptr<mesh::MeshTags<std::int32_t>> ft;
    std::shared_ptr<const common::IndexMap> index_map;
    std::shared_ptr<fem::FunctionSpace<T>> V;
    std::shared_ptr<fem::Function<T>> u, u_n, v_n, g, c0, rho0;

    std::shared_ptr<fem::Form<T>> a, L;
    std::shared_ptr<la::Vector<T>> m, b;

    std::span<T> g_, m_, b_, out;
    std::span<const T> _m, _b;

    std::shared_ptr<StiffnessSpectral2D<T, P>> stiff_op;
    std::shared_ptr<fem::Function<T>> coeff;
    std::span<T> c_;

    // Timing accumulators
    double time_window;
    double time_bc_update;
    double time_scatter;
    double time_copy;
    double time_stiffness;
    double time_assemble;
    double time_scatter_rev;
    double time_solve;
    double time_f1_total;
    
    // ✅ ADDED: Cached array views to avoid Python binding overhead
    std::span<T> u_n_array_cache;
    std::span<T> v_n_array_cache;
    std::span<T> g_array_cache;
};