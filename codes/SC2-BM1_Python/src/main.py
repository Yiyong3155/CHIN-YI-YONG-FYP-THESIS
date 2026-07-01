import sys
import time
import traceback
import numpy as np
from mpi4py import MPI
from dolfinx import mesh as dmesh, fem, io, geometry
import dolfinx.cpp.mesh
import basix
from fusx import LinearSpectral

try:
    comm     = MPI.COMM_WORLD
    mpi_rank = comm.Get_rank()

    if mpi_rank == 0:
        print("Starting...", flush=True)

    # USER SETTINGS
    source_frequency = 0.5e6
    source_amplitude = 60000.0
    speed_of_sound   = 1500.0
    density          = 1000.0
    domain_length    = 0.12
    degree_of_basis  = 4
    period           = 1.0 / source_frequency
    output_interval  = 100

    enable_region  = True
    region_min     = np.array([-0.08, -0.08, 0.0])
    region_max     = np.array([ 0.08,  0.08, 0.165])
    region_nx      = 61
    region_ny      = 61
    region_nz      = 61
    region_outfile = "region_data.csv"

    with io.XDMFFile(comm, "/home/shared/Nanobind_Linear_PlaneWave/codes/example_14/mesh.xdmf", "r") as fmesh:
        msh      = fmesh.read_mesh(dmesh.GhostMode.shared_facet, "planar_3d_0")
        msh.topology.create_connectivity(2, 3)
        mt_cell  = fmesh.read_meshtags(msh, "planar_3d_0_cells")
        mt_facet = fmesh.read_meshtags(msh, "planar_3d_0_facets")

    if mpi_rank == 0:
        print("Mesh loaded", flush=True)

    tdim         = msh.topology.dim
    num_cells    = msh.topology.index_map(tdim).size_local
    cell_ids     = np.arange(num_cells, dtype=np.int32)
    h_local      = dolfinx.cpp.mesh.h(msh._cpp_object, tdim, cell_ids)
    h_min_global = comm.allreduce(float(h_local.min()), op=MPI.MIN)

    V_DG = fem.functionspace(msh, ("DG", 0))
    c0   = fem.Function(V_DG)
    rho0 = fem.Function(V_DG)
    cells_1 = mt_cell.find(1)
    c0.x.array[cells_1]   = speed_of_sound
    rho0.x.array[cells_1] = density
    c0.x.scatter_forward()
    rho0.x.scatter_forward()

    CFL              = 0.65
    dt               = CFL * h_min_global / (speed_of_sound * degree_of_basis ** 2)
    steps_per_period = int(period / dt) + 1
    dt               = period / steps_per_period
    start_time       = 0.0
    final_time       = domain_length / speed_of_sound + 8.0 / source_frequency
    n_steps          = int((final_time - start_time) / dt) + 1

    sol_element = basix.create_element(
        basix.ElementFamily.P, basix.CellType.hexahedron, degree_of_basis,
        basix.LagrangeVariant.gll_warped, basix.DPCVariant.unset, False)._e

    if mpi_rank == 0:
        print(f"Degree           : {degree_of_basis}",  flush=True)
        print(f"Min mesh size    : {h_min_global:.2e}", flush=True)
        print(f"CFL              : {CFL}",               flush=True)
        print(f"dt               : {dt:.6e}",            flush=True)
        print(f"Total steps      : {n_steps}",           flush=True)

    model = LinearSpectral(
        sol_element, msh._cpp_object, mt_facet._cpp_object,
        c0._cpp_object, rho0._cpp_object,
        source_frequency, source_amplitude, speed_of_sound)

    if mpi_rank == 0:
        print(f"DOFs             : {model.number_of_dofs()}", flush=True)
        print("Initialising...", flush=True)

    model.init()

    valid_pts = valid_cells = valid_idx = None

    if enable_region:
        if mpi_rank == 0:
            print("Setting up region grid...", flush=True)

        xs = np.linspace(region_min[0], region_max[0], region_nx)
        ys = np.linspace(region_min[1], region_max[1], region_ny)
        zs = np.linspace(region_min[2], region_max[2], region_nz)
        gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
        region_grid_pts = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

        bb_tree         = geometry.bb_tree(msh, tdim)
        cell_candidates = geometry.compute_collisions_points(bb_tree, region_grid_pts)
        colliding_cells = geometry.compute_colliding_cells(msh, cell_candidates, region_grid_pts)

        n_pts      = len(region_grid_pts)
        local_cell = np.full(n_pts, -1, dtype=np.int32)
        for i in range(n_pts):
            c = colliding_cells.links(i)
            if len(c) > 0:
                local_cell[i] = c[0]

        rank_claim   = np.where(local_cell >= 0, mpi_rank, 999999).astype(np.int32)
        global_claim = np.empty_like(rank_claim)
        comm.Allreduce(rank_claim, global_claim, op=MPI.MIN)

        valid_pts, valid_cells, valid_idx = [], [], []
        for i in range(n_pts):
            if local_cell[i] >= 0 and global_claim[i] == mpi_rank:
                valid_pts.append(region_grid_pts[i])
                valid_cells.append(local_cell[i])
                valid_idx.append(i)

        valid_pts   = np.array(valid_pts,   dtype=np.float64)
        valid_cells = np.array(valid_cells, dtype=np.int32)
        valid_idx   = np.array(valid_idx,   dtype=np.int32)

        n_valid_global = comm.allreduce(len(valid_pts), op=MPI.SUM)
        if mpi_rank == 0:
            print(f"Region grid: {n_valid_global} / {n_pts} points found", flush=True)

    def write_region_csv(t_now, first=False):
        if not enable_region or valid_pts is None or len(valid_pts) == 0:
            return
        vals = np.zeros((len(valid_pts), 1), dtype=np.float64)
        model.u_sol().eval(valid_pts, valid_cells, vals)
        local_data = np.column_stack([valid_pts, vals[:,0], np.full(len(valid_pts), t_now)])
        all_data   = comm.gather(local_data, root=0)
        if mpi_rank == 0:
            combined = np.vstack(all_data)
            mode = "w" if first else "a"
            with open(region_outfile, mode) as f:
                if first:
                    f.write("x,y,z,pressure,time\n")
                np.savetxt(f, combined, delimiter=",", fmt="%.6e")

    if mpi_rank == 0:
        print("Running RK4...", flush=True)

    u_sol    = model.u_sol()
    vtx_full = io.VTXWriter(comm, "full_output.bp", [u_sol], "BP4")

    t_current = start_time
    seg_time  = dt * output_interval if output_interval > 0 else final_time
    t0_wall   = time.perf_counter()
    frame     = 0
    first_csv = True

    while t_current < final_time - 1e-14:
        t_end = min(t_current + seg_time, final_time)
        model.rk4(t_current, t_end, dt)
        t_current = t_end
        frame += 1
        vtx_full.write(t_current)
        write_region_csv(t_current, first=first_csv)
        first_csv = False
        if mpi_rank == 0:
            print(f"  Frame {frame}: t={t_current:.4e}", flush=True)

    solve_time = time.perf_counter() - t0_wall
    vtx_full.close()

    if mpi_rank == 0:
        print(f"Solve time       : {solve_time:.4f} s",           flush=True)
        print(f"Time per step    : {solve_time / n_steps:.6f} s", flush=True)
        if enable_region:
            print(f"Region CSV       : {region_outfile}",          flush=True)

        timings = dict(model.get_timings())
        n = int(timings["n_steps"])
        avg = lambda v: v / n if n > 0 else 0.0
        print(f"\n{'Process / Step':<30} {'Total (ms)':>12} {'Per step (ms)':>15}", flush=True)
        print("-" * 60, flush=True)
        for name, key in [
            ("Total f1() time",      "f1_total"),
            ("  Window computation", "window"),
            ("  Boundary condition", "bc_update"),
            ("  Scatter forward",    "scatter_fwd"),
            ("  Vector copy",        "vec_copy"),
            ("  Stiffness operator", "stiffness"),
            ("  Assembly (FEM)",     "assembly"),
            ("  Scatter reverse",    "scatter_rev"),
            ("  Solve (division)",   "solve"),
            ("  Estimated overhead", "overhead"),
        ]:
            print(f"{name:<30} {timings[key]:>12.2f} {avg(timings[key]):>15.4f}", flush=True)
        print("-" * 60, flush=True)
        print(f"{'Total solve time':<30} {timings['rk4_total']:>12.2f} {avg(timings['rk4_total']):>15.4f}", flush=True)
        print(f"{'Steps completed':<30} {n:>12}", flush=True)

except Exception:
    traceback.print_exc()
    sys.exit(1)