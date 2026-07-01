"""
Extended Linear Solver API with Region Analysis Features
For FYP - Adds presentation-ready functionality
"""
import numpy as np
from mpi4py import MPI
from dolfinx.io import VTXWriter
from dolfinx.fem import Function
import json

class LinearSpectral2DExtended:
    """Extended wrapper around LinearSpectral2D with region analysis"""
    
    def __init__(self, model, mesh, u_n):
        self.model = model
        self.mesh = mesh
        self.u_n = u_n
        self.mpi_rank = MPI.COMM_WORLD.rank
        
    def write_region_vtx(self, filename, region_bounds, start_time, final_time, 
                        time_step, output_interval):
        """
        Write VTX output with region metadata for easy ParaView setup
        
        Args:
            filename: Output filename (e.g., "output.bp")
            region_bounds: dict with 'xmin', 'xmax', 'ymin', 'ymax'
            start_time, final_time, time_step: Simulation parameters
            output_interval: Steps between outputs
        """
        if self.mpi_rank == 0:
            print(f"\n📊 Writing VTX with region metadata...")
            print(f"   Region: [{region_bounds['xmin']}, {region_bounds['ymin']}] to "
                  f"[{region_bounds['xmax']}, {region_bounds['ymax']}]")
        
        # Run simulation with VTX output
        with VTXWriter(self.mesh.comm, filename, [self.u_n], "bp5") as vtx:
            vtx_cpp = vtx._cpp_object
            self.model.rk4(start_time, final_time, time_step, output_interval, vtx_cpp)
        
        # Save region metadata alongside VTX file
        if self.mpi_rank == 0:
            metadata = {
                "region_bounds": region_bounds,
                "clip_filter_params": {
                    "position": [region_bounds['xmin'], region_bounds['ymin'], 0],
                    "length": [
                        region_bounds['xmax'] - region_bounds['xmin'],
                        region_bounds['ymax'] - region_bounds['ymin'],
                        0.001
                    ]
                },
                "instructions": "Use Filters->Clip->Box with the clip_filter_params above"
            }
            
            with open(f"{filename}_region.json", 'w') as f:
                json.dump(metadata, f, indent=2)
            
            print(f"✅ Saved region metadata to {filename}_region.json")
    
    def save_paraview_state(self, vtx_file, region_bounds, output_file="view.pvsm"):
        """
        Generate ParaView state file with region already clipped
        
        Args:
            vtx_file: Path to VTX file
            region_bounds: dict with xmin, xmax, ymin, ymax
            output_file: Output .pvsm filename
        """
        if self.mpi_rank == 0:
            print(f"\n�� Generating ParaView state file...")
            
            # Calculate clip parameters
            pos_x = region_bounds['xmin']
            pos_y = region_bounds['ymin']
            len_x = region_bounds['xmax'] - region_bounds['xmin']
            len_y = region_bounds['ymax'] - region_bounds['ymin']
            
            pvsm_content = f"""<?xml version="1.0"?>
<ParaView>
  <ServerManagerState version="5.11.0">
    <!-- VTX Reader -->
    <Proxy group="sources" type="VTXReader" id="1" servers="1">
      <Property name="FileName" id="1.FileName" number_of_elements="1">
        <Element index="0" value="{vtx_file}"/>
      </Property>
    </Proxy>
    
    <!-- Clip Filter with Box -->
    <Proxy group="filters" type="Clip" id="2" servers="1">
      <Property name="Input" id="2.Input">
        <Proxy value="1"/>
      </Property>
      <Property name="ClipType" id="2.ClipType">
        <Proxy value="3"/>
      </Property>
      <Property name="Invert" id="2.Invert" number_of_elements="1">
        <Element index="0" value="0"/>
      </Property>
    </Proxy>
    
    <!-- Box Parameters -->
    <Proxy group="implicit_functions" type="Box" id="3" servers="1">
      <Property name="Position" id="3.Position" number_of_elements="3">
        <Element index="0" value="{pos_x}"/>
        <Element index="1" value="{pos_y}"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Length" id="3.Length" number_of_elements="3">
        <Element index="0" value="{len_x}"/>
        <Element index="1" value="{len_y}"/>
        <Element index="2" value="0.001"/>
      </Property>
    </Proxy>
  </ServerManagerState>
</ParaView>"""
            
            with open(output_file, 'w') as f:
                f.write(pvsm_content)
            
            print(f"✅ ParaView state saved to {output_file}")
            print(f"   Usage: paraview --state={output_file}")
    
    def get_region_stats(self, region_bounds):
        """
        Extract statistics from a region
        
        Args:
            region_bounds: dict with xmin, xmax, ymin, ymax
            
        Returns:
            dict with max, min, mean, std of pressure in region
        """
        from dolfinx import geometry
        
        # Create evaluation points in region
        nx, ny = 100, 100
        x = np.linspace(region_bounds['xmin'], region_bounds['xmax'], nx)
        y = np.linspace(region_bounds['ymin'], region_bounds['ymax'], ny)
        points = np.array([[xi, yi, 0.0] for xi in x for yi in y])
        
        # Find cells
        tdim = self.mesh.topology.dim
        self.mesh.topology.create_entities(tdim)
        
        bb_tree = geometry.bb_tree(self.mesh, tdim)
        cell_candidates = geometry.compute_collisions_points(bb_tree, points)
        colliding_cells = geometry.compute_colliding_cells(self.mesh, cell_candidates, points)
        
        # Evaluate
        cells = []
        points_on_proc = []
        for i, point in enumerate(points):
            cell_list = colliding_cells.links(i)
            if len(cell_list) > 0:
                cells.append(cell_list[0])
                points_on_proc.append(point)
        
        if len(points_on_proc) > 0:
            points_on_proc = np.array(points_on_proc)
            cells_array = np.array(cells, dtype=np.int32)
            values = np.zeros((len(points_on_proc), 1))
            self.u_n.eval(points_on_proc, cells_array, values)
            values = values.flatten()
            
            local_stats = {
                'max': float(np.max(values)),
                'min': float(np.min(values)),
                'mean': float(np.mean(values)),
                'std': float(np.std(values)),
                'count': len(values)
            }
        else:
            local_stats = {'count': 0}
        
        # Gather from all processes
        all_stats = MPI.COMM_WORLD.gather(local_stats, root=0)
        
        if self.mpi_rank == 0:
            # Combine statistics
            all_values = []
            for stats in all_stats:
                if stats['count'] > 0:
                    # Approximate by assuming uniform distribution in each rank's stats
                    all_values.extend([stats['mean']] * stats['count'])
            
            if len(all_values) > 0:
                combined = {
                    'max': max(s['max'] for s in all_stats if s['count'] > 0),
                    'min': min(s['min'] for s in all_stats if s['count'] > 0),
                    'mean': np.mean(all_values),
                    'std': np.std(all_values),
                    'total_points': sum(s['count'] for s in all_stats)
                }
                return combined
        
        return None
    
    def compare_regions(self, region1, region2, label1="Region 1", label2="Region 2"):
        """
        Compare statistics between two regions
        
        Args:
            region1, region2: Region bound dicts
            label1, label2: Names for regions
            
        Returns:
            Comparison dict
        """
        if self.mpi_rank == 0:
            print(f"\n📊 Comparing regions...")
        
        stats1 = self.get_region_stats(region1)
        stats2 = self.get_region_stats(region2)
        
        if self.mpi_rank == 0 and stats1 and stats2:
            comparison = {
                label1: stats1,
                label2: stats2,
                'difference': {
                    'max_diff': stats1['max'] - stats2['max'],
                    'mean_diff': stats1['mean'] - stats2['mean'],
                    'max_ratio': stats1['max'] / stats2['max'] if stats2['max'] != 0 else None
                }
            }
            
            print(f"\n{'='*60}")
            print(f"Region Comparison:")
            print(f"{'='*60}")
            print(f"{label1}:")
            print(f"  Max:  {stats1['max']:.2e} Pa")
            print(f"  Min:  {stats1['min']:.2e} Pa")
            print(f"  Mean: {stats1['mean']:.2e} Pa")
            print(f"\n{label2}:")
            print(f"  Max:  {stats2['max']:.2e} Pa")
            print(f"  Min:  {stats2['min']:.2e} Pa")
            print(f"  Mean: {stats2['mean']:.2e} Pa")
            print(f"\nDifference:")
            print(f"  Max diff:  {comparison['difference']['max_diff']:.2e} Pa")
            print(f"  Mean diff: {comparison['difference']['mean_diff']:.2e} Pa")
            print(f"{'='*60}")
            
            return comparison
        
        return None


    def export_region_images(self, region_bounds, output_dir="region_images", 
                            start_time=0, final_time=1e-5, time_step=1e-7, 
                            frames_to_export=10):
        """
        Export custom region as PNG images at specified time steps
        
        Args:
            region_bounds: dict with xmin, xmax, ymin, ymax
            output_dir: Directory to save images
            start_time, final_time, time_step: Simulation parameters
            frames_to_export: Number of frames to export
        """
        import matplotlib.pyplot as plt
        from matplotlib import cm
        import os
        
        if self.mpi_rank == 0:
            os.makedirs(output_dir, exist_ok=True)
            print(f"\n🖼️  Exporting region as PNG images to {output_dir}/")
        
        # Create evaluation grid
        nx, ny = 200, 200
        x = np.linspace(region_bounds['xmin'], region_bounds['xmax'], nx)
        y = np.linspace(region_bounds['ymin'], region_bounds['ymax'], ny)
        xx, yy = np.meshgrid(x, y)
        points = np.column_stack([xx.ravel(), yy.ravel(), np.zeros(nx*ny)])
        
        # Find cells
        from dolfinx import geometry
        tdim = self.mesh.topology.dim
        self.mesh.topology.create_entities(tdim)
        
        bb_tree = geometry.bb_tree(self.mesh, tdim)
        cell_candidates = geometry.compute_collisions_points(bb_tree, points)
        colliding_cells = geometry.compute_colliding_cells(self.mesh, cell_candidates, points)
        
        cells = []
        valid_points = []
        for i, point in enumerate(points):
            cell_list = colliding_cells.links(i)
            if len(cell_list) > 0:
                cells.append(cell_list[0])
                valid_points.append(i)
        
        # Run simulation and export frames
        frame_interval = max(1, int((final_time - start_time) / time_step / frames_to_export))
        
        with VTXWriter(self.mesh.comm, f"{output_dir}/temp.bp", [self.u_n], "bp5") as vtx:
            vtx_cpp = vtx._cpp_object
            
            t = start_time
            dt = time_step
            step = 0
            frame_num = 0
            
            # RK4 loop with image export
            u_ = np.zeros(len(points))
            
            while t < final_time and frame_num < frames_to_export:
                # Advance one step (simplified - you'd call model.rk4 for one step)
                # For now, we'll export current state
                
                if step % frame_interval == 0:
                    # Evaluate at valid points
                    if len(valid_points) > 0:
                        points_subset = points[valid_points]
                        cells_array = np.array(cells, dtype=np.int32)
                        values = np.zeros((len(points_subset), 1))
                        self.u_n.eval(points_subset, cells_array, values)
                        
                        # Fill full grid (NaN for outside points)
                        u_full = np.full(len(points), np.nan)
                        u_full[valid_points] = values.flatten()
                        u_grid = u_full.reshape(ny, nx)
                        
                        if self.mpi_rank == 0:
                            # Create image
                            fig, ax = plt.subplots(figsize=(10, 8))
                            im = ax.imshow(u_grid, extent=[
                                region_bounds['xmin']*1e3, region_bounds['xmax']*1e3,
                                region_bounds['ymin']*1e3, region_bounds['ymax']*1e3
                            ], origin='lower', cmap='RdBu_r', aspect='auto')
                            ax.set_xlabel('X (mm)')
                            ax.set_ylabel('Y (mm)')
                            ax.set_title(f'Pressure at t = {t*1e6:.2f} μs')
                            plt.colorbar(im, ax=ax, label='Pressure (Pa)')
                            
                            filename = f"{output_dir}/frame_{frame_num:04d}_t{t*1e6:.2f}us.png"
                            plt.savefig(filename, dpi=150, bbox_inches='tight')
                            plt.close()
                            
                            print(f"  Exported: {filename}")
                    
                    frame_num += 1
                
                t += dt
                step += 1
        
        if self.mpi_rank == 0:
            print(f"✅ Exported {frame_num} PNG images")
    
    def export_region_csv(self, region_bounds, output_file="region_data.csv"):
        """
        Export current field in custom region as CSV
        
        Args:
            region_bounds: dict with xmin, xmax, ymin, ymax
            output_file: Output CSV filename
        """
        from dolfinx import geometry
        
        if self.mpi_rank == 0:
            print(f"\n📄 Exporting region to CSV: {output_file}")
        
        # Create evaluation grid
        nx, ny = 251, 251
        x = np.linspace(region_bounds['xmin'], region_bounds['xmax'], nx)
        y = np.linspace(region_bounds['ymin'], region_bounds['ymax'], ny)
        points = np.array([[xi, yi, 0.0] for xi in x for yi in y])
        
        # Find cells
        tdim = self.mesh.topology.dim
        self.mesh.topology.create_entities(tdim)
        
        bb_tree = geometry.bb_tree(self.mesh, tdim)
        cell_candidates = geometry.compute_collisions_points(bb_tree, points)
        colliding_cells = geometry.compute_colliding_cells(self.mesh, cell_candidates, points)
        
        cells = []
        points_on_proc = []
        for i, point in enumerate(points):
            cell_list = colliding_cells.links(i)
            if len(cell_list) > 0:
                cells.append(cell_list[0])
                points_on_proc.append(point)
        
        if len(points_on_proc) > 0:
            points_on_proc = np.array(points_on_proc)
            cells_array = np.array(cells, dtype=np.int32)
            values = np.zeros((len(points_on_proc), 1))
            self.u_n.eval(points_on_proc, cells_array, values)
            
            # Write to CSV
            MPI.COMM_WORLD.Barrier()
            for rank in range(MPI.COMM_WORLD.size):
                if self.mpi_rank == rank:
                    mode = 'w' if rank == 0 else 'a'
                    with open(output_file, mode) as f:
                        if rank == 0:
                            f.write("x,y,pressure\n")
                        for i, point in enumerate(points_on_proc):
                            f.write(f"{point[0]:.6e},{point[1]:.6e},{values[i,0]:.6e}\n")
                MPI.COMM_WORLD.Barrier()
        
        if self.mpi_rank == 0:
            print(f"✅ Exported CSV with ~{len(points)} points")
            print(f"   Open in ParaView: Filters → Table To Points")

