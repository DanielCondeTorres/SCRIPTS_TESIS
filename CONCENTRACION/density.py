import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.density import DensityAnalysis
import matplotlib
matplotlib.use('Agg')  # Disable graphical interface
import matplotlib.pyplot as plt
import argparse

def normalize_density(density):
    """Normalizes the density by dividing by the maximum value."""
    return density / np.max(density) if np.max(density) > 0 else density

def plot_density_map(pdb_file, xtc_file, output_dir, plane_type):
    # Load the trajectory
    u = mda.Universe(pdb_file, xtc_file)

    # Get the simulation box dimensions from the first frame
    box_dimensions = u.trajectory[0].dimensions[:3]  # Get the first frame box dimensions (X, Y, Z)
    print(f"Box dimensions: {box_dimensions}")
    # Get the simulation dt (in ps)
    dt = u.trajectory.dt  # Time step between frames (ps)

    # Calculate the frame at which 10 μs (10,000 ps) is reached
    frame_final = int(10000 / dt)

    # Select Ca+ ions, membrane, and BB atoms
    ca_ions = u.select_atoms("name CA+ or name NA+")
    membrane = u.select_atoms("resname POPC POPS")
    bb_atoms = u.select_atoms("name BB")

    # Calculate the center of mass of the membrane
    mem_COM = membrane.center_of_mass()

    # Configure density analysis using the box dimensions
    density_ca = DensityAnalysis(ca_ions, delta=1.0, gridcenter=mem_COM, xdim=box_dimensions[0], ydim=box_dimensions[1], zdim=box_dimensions[2])
    density_ca.run(start=0, stop=len(u.trajectory), step=100)

    density_membrane = DensityAnalysis(membrane, delta=1.0, gridcenter=mem_COM, xdim=box_dimensions[0], ydim=box_dimensions[1], zdim=box_dimensions[2])
    density_membrane.run(start=0, stop=len(u.trajectory), step=100)

    density_bb = DensityAnalysis(bb_atoms, delta=1.0, gridcenter=mem_COM, xdim=box_dimensions[0], ydim=box_dimensions[1], zdim=box_dimensions[2])
    density_bb.run(start=0,stop=len(u.trajectory), step=100)

    # Extract and normalize densities for each plane
    def generate_density(plane_type):
        if plane_type == "XY":
            return normalize_density(density_ca.results.density.grid.mean(axis=2)), normalize_density(density_membrane.results.density.grid.mean(axis=2)), normalize_density(density_bb.results.density.grid.mean(axis=2))
        elif plane_type == "XZ":
            return normalize_density(density_ca.results.density.grid.mean(axis=1)), normalize_density(density_membrane.results.density.grid.mean(axis=1)), normalize_density(density_bb.results.density.grid.mean(axis=1))
        elif plane_type == "YZ":
            return normalize_density(density_ca.results.density.grid.mean(axis=0)), normalize_density(density_membrane.results.density.grid.mean(axis=0)), normalize_density(density_bb.results.density.grid.mean(axis=0))
        else:
            raise ValueError("Invalid plane type. Choose 'XY', 'XZ', or 'YZ'.")

    # If 'all' is selected, generate maps for all planes
    if plane_type == "all":
        planes = ["XY", "XZ", "YZ"]
    else:
        planes = [plane_type]

    # Save images for each plane
    for plane in planes:
        ca_grid, mem_grid, bb_grid = generate_density(plane)
        title_size = 18
        label_size = 16
        ticks_size = 14
        bar_size = 14
        # Save image for Ca+ in the selected plane
        plt.figure(figsize=(8, 6))
        im = plt.imshow(ca_grid.T, origin='lower', cmap='Purples', alpha=0.6, vmin=0, vmax=1)
        cbar = plt.colorbar(im)
        # Ajustar el tamaño de la etiqueta del colorbar
        cbar.set_label("Normalized Density", fontsize=bar_size)
        # Ajustar el tamaño de los ticks del colorbar
        cbar.ax.tick_params(labelsize=ticks_size-2)
        plt.title(f"Ca+2 Density Map (Normalized) - {plane} Plane",fontsize=title_size)
        plt.xlabel(f"{plane[0]} (Å)",fontsize=label_size)  # X-axis (or relevant axis)
        plt.ylabel(f"{plane[1]} (Å)",fontsize=label_size)  # Y-axis (or relevant axis)
        plt.tick_params(axis='both', which='major', labelsize=ticks_size)
        plt.savefig(f"{output_dir}/ca_{plane.lower()}.png", dpi=300, bbox_inches='tight')
        plt.close()

        # Save image for Membrane in the selected plane
        plt.figure(figsize=(8, 6))
        im=plt.imshow(mem_grid.T, origin='lower', cmap='YlOrBr', alpha=0.6)
        cbar = plt.colorbar(im)

        # Ajustar el tamaño de la etiqueta del colorbar
        cbar.set_label("Normalized Density", fontsize=bar_size)

        # Ajustar el tamaño de los ticks del colorbar
        cbar.ax.tick_params(labelsize=ticks_size-2)
        plt.title(f"Membrane Density Map (Normalized) - {plane} Plane",fontsize=title_size)
        plt.xlabel(f"{plane[0]} (Å)",fontsize=label_size)  # X-axis (or relevant axis)
        plt.ylabel(f"{plane[1]} (Å)",fontsize=label_size)  # Y-axis (or relevant axis)
        plt.tick_params(axis='both', which='major', labelsize=ticks_size)
        plt.savefig(f"{output_dir}/membrane_{plane.lower()}.png", dpi=300, bbox_inches='tight')
        plt.close()

        # Save image for BB in the selected plane
        plt.figure(figsize=(8, 6))
        im=plt.imshow(bb_grid.T, origin='lower', cmap='Reds', alpha=0.6)
        cbar = plt.colorbar(im)

        # Ajustar el tamaño de la etiqueta del colorbar
        cbar.set_label("Normalized Density", fontsize=bar_size)

        # Ajustar el tamaño de los ticks del colorbar
        cbar.ax.tick_params(labelsize=ticks_size-2)
        plt.title(f"BB Density Map (Normalized) - {plane} Plane",fontsize=title_size)
        plt.xlabel(f"{plane[0]} (Å)",fontsize=label_size)  # X-axis (or relevant axis)
        plt.ylabel(f"{plane[1]} (Å)",fontsize=label_size)  # Y-axis (or relevant axis)
        plt.tick_params(axis='both', which='major', labelsize=ticks_size)
        plt.savefig(f"{output_dir}/bb_{plane.lower()}.png", dpi=300, bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate and save density maps for different components.")
    parser.add_argument("-pdb", type=str, required=True, help="PDB file of the simulation.")
    parser.add_argument("-xtc", type=str, required=True, help="XTC file of the simulation.")
    parser.add_argument("-o", "--output_dir", type=str, default=".", help="Output directory for images (default: current directory)")
    parser.add_argument("-plane", type=str, choices=["XY", "XZ", "YZ", "all"], required=True, help="Choose the plane for the density map ('XY', 'XZ', 'YZ', or 'all' for all planes).")
    print('Executing density.py script...')
    args = parser.parse_args()

    plot_density_map(args.pdb, args.xtc, args.output_dir, args.plane)
