#!/usr/bin/env python3
"""
compute_z_density_cg_centered.py

Compute and plot 1D density profiles along Z (nm), centered on the membrane COM,
normalized to 1, for CG ions (CA+), peptide backbone beads (BB), and membrane.

Usage:
    python compute_z_density_cg_centered.py topology.pdb traj.xtc

Dependencies:
    - MDAnalysis
    - NumPy
    - Matplotlib
"""

import sys
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt

def compute_z_profile_centered(u, sel_str, mem_group, nbins=200):
    """
    Compute a centered 1D density histogram along Z relative to membrane COM.

    Parameters
    ----------
    u : MDAnalysis.Universe
        Universe with topology + trajectory.
    sel_str : str
        Selection for target atoms.
    mem_group : AtomGroup
        AtomGroup for membrane to compute COM.
    nbins : int
        Number of bins along Z.

    Returns
    -------
    centers : np.ndarray
        Bin centers (nm) relative to membrane COM.
    H_norm : np.ndarray
        Normalized density profile.
    """
    ag = u.select_atoms(sel_str)
    all_z_rel = []

    for ts in u.trajectory:
        # compute membrane COM along Z (Å)
        com_mem_z = mem_group.center_of_mass()[2]
        # get Z positions of selection and subtract COM
        zs = ag.positions[:, 2] - com_mem_z
        all_z_rel.append(zs)

    all_z_rel = np.hstack(all_z_rel)

    # define histogram range symmetric around zero
    zmax = np.max(np.abs(all_z_rel))
    edges = np.linspace(-zmax, zmax, nbins + 1)

    H, _ = np.histogram(all_z_rel, bins=edges, density=False)
    H_norm = H / H.max()

    # centers in nm
    centers = (edges[:-1] + edges[1:]) / 2.0 / 10.0
    return centers, H_norm


def main(pdb, xtc, nbins=200):
    u = mda.Universe(pdb, xtc)

    # define selections
    ion_sel = "resname ION and name CA+ or name NA+"       # CG ions
    bb_sel  = "name BB"                         # peptide BB beads
    mem_sel = "resname POPC POPE and name PO4"  # membrane beads

    # preselect membrane group
    mem_group = u.select_atoms(mem_sel)

    # compute centered profiles
    try:
        z_ions, prof_ions = compute_z_profile_centered(u, ion_sel, mem_group, nbins)
    except:
        z_ions=0;prof_ions = 0
    z_bb,   prof_bb   = compute_z_profile_centered(u, bb_sel,  mem_group, nbins)
    z_mem,  prof_mem  = compute_z_profile_centered(u, mem_sel, mem_group, nbins)

    # plot
    plt.figure(figsize=(7,4))
    plt.plot(z_ions, prof_ions, label='Positive ions', linewidth=1.5)
    plt.plot(z_bb,   prof_bb,   label='Peptide BB', linewidth=1.5)
    plt.plot(z_mem,  prof_mem,  label='Membrane P', linewidth=1.5)
    plt.xlabel('Z-position (nm)',fontsize=20)
    plt.ylabel('Normalized density',fontsize=20)
    plt.legend(frameon=False)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.axvline(0, color='k', linestyle='--', linewidth=1)
    plt.tight_layout()
    plt.savefig("z_density_profiles_cg_centered.png", dpi=300)
    print("Saved: z_density_profiles_cg_centered.png")

    # save data
    np.savez("z_density_profiles_cg_centered.npz",
             z=z_ions,
             ions=prof_ions,
             bb=prof_bb,
             mem=prof_mem)
    print("Saved: z_density_profiles_cg_centered.npz")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compute_z_density_cg_centered.py topology.pdb traj.xtc")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

