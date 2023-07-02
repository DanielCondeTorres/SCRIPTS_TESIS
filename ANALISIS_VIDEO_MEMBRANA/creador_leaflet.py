import sys, string, os, subprocess
import numpy as np
import MDAnalysis
import MDAnalysis.analysis
from   MDAnalysis.analysis.leaflet import LeafletFinder
import MDAnalysis.selections.gromacs

# gro file in input
input_filename = str(sys.argv[1])
# selection to use to divide the 2 leaflet in general headgroup is the best
selection = str(sys.argv[2])
# name for the output leaflet1 no need of extension, the gro extension will be added automatically
output_filename = str(sys.argv[3])
# name for the output leaflet2 no need of extension, the gro extension will be added automatically
output_filename2 = str(sys.argv[4])

# store corrdinates, etc in the universe (U)
U = MDAnalysis.Universe(input_filename)
print('UNIVERSO:::::::::::',U)
# use the selection to define the 2 leaflets
L = LeafletFinder(U,selection)

# store the info of each leaflet in a variable
leaflet0 = L.groups(0)
leaflet1 = L.groups(1)

# create the gro files
W0_gro = MDAnalysis.coordinates.GRO.GROWriter(output_filename+".gro")
W0_gro.write(leaflet0.residues.atoms)

W1_gro = MDAnalysis.coordinates.GRO.GROWriter(output_filename2+".gro")
W1_gro.write(leaflet1.residues.atoms)

# create the ndx files - useful to separate the leaflets in an xtc file
W0_ndx = MDAnalysis.selections.gromacs.SelectionWriter(output_filename+".ndx")
W0_ndx.write(leaflet0.residues.atoms)

W1_ndx = MDAnalysis.selections.gromacs.SelectionWriter(output_filename2+".ndx")
W1_ndx.write(leaflet1.residues.atoms)
