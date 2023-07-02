# Stream_l0_Vs_l1.py python script to superimpose streamlines of 2 leaflets using the matplotlib streamplot function
# Made by Dr Tyler Reddy and Dr Matthieu Chavent
# Oxford University, March 2014
# if you use this script for a scientific paper please cite the associated Faraday Discussions paper:
# Methodologies for the Analysis of Instantaneous Lipid Diffusion in MD Simulations of Large Membrane Systems
# http://pubs.rsc.org/en/content/articlelanding/2014/fd/c3fd00145h#!divAbstract

import string, os
import matplotlib
import MDAnalysis, MDAnalysis.visualization.streamlines
import matplotlib.pyplot, sys,numpy
import matplotlib.gridspec 
#import cPickle as pickle
import _pickle as cP

# define the different arguments of the script (optional or not) -------------
input_prot_gro = str(sys.argv[1]) 			# gro file for the proteins  
input_prot_xtc = str(sys.argv[2]) 			# xtc file for the proteins
input_leaflet1_gro = str(sys.argv[3])		# gro file for the lipids in leaflet1
input_leaflet1_xtc = str(sys.argv[4])		# xtc file for the lipids in leaflet1
input_leaflet2_gro = str(sys.argv[5])		# gro file for the lipids in leaflet2
input_leaflet2_xtc = str(sys.argv[6])		# xtc file for the lipids in leaflet2
input_selection = str(sys.argv[7])			# selection for the MDAnalysis by default  'name PO4' (see https://mdanalysis.googlecode.com/git/package/doc/html/documentation_pages/selections.html )
grid = int(sys.argv[8])						# size of the grid - ie resolution
xmin_in = int(sys.argv[9])					# xmin: in Angstrom
xmax_in = int(sys.argv[10])					# xmax: in Angstrom
ymin_in = int(sys.argv[11])					# ymin: in Angstrom
ymax_in = int(sys.argv[12])					# ymax: in Angstrom
startf = int(sys.argv[13])					# first frame
endf = int(sys.argv[14])					# last frame
dt = int(sys.argv[15])						# choose the time step, in general 1
filename = str(sys.argv[16])				# output filename for the picture
extension = str(sys.argv[17])				# extension of the picture (.png, .svg, etc ..), need to be managed by Matplotlib
plot_prot = str(sys.argv[18])				# True or Flase: to display proteins on the plot		
n_cores = int(sys.argv[19])					# number of cores you want to use

# for the example files just use this line: 
# python stream_l0_Vs_l1.py prot.gro prot_800_810.xtc leaflet0.gro l0_800_810.xtc leaflet1.gro l1_800_810.xtc "name PO4" 20 0 1160 0 1160 9 10 1 stream_test .png false 4
#---------------------------------------------------------------------------

if plot_prot == 'true': 
	# load protein traj in MDAnalysis
    U = MDAnalysis.Universe(input_prot_gro,input_prot_xtc)
    protein = U.select_atoms('protein')

# load lipid traj in MDAnalysis	
U2 = MDAnalysis.Universe(input_leaflet1_gro,input_leaflet1_xtc)
PO4_1 = U2.select_atoms('all')

U3 = MDAnalysis.Universe(input_leaflet2_gro,input_leaflet2_xtc)
PO4_2 = U3.select_atoms('all')

# calculate number of cells in x and y
nb_squares_x = (xmax_in-xmin_in)/grid - 1
nb_squares_y = (ymax_in-ymin_in)/grid - 1
salto_ticks=20
# loop on the trajectory frames
for i in range(startf, endf,dt):

	# define the next frame to calculate distamnce with.
    next_frame = i+dt
	
	# calculate displacement following x (u1) for each vector, displacement following y (v1) for each vector, average of the displacementand standard deviation of the displacement for U2. More info, http://mdanalysis.googlecode.com/git-history/develop/package/doc/html/documentation_pages/visualization/streamlines.html
    u1, v1, average_displacement1,standard_deviation_of_displacement1 = MDAnalysis.visualization.streamlines.generate_streamlines(input_leaflet1_gro,input_leaflet1_xtc, grid_spacing = grid, MDA_selection =input_selection, start_frame=startf, end_frame=next_frame, xmin = xmin_in, xmax = xmax_in, ymin = ymin_in , ymax = ymax_in, maximum_delta_magnitude = 5.0, num_cores=n_cores)

	# calculate displacement following x (u2) for each vector, displacement following y (v2) for each vector, average of the displacementand standard deviation of the displacement for U3. More info, http://mdanalysis.googlecode.com/git-history/develop/package/doc/html/documentation_pages/visualization/streamlines.html
    u2, v2, average_displacement2,standard_deviation_of_displacement2 = MDAnalysis.visualization.streamlines.generate_streamlines(input_leaflet2_gro,input_leaflet2_xtc, grid_spacing = grid, MDA_selection =input_selection, start_frame=startf, end_frame=next_frame, xmin = xmin_in, xmax = xmax_in, ymin = ymin_in , ymax = ymax_in, maximum_delta_magnitude = 5.0, num_cores=n_cores)

	# retrieve the data in the trajectory and store them in an intermediate file called streamline_test.p --------------------------------------------	

    if plot_prot == 'true': 
		# store the data from protein trajectory
        for ts in U.trajectory:  
            if ts.frame == i:
                tiempo=U.trajectory.time
                protein_coords = protein.positions 
                cP.dump([u1,v1,u2,v2,average_displacement1,standard_deviation_of_displacement1,average_displacement2,standard_deviation_of_displacement2, protein_coords],open('streamline_store.p','wb')) 
                break
    else:
        cP.dump([u1,v1,u2,v2,average_displacement1,standard_deviation_of_displacement1,average_displacement2,standard_deviation_of_displacement2],open('streamline_store.p','wb')) 

	# end of the retrieving part ------------------------------------------------------------------------------------------------------------------------

    nb_squares_x=int(round(nb_squares_x))
    nb_squares_y=int(round(nb_squares_y))
	# return evenly spaced numbers over a specified interval for x and y
    x = numpy.linspace(xmin_in,xmax_in,nb_squares_x)
    y = numpy.linspace(ymin_in,ymax_in,nb_squares_y)
 

	# read the data stored in the file streamline_test.p
    pickled_data = cP.load(open('streamline_store.p','rb'))
    if plot_prot == 'true': 		
        u1,v1,u2,v2,average_displacement1,standard_deviation_of_displacement1,average_displacement2,standard_deviation_of_displacement2, protein_coords = pickled_data #reconstitute coordinates (with the protein coordinates added) 
    else:
        u1,v1,u2,v2,average_displacement1,standard_deviation_of_displacement1,average_displacement2,standard_deviation_of_displacement2 = pickled_data #reconstitute coordinates (without the protein coordinates) 

	# print few info about average displacement and standard deviation 
    print( "------")
    print( "avg_disp1: ",str(average_displacement1))
    print("std_dev1: ",str(standard_deviation_of_displacement1))
    print("------")
    print("avg_disp2: ",str(average_displacement2))
    print("std_dev2: ",str(standard_deviation_of_displacement2))
    print("------")

	# calculate the speed for each vector
    speed1 = numpy.sqrt(u1*u1 + v1*v1)
    speed2 = numpy.sqrt(u2*u2 + v2*v2)	

	# define matplotlib figure and axes
    fig = matplotlib.pyplot.figure()


	# use the matplotlib streamplot (more info here: http://matplotlib.org/examples/images_contours_and_fields/streamplot_demo_features.html)
	# ie the density can be changed: increasing the number will result in increasing the number of streamlines

	# the big graph --------------
	# define position of the graph and create a squarred graph
    gs1 = matplotlib.gridspec.GridSpec(2, 2)
    gs1.update(left=0.1, right=0.58, wspace=0.02)
    ax = matplotlib.pyplot.subplot(gs1[:,:], aspect='equal');ax.set_title(f't: {tiempo/1000} ns ')
	# calculate the streamlines using matplotlib streamplot
    s1= ax.streamplot(x,y,u2,v2,density=(8,8),color=speed2,linewidth=3.0*speed2/speed2.max(),cmap=matplotlib.pyplot.cm.PuBu)
    s2= ax.streamplot(x,y,u1,v1,density=(8,8),color=speed1,linewidth=3.0*speed1/speed1.max(),cmap=matplotlib.pyplot.cm.OrRd)
	# set the transparency of the 2nd plot to better see the one in background	
    s2.lines.set_alpha(0.8)	
	# define the labels of the graph
    ax.set_ylabel('X [A]', size=15)
    ax.set_xlabel('Y [A]', size=15)
	# define the boundaries of the graph          
    ax.set_xlim([xmin_in,xmax_in])		
    ax.set_ylim([ymin_in,ymax_in])
    ax.set_xticks(numpy.arange(xmin_in,xmax_in,salto_ticks))
    ax.set_yticks(numpy.arange(ymin_in,ymax_in,salto_ticks))
	# end big graph --------------

	# the 2 small graphs ----------------------------
	# define position of the graph and create a squarred graph
    gs2 = matplotlib.gridspec.GridSpec(4, 2)
    gs2.update(left=0.63, right=0.95, hspace=1.2,wspace=0.1)

	# small graph1	
    ax1 = matplotlib.pyplot.subplot(gs2[:2, :],aspect='equal')
	# calculate the streamlines using matplotlib streamplot
    im1 = ax1.streamplot(x,y,u1,v1,density=(8,8),color=speed1,linewidth=3.0*speed1/speed1.max(),cmap=matplotlib.pyplot.cm.OrRd)
	# define the colorbar for graph1
    cbar1 = fig.colorbar(im1.lines,ax=ax1,orientation='vertical',format="%.1f")
    cbar1.set_label('Displacement [A] leaflet1',size=10)
	# define the labels of the graph
    ax1.tick_params(axis='both', labelsize=10)
	# define the boundaries of the graph          
    ax1.set_xlim([xmin_in,xmax_in])		
    ax1.set_ylim([ymin_in,ymax_in])
    ax1.set_xticks(numpy.arange(xmin_in,xmax_in,salto_ticks))
    ax1.set_yticks(numpy.arange(ymin_in,ymax_in,salto_ticks))
	
	#small graph2
    ax2 = matplotlib.pyplot.subplot(gs2[-2:, :],aspect='equal')
	# calculate the streamlines using matplotlib streamplot
    im2 = ax2.streamplot(x,y,u2,v2,density=(8,8),color=speed2,linewidth=3.0*speed2/speed2.max(),cmap=matplotlib.pyplot.cm.PuBu)
	# define the colorbar for graph2	
    cbar2 = fig.colorbar(im2.lines,ax=ax2,orientation='vertical',format="%.1f")
    cbar2.set_label('Displacement [A] leaflet2',size=10)
	# define the labels of the graph
    ax2.tick_params(axis='both', labelsize=10)
	# define the boundaries of the graph 
    ax2.set_xlim([xmin_in,xmax_in])		
    ax2.set_ylim([ymin_in,ymax_in])
    ax2.set_xticks(numpy.arange(xmin_in,xmax_in,salto_ticks))
    ax2.set_yticks(numpy.arange(ymin_in,ymax_in,salto_ticks))
	# end 2 small graphs --------------------------


	# plot the protein coordinates with small purple circles
    if plot_prot == 'true':
        prot_x = protein_coords[...,0]
        prot_y = protein_coords[...,1]
        ax.plot(prot_x,prot_y,'o',markerfacecolor='#660066',markeredgecolor='#660066',markersize=0.8,alpha=0.5)
        ax1.plot(prot_x,prot_y,'o',markerfacecolor='#660066',markeredgecolor='#660066',markersize=0.8,alpha=0.5)
        ax2.plot(prot_x,prot_y,'o',markerfacecolor='#660066',markeredgecolor='#660066',markersize=0.8,alpha=0.5)
	# definition of the figure name: 
    aligned_i = "%05d" % i # useful to perform a movie
    figure_name = filename+aligned_i+extension
    fig.savefig(figure_name,dpi=200)



