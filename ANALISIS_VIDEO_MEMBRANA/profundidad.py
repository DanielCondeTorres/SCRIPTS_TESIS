# Stream_Vs_z_2D.py python script to to display streamlines using the matplotlib streamplot function
# and superimpose it with heightmap usinf matplotlib contourf function
# Made by Dr Tyler Reddy and Dr Matthieu Chavent
# Oxford University, March 2014
# if you use this script for a scientific paper please cite the associated Faraday Discussions paper:
# Methodologies for the Analysis of Instantaneous Lipid Diffusion in MD Simulations of Large Membrane Systems
# http://pubs.rsc.org/en/content/articlelanding/2014/fd/c3fd00145h#!divAbstract

import string, os
import matplotlib
import MDAnalysis, MDAnalysis.visualization.streamlines
import matplotlib.pyplot, sys,numpy
import _pickle as cP
import numpy as np
# define the different arguments of the script (optional or not) -------------
input_prot_gro = str(sys.argv[1]) 			# gro file for the proteins  
input_prot_xtc = str(sys.argv[2]) 			# xtc file for the proteins
input_leaflet_gro = str(sys.argv[3])		# gro file for the lipids
input_leaflet_xtc = str(sys.argv[4])		# xtc file for the lipids
input_selection = str(sys.argv[5])			# selection for the MDAnalysis by default  'name PO4' (see https://mdanalysis.googlecode.com/git/package/doc/html/documentation_pages/selections.html )
grid = int(sys.argv[6])						# size of the grid - ie resolution
xmin_in = int(sys.argv[7])					# xmin: in Angstrom
xmax_in = int(sys.argv[8])					# xmax: in Angstrom
ymin_in = int(sys.argv[9])					# ymin: in Angstrom
ymax_in = int(sys.argv[10])					# ymax: in Angstrom
zmin_in = int(sys.argv[11])					# zmin: in Angstrom, give the minimum value for the colorbar
zmax_in = int(sys.argv[12])					# zmax: in Angstrom, give the maximum value for the colorbar
startf = int(sys.argv[13])					# first frame
endf = int(sys.argv[14])					# last frame
dt = int(sys.argv[15])						# choose the time step, in general 1
filename = str(sys.argv[16])				# output filename for the picture
extension = str(sys.argv[17])				# extension of the picture (.png, .svg, etc ..), need to be managed by Matplotlib
plot_prot = str(sys.argv[18])				# True or Flase: to display proteins on the plot		
n_cores = int(sys.argv[19])					# number of cores you want to use

# for the example files just use this line: 
# python stream_Vs_z_2D.py prot.gro prot_800_810.xtc leaflet0.gro l0_800_810.xtc "name PO4" 20 0 1160 0 1160 30 90 9 10 1 test .png false 4
#---------------------------------------------------------------------------

if plot_prot == 'true': 
	# load protein traj in MDAnalysis
	U = MDAnalysis.Universe(input_prot_gro,input_prot_xtc)
	protein = U.select_atoms('protein')

# load lipid traj in MDAnalysis	
U2 = MDAnalysis.Universe(input_leaflet_gro,input_leaflet_xtc)
PO4 = U2.select_atoms('all')

# calculate number of cells in x and y
nb_squares_x = (xmax_in-xmin_in)/grid - 1
nb_squares_y = (ymax_in-ymin_in)/grid - 1

# loop on the trajectory frames
for i in range(startf, endf,dt):
     
	# define the next frame to calculate distamnce with.
	next_frame = i+dt;u1, v1, average_displacement,standard_deviation_of_displacement = MDAnalysis.visualization.streamlines.generate_streamlines(input_leaflet_gro,input_leaflet_xtc, grid_spacing = grid, MDA_selection =input_selection, start_frame=startf, end_frame=next_frame, xmin = xmin_in, xmax = xmax_in, ymin = ymin_in , ymax = ymax_in, maximum_delta_magnitude = 5.0, num_cores=n_cores)


	# retrieve the data in the trajectory and store them in an intermediate file called streamline_test.p --------------------------------------------	

	# store the data from lipid trajectory
	for ts in U2.trajectory:  
		if ts.frame == i:
			PO4_coords = PO4.positions;tiempo=U2.trajectory.time
			cP.dump([u1,v1,average_displacement,standard_deviation_of_displacement, PO4_coords],open('streamline_store.p','wb'))
			break

	if plot_prot == 'true': 
		# store the data from protein trajectory
		for ts in U.trajectory:  
			if ts.frame == i:
				protein_coords = protein.positions #apres code multicore
				cP.dump([u1,v1,average_displacement,standard_deviation_of_displacement, PO4_coords, protein_coords],open('streamline_store.p','wb')) 
				break
	
	# end of the retrieving part ------------------------------------------------------------------------------------------------------------------------


	# return evenly spaced numbers over a specified interval for x and y
	print(xmin_in,xmax_in,nb_squares_x);x = numpy.linspace(xmin_in,xmax_in,round(nb_squares_x));nb_squares_x=round(nb_squares_x);nb_squares_y=round(nb_squares_y)
	y = numpy.linspace(ymin_in,ymax_in,int(nb_squares_y));pickled_data = cP.load(open('streamline_store.p','rb'))
	if plot_prot == 'true': 		
		u1, v1, average_displacement, standard_deviation_of_displacement, PO4_coords, protein_coords = pickled_data #reconstitute coordinates (with the protein coordinates added) 
	else:
		u1, v1, average_displacement, standard_deviation_of_displacement, PO4_coords = pickled_data #reconstitute coordinates (without the protein coordinates) 


	# calculate the speed for each vector
	speed = numpy.sqrt(u1*u1 + v1*v1)
	
	# define matplotlib figure and axes
	fig = matplotlib.pyplot.figure()
	ax = fig.add_subplot(111)

	# label matplotlib axes	
	matplotlib.pyplot.xlabel('X [A]', size=15)
	matplotlib.pyplot.ylabel('Y [A]', size=15)

	# fill the matrices by zeros
	population=numpy.zeros((nb_squares_x,nb_squares_y))
	z_avg_tmp=numpy.zeros((nb_squares_x,nb_squares_y))

	# fill population matrix
	population +=  numpy.histogram2d(PO4_coords[:,0], PO4_coords[:,1], bins=(nb_squares_x,nb_squares_y), range=[[xmin_in, xmax_in], [ymin_in, ymax_in]])[0]
	# fill z_avg matrix
	z_avg_tmp +=  numpy.histogram2d(PO4_coords[:,0], PO4_coords[:,1],bins=(nb_squares_x,nb_squares_y), range=[[xmin_in, xmax_in], [ymin_in, ymax_in]], weights= PO4_coords[:,2])[0]

 	#this create an array with "1" where there's no thickness value, which allows the division below
	pop = population < 1  
	z_avg_tmp = z_avg_tmp/(population + pop)

	# define the colour surface using matplotlib countourf function (see: http://matplotlib.org/examples/pylab_examples/contourf_demo.html)
	levels = matplotlib.ticker.MaxNLocator(nbins=100).tick_values(zmin_in, zmax_in)
	matplotlib.pyplot.contourf(x, y, numpy.transpose(z_avg_tmp), cmap=matplotlib.pyplot.cm.jet, levels=levels, vmin = zmin_in , vmax = zmax_in)	
	
	# colorbar for the z_value and label of the colorbar
	cbar = matplotlib.pyplot.colorbar()
	cbar.set_label('z_value [A]',size=15);print('al palo');print('SHAPES',np.shape(x),np.shape(y),np.shape(u1),np.shape(v1))

	# use the matplotlib streamplot (more info here: http://matplotlib.org/examples/images_contours_and_fields/streamplot_demo_features.html)
	# ie the density can be changed: increasing the number will result in increasing the number of streamlines
	im = ax.streamplot(x,y,u1,v1,density=(8,8),color=speed,linewidth=1.5*speed/speed.max(),cmap=matplotlib.pyplot.cm.bone);print('llega')	
	# plot the protein coordinates with small purple circles
	if plot_prot == 'true':
		prot_x = protein_coords[...,0]
		prot_y = protein_coords[...,1];ax.set_title(f't: {tiempo/1000} ns ')#ax.text(np.median(prot_x), max(prot_y) + 2, f't: {tiempo/1000} ns ', fontsize=12, ha='center')
		ax.plot(prot_x,prot_y,'o',markerfacecolor='#660066',markeredgecolor='#660066',markersize=0.8,alpha=0.5)
		
	
	# definition of axes for the graph
	matplotlib.pyplot.axis([xmin_in,xmax_in,ymin_in,ymax_in]);ax.set_xlim([xmin_in,xmax_in]);ax.set_ylim([ymin_in,ymax_in]);ax.set_xticks(numpy.arange(xmin_in,xmax_in,20));ax.set_yticks(numpy.arange(ymin_in,ymax_in,20))

	# definition of the figure name: 
	aligned_i = "%05d" % i # useful to perform a movie
	figure_name = filename+aligned_i+extension
	fig.savefig(figure_name,dpi=200)
