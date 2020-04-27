#IHI L = 6., B = 0.6, D=0.375
#SRI L = 4., B = 0.4, D=0.25
#UT  L = 2.5 B = 0.25 D=0.156
#YNU L=2.0,	B=0.25, D=0.125 
#B IS THE FULL HULL WIDTH!!!!!!!

##############################################
#SHIP DIMENSIONS
##############################################
L = 2.5			#Ship length, x ranges from -0.5*L to 0.5*L
T = 0.156 		#Depth below waterline
B = 0.25*0.5 		#Breadth y ranges from 0 to B such that half hull is created


##############################################
#DOMAIN DIMENSIONS
##############################################
Lx_before_hull = 2 * L
Lx_after_hull = 4.0 * L
Lz_over_hull = 1 * L
Lz_below_hull = 2.0 * L
Ly = 2.0 * L
Ly_boundary_layer = 0.003 * L #this is the thickness of the boundary layer mesh


##############################################
#CELL COUNTS
##############################################
multiplier = 1.0
Nx_hull = 80*multiplier		#total amount of cells per ship length
Nx_before_hull = 20*multiplier
Nx_after_hull = 80*multiplier

Nz_over_hull = 20*multiplier	#this is the cell count in air
Nz_hull = 30*multiplier
Nz_below_hull = 20*multiplier
Ny = 110*multiplier
Ny_boundary_layer = 5 #number of boundary layer cells 


##############################################
#GRADING
##############################################
SetOptimalGrading = True #if this is true all but Gy and Gz_hull are overwritten
Gx_before_hull = 5.0
Gx_after_hull = 7.816749006
Gz_over_hull = 20.806
Gz_below_hull = 50.

Gz_hull = 1. 





n_blocks_x = 10 #number of longitudal blocks per ship length


outputfile = "system/blockMeshDict"
#outputfile = "constant/polyMesh/blockMeshDict"

num_spline_points = 10

