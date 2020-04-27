import numpy as np

"""
Post processing files
"""
infile="Averages_57_avg.fld"
casefile="base/pin.xml"
meshfile="base/mesh.xml"




"""
Input geometry
"""

BL_arc = False

#pin diameter
D = 1.0

#number of fins
N_x = 6 #streamwise
N_y = 5 

#center to center distances
S_y = 3*D
S_x = np.sqrt(3.0)*0.5*S_y

#height of the channel
H = 18*D

#lengths before and after the fin array
L_before = 2*H
L_after = 2.5*H

#cylinder boundary layer height
#R_bl = 1.5*(0.5*D) #if arcs are used

R_bl = 2.0*(0.5*D) #if arcs are used

#height of the channel boundary layer mesh
bl_height = 0.1*D


"""
Cell counts
"""
#number of cells in the cylinder vicinity
Ncells_theta = 5 #has to be almost 15 if first order mesh is used
Ncells_bl = 4
grading_bl = 1.3

#number of cells between the cylinders
Ncells_sy = 5
Ncells_sx = 4

#number of cells before and after the fin arrays (x-dir)
Ncells_x_before = 15
Ncells_x_after = 19

grading_x_before = 1.13 #1.05
grading_x_after =  1.0/1.13

#number of cells in the channel boundary layer
Ncells_ch_bl = 3
#number of cells y-dir in the next block of channel bl
Ncells_ch_bl_next = 8
#number of cells x-dir immediately after the fin array
Ncells_x_array_after = Ncells_sx + 2


#base size for unstructured mesh NOT USED
base_size = D/1.5





