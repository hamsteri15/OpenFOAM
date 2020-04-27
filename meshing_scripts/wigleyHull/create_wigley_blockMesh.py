
import numpy as np
import params as p #params.py file should be located in the folder where you run this script (usually the case folder). 
from scipy.optimize import fsolve #to set grading

def get_grading(dx_start, ncells, total_L):
	
	func = lambda r : dx_start*(1.0 - r**(ncells))/(1.0 - r) - total_L
	
	init_guess = 1.4 
	cell_to_cell = fsolve(func, init_guess)
	grading = cell_to_cell**(ncells - 1)
	return grading[0]

L = p.L			
T = p.T 			
B = p.B 		
Lx_before_hull = p.Lx_before_hull
Lx_after_hull = p.Lx_after_hull
Lz_over_hull = p.Lz_over_hull
Lz_below_hull = p.Lz_below_hull
Ly = p.Ly
Nx_hull = int(round(p.Nx_hull))
Nx_before_hull = int(round(p.Nx_before_hull))
Nx_after_hull = int(round(p.Nx_after_hull))
Nz_over_hull = int(round(p.Nz_over_hull))
Nz_hull = int(round(p.Nz_hull))
Nz_below_hull = int(round(p.Nz_below_hull))
Ny = int(round(p.Ny))




Ly_boundary_layer = p.Ly_boundary_layer
Ny_boundary_layer = int(round(p.Ny_boundary_layer))


Gz_hull = p.Gz_hull



if (p.SetOptimalGrading):
	dx_hull = (L/Nx_hull)#*0.5
	Gx_before_hull = get_grading(dx_hull, Nx_before_hull, Lx_before_hull) #p.Gx_before_hull
	Gx_after_hull = get_grading(dx_hull, Nx_after_hull, Lx_after_hull)
	
	dz_hull = T/Nz_hull
	Gz_over_hull = get_grading(dz_hull, Nz_over_hull, Lz_over_hull)
	Gz_below_hull = get_grading(dz_hull, Nz_below_hull, Lz_below_hull)
	
	
	
else:

	Gx_before_hull = p.Gx_before_hull	
	Gx_after_hull = p.Gx_after_hull
	Gz_over_hull = p.Gz_over_hull
	Gz_below_hull = p.Gz_below_hull



n_blocks_x = p.n_blocks_x
outputfile = p.outputfile


num_spline_points = p.num_spline_points

def points_imax():
	return n_blocks_x + 1 + 2

#constants for now
def points_jmax():
	return 2

"""
def points_kmax():
	return 5
"""	
def points_kmax():
	return 4




class Writer(object):
	
	def __init__(self, o_file):
	
		self.f = open(o_file,"w+")
		


	def write_header(self):
		

		string = "/**/\n\
		FoamFile\n\
		{\n\
		\t version     2.0;\n\
		\t format      ascii;\n\
		\t class       dictionary;\n\
		\t object      blockMeshDict;\n\
		}\n"
			
		self.f.write(string + 10*'\n')
	
		string = "convertToMeters 1;"
		self.f.write(string)


	def write_points(self, points):
		
		self.f.write(5*"\n")
		self.f.write("vertices\n")
		self.f.write("(\n")
		imax = points_imax(); jmax = points_jmax(); kmax = points_kmax();
		
		for i in range(imax):
			for j in range(jmax):
				for k in range(kmax):
		
					string = "( {} {} {} )	//{}\n".format(points[i,j,k].x, points[i,j,k].y, points[i,j,k].z, points[i,j,k].n)
					self.f.write(string)
		self.f.write(");")
	
	def write_blocks(self, blocks):

		

		self.f.write(5*"\n")
		self.f.write("blocks\n")
		self.f.write("(\n")
	
		for block in blocks:
			string = block.get_string()
			self.f.write(string + "\n")
		
		self.f.write("\n);")			

	def write_edges(self, edges):
		self.f.write(5*"\n")
		self.f.write("edges\n")
		self.f.write("(\n")
		for edge in edges:
			self.f.write(edge.get_string() + "\n")
		self.f.write(");")

	def write_patches(self, patches):
		self.f.write(5*"\n")
		self.f.write("boundary\n")
		self.f.write("(\n")
	
		for patch in patches:
			self.f.write(patch.get_string())
			self.f.write(5*"\n")
	
	
		self.f.write(");\n")			
	
		self.f.write("mergePatchPairs ();")


class Block(object):
	
	
	
	def __init__(self, p1, p2, p3 , p4, p5, p6, p7, p8, num, identifier, nx=10, ny=10, nz=10):
		
		self.p1 = p1; self.p2 = p2; self.p3 = p3; self.p4 = p4;
	 	self.p5 = p5; self.p6 = p6; self.p7 = p7; self.p8 = p8;
		
		self.nx = nx; self.ny = ny; self.nz = nz;
		

		self.edge_grading = []
		

		self.n = num
		
		self.identifier = identifier

		
				
	def neg_ynormal_face(self):
		return [self.p4.n, self.p3.n, self.p7.n, self.p8.n]

	def pos_ynormal_face(self):
		return [self.p1.n, self.p2.n, self.p6.n, self.p5.n]
	
	def neg_xnormal_face(self):
		return [self.p1.n, self.p4.n, self.p8.n, self.p5.n]

	def pos_xnormal_face(self):
		return [self.p3.n, self.p2.n, self.p6.n, self.p7.n]

	def neg_znormal_face(self):
		return [self.p6.n, self.p7.n, self.p8.n, self.p5.n]




	def pos_znormal_face(self):
		return [self.p1.n, self.p4.n, self.p3.n, self.p2.n]


	def get_y_grading(self):
		"""
		sets the boundary layer cells and the grading after the BL
		to match the BL_width
		"""
	
		length_ratio = Ly_boundary_layer / Ly
		remaining_length = 1.0 - length_ratio 
	
		cell_ratio = float(Ny_boundary_layer) / Ny
		remaining_cells = 1.0 - cell_ratio
	
		dy_bl =  Ly_boundary_layer / Ny_boundary_layer
	
		
		remaining_cells_real = float(Ny - Ny_boundary_layer)
		remaining_length_real = Ly * remaining_length
		
		grading_after = get_grading(dy_bl, remaining_cells_real, remaining_length_real)

	
		temp = " ( ({} {} 1.0) ({} {} {}) ) ".format(length_ratio, cell_ratio, remaining_length, remaining_cells, grading_after)
		return temp * 4 #all 4 edges in y-direction have the same grading


	def get_z_grading_air(self):

		thickness = 0.25*T
		length_ratio = thickness / Lz_over_hull
		remaining_length = 1.0 - length_ratio
		dz_hull = T /Nz_hull
		ncells_close = round(thickness / dz_hull)

		cell_ratio = float(ncells_close) / Nz_over_hull
		remaining_cells = 1.0 - cell_ratio
		dy_bl =  Ly_boundary_layer / Ny_boundary_layer

		remaining_cells_real = Nz_over_hull - ncells_close #Nz_over_hull * remaining_cells
		remaining_length_real = Lz_over_hull * remaining_length
		grading_after = get_grading(dz_hull, remaining_cells_real, remaining_length_real)

		temp = " ( ({} {} 1.0) ({} {} {}) ) ".format(length_ratio, cell_ratio, remaining_length, remaining_cells, grading_after)
		return temp * 4 #all 4 edges in y-direction have the same grading

	def get_string(self):

		string = "hex ({} {} {} {} {} {} {} {})".format(self.p1.n, \
										 self.p2.n, \
										 self.p3.n, \
										 self.p4.n, \
										 self.p5.n, \
										 self.p6.n, \
										 self.p7.n, \
										 self.p8.n)
		
		string+=" ({} {} {})".format(self.nx, self.ny, self.nz)
		
		
		
		string+=" edgeGrading ( "
		
		#x grading
		for i in range(4):
			string += str(self.edge_grading[i]) + " "
		
		#y-grading
		string += self.get_y_grading()
		
		#z-grading
		if ("air" in self.identifier):
			string += self.get_z_grading_air()
		else:    
			for i in range(8, 12):
				string += str(self.edge_grading[i]) + " "
		
		
		
		 
		string+=")" 

		return string


class Patch(object):

	def __init__(self, name, type, block_faces):

		self.name = name
		self.type = type
		self.block_faces = block_faces
		

	def get_string(self):

		string = self.name + "\n {"
		string+= "	type {};\n".format(self.type)
		string+= "	faces \n"
		string+= "	(\n"
	
		for face in self.block_faces:
			string+="({} {} {} {})".format(face[0], face[1], face[2], face[3]) + "\n"
	
		string+="	);\n"
		string+=" }\n"
		return string
		


class Edge(object):
	
	def __init__(self, p1, p2, num):
		
		self.p1 = p1
		self.p2 = p2
		self.num = num

		self.spline_points = []
		
		self.create_spline_points()
	
	def create_spline_points(self):

		x_coords = np.linspace(self.p1.x, self.p2.x, num_spline_points)
		z_coords = np.linspace(self.p1.z, self.p2.z, num_spline_points)
		
		
		for idx in range(num_spline_points):
			y = hull_shape(x_coords[idx], z_coords[idx])	
			self.spline_points.append(Point(x_coords[idx], y, z_coords[idx], 100))
		
	def get_string(self):
	
		string = "spline {} {} ( ".format(self.p1.n, self.p2.n)
		for idx in range(num_spline_points):
			x = self.spline_points[idx].x
			y = self.spline_points[idx].y
			z = self.spline_points[idx].z
			string += "( " + str(x) + " " + str(y) + " " + str(z) + " ) \n"
		string += ")"
	
		return string


class Point(object):
	
	
	
	def __init__(self, x, y, z, num):
		
		self.x = x; self.y = y; self.z = z;
		self.n = num
	
		



def main():

	create_mesh()
	




def create_hull_points():
	
	#x_coords = symmetry_x_coords()
	x_coords = create_x_coords()
	z_coords = all_z_coords()
	
	points = []
	idx = 0
	for i in range(points_imax()):
		for k in range(points_kmax()):
					
			
			x = x_coords[i]
			z = z_coords[k]
			y = hull_shape(x,z)

			
			#points on hull surface and atmosphere
			if (x < 0.5*L and x > -0.5*L and z > -T):
			
				y = hull_shape(x,z)	
			
			else:	
				y = 0.
			
			point = Point(x, y, z, idx)
			points.append(point)
			idx+=1
	
		
	return points
	




	
def create_far_points():
	
	
	x_coords = create_x_coords()

	
	z_coords = all_z_coords()
	
	
	#z_coords[2] *= 1.2 
	points = []
	
	
	i = 0
	for x in x_coords:	
		for z in z_coords:
			y = Ly
			#if ( z == -T):
			#	z *= 2.5	
			point = Point(x, y, z, i)
			points.append(point)
			i+=1
	
	
	return points


def create_x_coords():
	
	
	x_coords = np.linspace(-0.5*L, 0.5*L, n_blocks_x + 1) #this should not be linspace
	x_coords = np.insert(x_coords, 0, -0.5*L - Lx_after_hull)
	x_coords = np.append(x_coords, 0.5*L + Lx_before_hull)
	
	
	
	return x_coords


	

def all_z_coords():

	z_coords = np.array([Lz_over_hull, 0., -T, -Lz_below_hull])		
	return z_coords
	

	

def renumber_points(hull_points, far_points):

	#create a 3D array for easier indexing
	points = np.empty((points_imax(), points_jmax(), points_kmax()), dtype=object)
	
	for i in range(points_imax()):
		for k in range(points_kmax()):
		
			j=0
			points[i,j,k] = hull_points[k + i*points_kmax()]
			j=1
			points[i,j,k] = far_points[k + i*points_kmax()]
			
			
	
	#renumber the points
	num = 0
	imax = points_imax(); jmax = points_jmax(); kmax = points_kmax();
	
	for i in range(imax):
		for j in range(jmax):	
			for k in range(kmax):
				points[i,j,k].n = num
				num+=1



	return points


	



def create_blocks(points):

	 
	blocks = []
		
	
	for i in range(0, points_imax() - 1):
		for k in range(points_kmax() - 1):
		
			#all blocks
			j = 0
		
			p1 = points[i,j,k + 1]
			p2 = points[i+1,j,k + 1]
			p3 = points[i+1,j+1,k + 1]
			p4 = points[i,j+1,k + 1]

			p5 = points[i,j,k]
			p6 = points[i+1,j,k]
			p7 = points[i+1,j+1,k]
			p8 = points[i,j + 1,k]
		
			id = ""
			block = Block(p1, p2, p3, p4, p5, p6, p7, p8, 0, id)
			blocks.append(block)

	for block in blocks:
	
		set_block_identifier(block)
	
	
	return blocks
	
def set_block_identifier(block):
	
	x0_hull = 0.5*L
	x1_hull = -0.5*L
	
	waterline_z = 0.
	
	block_x0 = block.p1.x 
	block_z0 = block.p1.z #lowest point of block
	
	
	
	if (block_z0 >= waterline_z):
		id_z = "air"
	elif (block_z0 < -T):
		id_z = "water"
	else:
		id_z = "hull"
		
		
	if (block_x0 >= x0_hull): 		
		id_x = "inlet"

	elif (block_x0 < x1_hull):
		id_x = "outlet"
	
	else:
		id_x = "hull"
	
	
	full_id = "x_{}_z_{}".format(id_x, id_z)
	block.identifier = full_id

def create_edges(points):
	
	edges = []
	
	
	#only for hull blocks
	start_i = 1
	end_i = start_i + n_blocks_x + 1	
	
	#z-direction
	j = 0
	k = 1
	for i in range(start_i, end_i):
			
		p1 = points[i,j,k + 1]
		p2 = points[i,j,k]
		edge = Edge(p1,p2,1)
		edges.append(edge)

	
	
	#x-direction curvature
	start_i = 1
	end_i = start_i + n_blocks_x
	start_k = 0 #also curve the upper domain
	end_k = 3	

	j = 0
	for i in range(start_i, end_i):
		for k in range(start_k, end_k):
			p1 = points[i,j,k]
			p2 = points[i+1,j,k]
			edge = Edge(p1,p2,1)
			edges.append(edge)



	return edges


def create_patches(blocks):

	
	
	faces_far = []
	faces_inlet_water = []
	faces_inlet_air = []
	faces_hull = []
	faces_symmetry = []
	faces_outlet = []
	faces_atm = []
	faces_bottom = []
	faces_freeboard = []
	
	for block in blocks:
		
		#all blocks have the farfield y-normal face
		faces_far.append(block.neg_ynormal_face())
	
		if ((block.identifier == "x_inlet_z_water") or block.identifier == "x_inlet_z_hull"):
			faces_inlet_water.append(block.pos_xnormal_face())	
	
		if (block.identifier == "x_inlet_z_air"):
			faces_inlet_air.append(block.pos_xnormal_face())	
	
		if (block.identifier == "x_hull_z_hull"):
			faces_hull.append(block.pos_ynormal_face())	
	
		if (block.identifier == "x_hull_z_air"):
			faces_freeboard.append(block.pos_ynormal_face())	
	
		if ("x_outlet" in block.identifier):
			faces_outlet.append(block.neg_xnormal_face())
			
		if (("x_outlet" in block.identifier) or ("x_inlet" in block.identifier) or (block.identifier == "x_hull_z_water") ):
			faces_symmetry.append(block.pos_ynormal_face())
	
		 
		if ("z_air" in block.identifier):
			faces_atm.append(block.neg_znormal_face())
	
		if ("z_water" in block.identifier):
			faces_bottom.append(block.pos_znormal_face())
	
			
		
	patches = []
	
		
	
	patches.append(Patch("far", "symmetryPlane", faces_far))
	patches.append(Patch("inlet_water", "patch", faces_inlet_water))
	patches.append(Patch("inlet_air", "patch", faces_inlet_air))
	patches.append(Patch("hull", "wall", faces_hull))
	patches.append(Patch("freeboard", "wall", faces_freeboard))
	patches.append(Patch("outlet", "patch", faces_outlet))
	patches.append(Patch("centerplane", "symmetryPlane", faces_symmetry))
	patches.append(Patch("atmosphere", "symmetryPlane", faces_atm))
	patches.append(Patch("bottom", "symmetryPlane", faces_bottom))
	

	return patches



def set_cell_counts(blocks):
	
	
	for block in blocks:
		
		
		block.ny = Ny
		
		#X-DIRECTION
		if ("inlet" in block.identifier):
			block.nx = Nx_before_hull
	
		elif ("outlet" in block.identifier):
			block.nx = Nx_after_hull

	
		#hull blocks
		else:
		
			block.nx = int(Nx_hull/n_blocks_x)
			
		
		
		#Z-DIRECTION
		if ("z_air" in block.identifier):
			block.nz = Nz_over_hull
		
		
		elif("z_water" in block.identifier):
			block.nz = Nz_below_hull
			
		#hull
		else: 
			block.nz = Nz_hull
		
		

def set_gradings(blocks):
	
	
	
	for block in blocks:
		
		
		gradings = np.zeros(12)
		gradings[:] = 1.0	
		gradings[4:8] = 1.0 #y-grading is set under the Block object
		
		
		if ("x_inlet" in block.identifier):
			gradings[0:4] = Gx_before_hull
	
		if ("x_outlet" in block.identifier):
			gradings[0:4] = 1.0 / Gx_after_hull
	
	
		if ("z_air" in block.identifier):
			gradings[8:] = Gz_over_hull
		
		if("z_water" in block.identifier):
			gradings[8:] = 1.0/Gz_below_hull
		
		
		e8_length = distance_between_points_on_hull(block.p5, block.p1)
		e9_length = distance_between_points_on_hull(block.p6, block.p2)
		
	
		dz_hull = T/Nz_hull
		grading_e8 = get_grading(dz_hull, Nz_hull, e8_length)
		grading_e9 = get_grading(dz_hull, Nz_hull, e9_length)
	
		
		grading_e8 = 1./grading_e8
		grading_e9 = 1./grading_e9
		
		
		
		if (block.identifier == "x_hull_z_hull"):
			gradings[8] = grading_e8
			gradings[9] = grading_e9
		
		if (block.identifier == "x_inlet_z_hull"):
			gradings[8] = grading_e8
			
		if (block.identifier == "x_outlet_z_hull"):	
			gradings[9] = grading_e9
	
	
		block.edge_grading = gradings
	
	
def distance_between_points_on_hull(p1, p2):
	"""
	computes the distance between points p1 and p2 along the hull
	
	"""
	
	x_coords = np.linspace(p1.x, p2.x, 20)
	z_coords = np.linspace(p1.z, p2.z, 20)
	y_coords = B * (1.0 - (z_coords/T)**2) * (1.0 - (2.0*x_coords/L)**2)
	
	
	x1 = x_coords[0:-1]
	x2 = x_coords[1:]
	
	y1 = y_coords[0:-1]
	y2 = y_coords[1:] 
	
	z1 = z_coords[0:-1]
	z2 = z_coords[1:]
	
	
	
	distance =  np.sum( np.sqrt(  (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2 ))
	
	return distance
	
	

	


def create_mesh():

	
	
	
	points_hull = create_hull_points()
	points_far = create_far_points()
	
	points = renumber_points(points_hull, points_far)

	blocks = create_blocks(points)
	edges = create_edges(points)
	patches = create_patches(blocks)
	set_cell_counts(blocks)
	set_gradings(blocks)
	

	#patches = []

	writer = Writer(outputfile)
	writer.write_header()
	writer.write_points(points)
	writer.write_blocks(blocks)
	writer.write_edges(edges)
	writer.write_patches(patches)
	

def hull_shape(x,z):

	
	#below waterline
	if (z <= 0. and z >= -T):
		
		y = B * (1.0 - (z/T)**2) * (1.0 - (2.0*x/L)**2)


	#freeboard
	else:
		y = B * (1.0 - 0) * (1.0 - (2.0*x/L)**2)

	

	return y




	



main()
