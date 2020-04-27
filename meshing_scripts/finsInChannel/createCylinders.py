

import params
import numpy as np
import pylab as pl
import sys

"""
Input
"""
#channel height
H = params.H

L_before = params.L_before
L_after = params.L_after


D = params.D



S_y = params.S_y
S_x = params.S_x

N_x = params.N_x + 2
N_y = params.N_y + 2


bl_height = params.bl_height


N_theta = params.Ncells_theta
N_bl = params.Ncells_bl
g_bl = params.grading_bl

N_sy = params.Ncells_sy
N_sx = params.Ncells_sx

Ncells_x_before = params.Ncells_x_before
Ncells_x_after = params.Ncells_x_after

Ncells_ch_bl = params.Ncells_ch_bl
Ncells_ch_bl_next = params.Ncells_ch_bl_next

grading_x_before = params.grading_x_before
grading_x_after = params.grading_x_after


Ncells_x_array_after = params.Ncells_x_array_after

base_size = params.base_size


"""
Computed
"""
R = 0.5*D
sin45 = 1./np.sqrt(2.0)
cos45 = 1./np.sqrt(2.0)
R_bl = params.R_bl







class Point(object):
	
	
	def __init__(self, x, y, z, num):
	
		self.x = x
		self.y = y
		self.z = z
		self.num = num
		
	
	def Print(self):
	
		print "Point({}) = {{ {}, {}, {} , {}}};".format(self.num, self.x, self.y, self.z, base_size)
		
	

class Line(object):
	
	def __init__(self, p1, p2, num, id, N, grading, transfinite=True):


		self.p1 = p1; self.p2 = p2;
		self.num = num;
		self.id = id;
		self.N = N
		self.grading = grading

		self.transfinite = transfinite

		
		

	def Print(self):
		
		string = "Line({}) = {{ {}, {} }};".format(self.num, self.p1.num, self.p2.num)

		if self.transfinite:
			string += "Transfinite Line {{{}}} = {} Using Progression {};".format(self.num, self.N, self.grading)

		
		print string

	

class Arc(object):
	
	def __init__(self, p1, p2, p3, num, id, N, grading, transfinite=True):


		self.p1 = p1; self.p2 = p2; self.p3 = p3
		self.num = num;
		self.id = id;
		
		self.N = N
		self.grading = grading

		self.transfinite = transfinite


	def Print(self):
		
		string = "Circle({}) = {{ {}, {}, {} }};".format(self.num, self.p1.num, self.p2.num, self.p3.num)
		
		if self.transfinite:
			string += "Transfinite Line {{{}}} = {} Using Progression {};".format(self.num, self.N, self.grading)		


		print string

class Surface(object):


	def __init__(self, l1,l2,l3,l4, num, id, transfinite=True):

		self.l1 = l1; self.l2 = l2; self.l3 = l3; self.l4 = l4;
		self.num = num;
		
		self.id = id
		self.transfinite = transfinite

	def Print(self):
		string = ""
		if "connect" in self.id:
			string += "Line Loop({}) = {{{},{},{},{}}};".format(self.num, self.l1.num, -self.l2.num, -self.l3.num, -self.l4.num)
			
		elif "SQUARE" in self.id:
			string += "Line Loop({}) = {{{},{},{},{}}};".format(self.num, self.l1.num, self.l2.num, self.l3.num, self.l4.num)	
			
		else:	
			string += "Line Loop({}) = {{{},{},{},{}}};".format(self.num, self.l1.num, self.l2.num, -self.l3.num, -self.l4.num)
			
			
		string += " Plane Surface({}) = {{ {} }};".format(self.num, self.num)

		if self.transfinite:
			string += " Transfinite Surface {{ {} }};".format(self.num)
			string += " Recombine Surface {{ {} }};".format(self.num)

		print string


class Square(object):


	def __init__(self, point, num, id, height, width):

		self.center = point
		self.num = num;
		
		self.id = id
		self.height = height
		self.width = width
		
		

		"""
		POINTS
		"""
				
		#points on bl
		self.NW_bl, self.NE_bl, self.SE_bl, self.SW_bl = self.createBlPoints()
		
		"""
		LINES
		"""
		
		#BL lines
		self.line_N, self.line_E, self.line_S, self.line_W = self.createBlLines()
		
		

		"""
		SURFACES
		"""
		self.surf = self.createAllSurfaces()

	
	def createAllSurfaces(self):
		surf = Surface(self.line_N, self.line_E, self.line_S, self.line_W, self.num+4, "Surface_SQUARE", transfinite=True)
		
		return surf
		
	
		
	def createBlLines(self):
		
		nCells_x = N_theta
		nCells_y = N_theta
		grading_x = 1.0
		grading_y = 1.0
		#idd += "_chBL_ext_before"
		if ("chBL" in self.id):
			
			nCells_y = Ncells_ch_bl
			grading_y = 1.0
		
		if ("ext_before" in self.id):
			
			nCells_x = Ncells_x_before
			grading_x = grading_x_before
			
		if ("ext_after" in self.id):
			
			nCells_x = Ncells_x_after	
			grading_x = grading_x_after
		
		
			
			
		line_N = Line(self.NW_bl, self.NE_bl, self.num + 5, "Line_NORTH", nCells_x, 1./grading_x, transfinite=True)
		line_E = Line(self.NE_bl, self.SE_bl, self.num + 6, "Line_EAST", nCells_y, grading_y, transfinite=True)
		line_S = Line(self.SE_bl, self.SW_bl, self.num + 7, "Line_SOUTH", nCells_x, grading_x, transfinite=True)
		line_W = Line(self.SW_bl, self.NW_bl, self.num + 8, "Line_WEST", nCells_y, grading_y, transfinite=True)	

		return line_N, line_E, line_S, line_W	
		

		
	def createBlPoints(self):
		
		NW = Point(-0.5*self.width + self.center.x, 0.5*self.height + self.center.y, 0., self.num + 5)
		NE = Point(0.5*self.width + self.center.x, 0.5*self.height + self.center.y, 0., self.num + 6)
		SE = Point(0.5*self.width + self.center.x, -0.5*self.height + self.center.y, 0., self.num + 7)
		SW = Point(-0.5*self.width + self.center.x, -0.5*self.height + self.center.y, 0., self.num + 8)	
		
		return NW, NE, SE, SW
		
		
		
	def printPoints(self):
		self.center.Print()
		
		self.NW_bl.Print(); self.NE_bl.Print(); self.SE_bl.Print(); self.SW_bl.Print();	
		

	def printLines(self):
		
		
		self.line_N.Print(); self.line_E.Print(); self.line_S.Print(); self.line_W.Print();
		


	def printSurfaces(self):
		self.surf.Print();
		

	def Print(self):
		self.printPoints()	
		self.printLines()
		self.printSurfaces()
		
		
	def plotCenterPoint(self):
		
		pl.scatter(self.center.x, self.center.y)	


	



class Cylinder(object):


	def __init__(self, point, num, id, R, R_bl):

		self.center = point
		self.num = num;
		
		self.id = id
		self.R = R
		self.R_bl = R_bl
		

		"""
		POINTS
		"""
		#points on arc circles
		self.NW, self.NE, self.SE, self.SW = self.createCircleArcPoints()		
		#points on bl
		self.NW_bl, self.NE_bl, self.SE_bl, self.SW_bl = self.createBlPoints()
		
		"""
		LINES
		"""
		#circle arc lines
		self.arc_N, self.arc_E, self.arc_S, self.arc_W = self.createCircleArcs()
		
		#BL lines
		self.line_N, self.line_E, self.line_S, self.line_W = self.createBlLines()
		#BL line connectors
		self.line_NW, self.line_NE, self.line_SE, self.line_SW = self.createBlLineConnectors()
		

		"""
		SURFACES
		"""
		self.surf_N, self.surf_E, self.surf_S, self.surf_W = self.createAllSurfaces()

	
	def createAllSurfaces(self):
		surf_N = Surface(self.line_NW, self.line_N, self.line_NE, self.arc_N, self.num+1, "Surface_NORTH", transfinite=True)
		surf_E = Surface(self.line_NE, self.line_E, self.line_SE, self.arc_E, self.num+2, "Surface_EAST", transfinite=True)
		surf_S = Surface(self.line_SE, self.line_S, self.line_SW, self.arc_S, self.num+3, "Surface_SOUTH", transfinite=True)	
		surf_W = Surface(self.line_SW, self.line_W, self.line_NW, self.arc_W, self.num+4, "Surface_WEST", transfinite=True)
	
		return surf_N, surf_E, surf_S, surf_W
		
	def createCircleArcs(self):
		
		
		arc_N = Arc(self.NW, self.center, self.NE, self.num + 1, "Circle_NORTH", N_theta, 1.0, transfinite=True)
		arc_E = Arc(self.NE, self.center, self.SE, self.num + 2, "Circle_EAST", N_theta, 1.0, transfinite=True)
		arc_S = Arc(self.SE, self.center, self.SW, self.num + 3, "Circle_SOUTH", N_theta, 1.0, transfinite=True)
		arc_W = Arc(self.SW, self.center, self.NW, self.num + 4, "Circle_WEST", N_theta, 1.0, transfinite=True)	
		
		return arc_N, arc_E, arc_S, arc_W
		
	def createBlLines(self):
		
		
		if (params.BL_arc):
			line_N = Arc(self.NW_bl, self.center, self.NE_bl, self.num + 5, "Line_NORTH", N_theta, 1.0, transfinite=True)
			line_E = Arc(self.NE_bl, self.center, self.SE_bl, self.num + 6, "Line_EAST", N_theta, 1.0, transfinite=True)
			line_S = Arc(self.SE_bl, self.center, self.SW_bl, self.num + 7, "Line_SOUTH", N_theta, 1.0, transfinite=True)
			line_W = Arc(self.SW_bl, self.center, self.NW_bl, self.num + 8, "Line_WEST", N_theta, 1.0, transfinite=True)	
		
		else:
			line_N = Line(self.NW_bl, self.NE_bl, self.num + 5, "Line_NORTH", N_theta, 1.0, transfinite=True)
			line_E = Line(self.NE_bl, self.SE_bl, self.num + 6, "Line_EAST", N_theta, 1.0, transfinite=True)
			line_S = Line(self.SE_bl, self.SW_bl, self.num + 7, "Line_SOUTH", N_theta, 1.0, transfinite=True)
			line_W = Line(self.SW_bl, self.NW_bl, self.num + 8, "Line_WEST", N_theta, 1.0, transfinite=True)
		
		return line_N, line_E, line_S, line_W	
		

	def createBlLineConnectors(self):
		
		


		line_NW = Line(self.NW, self.NW_bl, self.num + 9, "Line_NORTHWEST", N_bl, g_bl, transfinite=True)
		line_NE = Line(self.NE, self.NE_bl, self.num + 10, "Line_NORTHEAST", N_bl, g_bl, transfinite=True)
		line_SE = Line(self.SE, self.SE_bl, self.num + 11, "Line_SOUTHEAST", N_bl, g_bl, transfinite=True)
		line_SW = Line(self.SW, self.SW_bl, self.num + 12, "Line_SOUTHWEST", N_bl, g_bl, transfinite=True)		
		
		return line_NW, line_NE, line_SE, line_SW

	def createCircleArcPoints(self):
		
				
	

		NW = Point(-cos45*self.R + self.center.x, sin45*self.R + self.center.y, 0., self.num + 1)
		NE = Point(cos45*self.R + self.center.x, sin45*self.R + self.center.y, 0., self.num + 2)
		SE = Point(cos45*self.R + self.center.x, -sin45*self.R + self.center.y, 0., self.num + 3)
		SW = Point(-cos45*self.R + self.center.x, -sin45*self.R + self.center.y, 0., self.num + 4)
		
		return NW, NE, SE, SW	
		
	def createBlPoints(self):
		
		NW = Point(-cos45*self.R_bl + self.center.x, sin45*self.R_bl + self.center.y, 0., self.num + 5)
		NE = Point(cos45*self.R_bl + self.center.x, sin45*self.R_bl + self.center.y, 0., self.num + 6)
		SE = Point(cos45*self.R_bl + self.center.x, -sin45*self.R_bl + self.center.y, 0., self.num + 7)
		SW = Point(-cos45*self.R_bl + self.center.x, -sin45*self.R_bl + self.center.y, 0., self.num + 8)	
		
		return NW, NE, SE, SW
		
		
		
	def printPoints(self):
		self.center.Print()
		self.NW.Print(); self.NE.Print(); self.SE.Print(); self.SW.Print();
		self.NW_bl.Print(); self.NE_bl.Print(); self.SE_bl.Print(); self.SW_bl.Print();	
		

	def printLines(self):
		
		self.arc_N.Print(); self.arc_E.Print(); self.arc_S.Print(); self.arc_W.Print();
		self.line_N.Print(); self.line_E.Print(); self.line_S.Print(); self.line_W.Print();
		self.line_NW.Print(); self.line_NE.Print(); self.line_SE.Print(); self.line_SW.Print();	


	def printSurfaces(self):
		self.surf_N.Print(); self.surf_E.Print(); self.surf_S.Print(); self.surf_W.Print();
		

	def Print(self):
		self.printPoints()	
		self.printLines()
		self.printSurfaces()
		
		
	def plotCenterPoint(self):
		
		pl.scatter(self.center.x, self.center.y)	


	def plotArcPoints(self):
		pl.scatter(self.NW.x, self.NW.y)
		pl.scatter(self.NE.x, self.NE.y)
		pl.scatter(self.SE.x, self.SE.y)
		pl.scatter(self.SW.x, self.SW.y)
		
		
	def plotBlPoints(self):
		pl.scatter(self.NW_bl.x, self.NW_bl.y)
		pl.scatter(self.NE_bl.x, self.NE_bl.y)
		pl.scatter(self.SE_bl.x, self.SE_bl.y)
		pl.scatter(self.SW_bl.x, self.SW_bl.y)
		

def main():

	
	
	cylinders = createCylinders()
	
	for cylinder in cylinders:
		cylinder.Print()
	
	
	#these are arbitrary
	nLines = 10000
	nSurfaces = 10000
	
	
	nLines, nSurfaces, lines_left, lines_right = connectUp(cylinders, nLines, nSurfaces)
	nLines, nSurfaces, lines_down, lines_up = connectRight(cylinders,nLines,nSurfaces)
	
	nSurfaces = connectMiddle(lines_left, lines_right, lines_down, lines_up, nSurfaces)


	createBackground(cylinders, nSurfaces, lines_up, lines_down, lines_left, lines_right)
	
def createBackground(cylinders, nSurfaces, lines_up, lines_down, lines_left, lines_right):

		
	#inlet
	string_inlet =  "Physical Line(1) = { "
	
	for cylinder in cylinders:
		
		if "ext_before" in cylinder.id:
			string_inlet += "{}, ".format(cylinder.line_W.num)
	
	
	for j in range(N_y-1):
		string_inlet += "{}, ".format(lines_left[0 + N_x*j].num)
			
	string_inlet = string_inlet[0:-2] + "};"
	print string_inlet
	
	
	
	#channel walls
	string_wall = "Physical Line(2) = { "
	for cylinder in cylinders:
	
		if ("chBL" in cylinder.id):
			if ("y0" in cylinder.id):
				string_wall += "{}, ".format(cylinder.line_S.num)
			if ("y1" in cylinder.id):
				string_wall += "{}, ".format(cylinder.line_N.num)
	
	
	for i in range(N_x-1):
		j = 0
		string_wall += "{}, ".format(lines_down[i + N_x*j].num)
	
	for i in range(N_x-1):
		j = N_y-1
		string_wall += "{}, ".format(lines_up[i + (N_x - 1)*j].num)
	
	
	string_wall = string_wall[0:-2] + "};"
	print string_wall
	
	
	#outlet
	string_outlet =  "Physical Line(3) = { "
	for cylinder in cylinders:
		
		if "ext_after" in cylinder.id:
			string_outlet += "{}, ".format(cylinder.line_E.num)
	
	
	for j in range(N_y-1):
		string_outlet += "{}, ".format(lines_right[N_x-1 + N_x*j].num)
	
	string_outlet = string_outlet[0:-2] + "};"
	print string_outlet
	
	
	#cylinder walls
	string_cyls = "Physical Line (4) = { "
	for cylinder in cylinders:
		if (not "null" in cylinder.id):
			string_cyls += "{}, {}, {}, {}, ".format(cylinder.arc_N.num, cylinder.arc_E.num, cylinder.arc_S.num, cylinder.arc_W.num)
	
	string_cyls = string_cyls[0:-2] + "};"
	print string_cyls
	
	
	#all surfaces physical
	string = "Physical Surface(10) = { "
	for cylinder in cylinders:
		if ("null" in cylinder.id):
			string += "{}, ".format(cylinder.surf.num)
		else:	
			string += "{}, {}, {}, {}, ".format(cylinder.surf_N.num, cylinder.surf_E.num, cylinder.surf_S.num, cylinder.surf_W.num)
		
	for i in range(10000,nSurfaces):
		string += "{}, ".format(i) 	
	
	string = string[0:-2] + "};"	
	print string
	

def connectMiddle(lines_left, lines_right, lines_down, lines_up, nSurfaces):
	

	for j in range(N_y-1):
		for i in range(N_x-1):
		
			
			idx = i + N_x*j
			l1 = lines_right[idx]
			
			idx = i + (N_x - 1)*(j+1)
			l2 = lines_down[idx]
		
			idx = (i+1) + N_x*j
			l3 = lines_left[idx]
			
			idx = i + (N_x - 1)*j
			l4 = lines_up[idx]
		
			surface = Surface(l1, l2, l3, l4, nSurfaces, "middle", transfinite=True)
			
			surface.Print()		
			
			nSurfaces += 1
			
			
	return nSurfaces		

def connectRight(cylinders, nLines, nSurfaces):
	
	

	lines_down = []
	lines_up = []
	
	for j in range(N_y):
		for i in range(N_x-1):
	
			idx = i + N_x*j
			idx_right = (i + 1) + N_x*j
	
			p1_down = cylinders[idx].SE_bl
			p2_down = cylinders[idx_right].SW_bl
			
			
			nCells = N_sx
			
			if (i == 0 or i == N_x - 2):
				nCells = Ncells_x_array_after
			
			
			line_down = Line(p1_down, p2_down, nLines, "line_connect_down_right", nCells, 1.0, transfinite=True)
			line_down.Print() 
			
			lines_down.append(line_down)
			
			
			p1_up = cylinders[idx].NE_bl
			p2_up = cylinders[idx_right].NW_bl
			
			line_up = Line(p1_up, p2_up, nLines+1, "line_connect_down_right", nCells, 1.0, transfinite=True)
			line_up.Print() 
			
			lines_up.append(line_up)
			
			l1_right = line_up
			l2_right = cylinders[idx_right].line_W
			l3_right = line_down
			l4_right = cylinders[idx].line_E
			
			surface_right = Surface(l1_right, l2_right, l3_right, l4_right, nSurfaces, "surface_connect_right", transfinite=True)
			
			surface_right.Print()
			
			
			nLines+=2
			nSurfaces+=1

	return nLines, nSurfaces, lines_down, lines_up

def connectUp(cylinders, nLines, nSurfaces):
	
	lines_left = []
	lines_right = []
	for j in range(N_y - 1):
		for i in range(N_x):
	
			idx = i + N_x*j
			idx_up = i + N_x*(j+1)
			
			nCells = N_sy
			
			if (j == 0 or j == N_y - 2):
				nCells = Ncells_ch_bl_next
			
			
			
			p1_left = cylinders[idx].NW_bl
			p2_left = cylinders[idx_up].SW_bl
			
			line_left = Line(p1_left, p2_left, nLines, "line_connect_up_left", nCells, 1.0, transfinite=True)
			line_left.Print()
			lines_left.append(line_left)
			
			
			p1_right = cylinders[idx].NE_bl
			p2_right = cylinders[idx_up].SE_bl
			
			line_right = Line(p1_right, p2_right, nLines + 1, "line_connect_up_right", nCells, 1.0, transfinite=True)
			line_right.Print()
			
			lines_right.append(line_right)
			
			
			
			l1_up = line_left
			l2_up = cylinders[idx_up].line_S
			l3_up = line_right
			l4_up = cylinders[idx].line_N
			
			surface_up = Surface(l1_up, l2_up, l3_up, l4_up, nSurfaces, "surface_connect_up", transfinite=True)
			
			surface_up.Print()
		
			nLines+=2
			nSurfaces+=1
		
	return nLines, nSurfaces, lines_left, lines_right



def createCylinders():

	#Returns a list of cylinder objects
	
	L_x = (N_x - 1) * S_x
	L_y = (N_y - 1) * S_y
	

	y_coords = np.linspace(-0.5*L_y, 0.5*L_y, N_y)
	x_coords = np.linspace(0., L_x, N_x)
	
	

	
	y_coords[0] = -0.5*H + 0.5*bl_height
	y_coords[-1] = 0.5*H - 0.5*bl_height
	
	
	
	
	x_coords[0] = -0.5*L_before
	x_coords[-1] = 0 + L_x + 0.5*L_after
	
	x_coords -= S_x
	
	
	ext_width1 = L_before
	ext_width2 = L_after
	
	cylinders = []
	
	
	
	
	counter = 1
	for j in range(N_y):
		
		for i in range(N_x):
			
			y = y_coords[j]
			x = x_coords[i]
			z = 0.
			new_p = Point(x,y,z,counter)
			
			if (i == 0 or i == N_x - 1 or j == 0 or j == N_y - 1):
			
			
				idd = "pin_null"
			
				if (i == 0 and j == 0):
					
					height = bl_height
					width = ext_width1
					idd += "_chBL_ext_before_y0"
					
				elif (i == 0 and j == N_y - 1):
				
					height = bl_height
					width = ext_width1
					idd += "_chBL_ext_before_y1"
					
				elif (i == N_x - 1 and j == N_y - 1):
				
					height = bl_height
					width = ext_width2
					idd += "_chBL_ext_after_y1"
			
			
				elif (i == N_x - 1 and j == 0):
			
					height = bl_height
					width = ext_width2
					idd += "_chBL_ext_after_y0"
					
			
				elif (i == 0):
					height = 2*cos45*R_bl
					width = ext_width1
					idd += "_ext_before"
			
				elif (i == N_x - 1):
					height = 2*cos45*R_bl
					width = ext_width2
					idd += "_ext_after"
				
			
				elif (j == 0):
					height = bl_height
					width = 2*cos45*R_bl
					idd += "_chBL_y0"
			
				elif (j == N_y - 1):
					height = bl_height
					width = 2*cos45*R_bl
					idd += "_chBL_y1"
			
				#print idd
			
				new_cyl = Square(new_p, counter, idd, height, width)
			
				
				
				
			else:			
				new_cyl = Cylinder(new_p, counter, "pin", R, R_bl)
				
				
			cylinders.append(new_cyl)
			counter +=15
	
	#sys.exit()
	return cylinders
	
	
	
	
main()			
	
	
	
		
		
		
