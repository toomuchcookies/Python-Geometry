from math import pi, sqrt, acos, cos, sin
import numpy as np

class Point:
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
	
	def abs(self):
		return sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
	
	def dist(self, other):
		temp = other - self
		return sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z)
	
	def cross(self, other):
		return Point(self.y*other.z - self.z*other.y, self.z*other.x - self.x*other.z, self.x*other.y - self.y*other.x)
	
	def dot(self, other):
		return (self.x*other.x + self.y*other.y + self.z*other.z)
		
	def transform(self, T):
		t = [[self.x], [self.y], [self.z], [1]]
		t = T*t
		return Point(t[0,0],t[1,0],t[2,0])
		
	def aslist(self):
		return [self.x, self.y, self.z]
	
	def __sub__(self,other):
		return Point(self.x - other.x, self.y - other.y, self.z - other.z)
	
	def __add__(self,other):
		return Point(self.x + other.x, self.y + other.y, self.z + other.z)
	
	def __neg__(self):
		return Point(-self.x, -self.y, -self.z)
	
	def __mul__(self, other):
		if type(other) == float or type(other) == int:
			return Point(self.x * other, self.y * other, self.z * other)
		else:
			raise TypeError, 'The arguments passed must be Numerical'

	def __rmul__(self, other):
		if type(other) == float or type(other) == int:
			return Point(self.x * other, self.y * other, self.z * other)
		else:
			raise TypeError, 'The arguments passed must be Numerical'
			
	def __div__(self, other):
		if type(other) == float or type(other) == int:
			return Point(self.x/other, self.y/other, self.z/other)
		else:
			raise TypeError, 'The arguments passed must be Numerical'
	
	def __gt__(self, other):
		if isinstance(other, Point):
			return (self.abs() > other.abs())
		elif type(other) == numbers.Real:
			return (self.abs() > other)

	def __lt__(self, other):
		if isinstance(other, Point):
			return (self.abs() < other.abs())
		elif type(other) == numbers.Real:
			return (self.abs() < other)
			
	def __repr__(self):
		return "Point(" + repr(self.x) + ", " + repr(self.y) + ", " + repr(self.z) + ")"

class Line:
	def __init__(self, P, vec=Point(0,0,0)):
		if type(P) == list and len(P) > 1:
			if isinstance(P[0],Point) and isinstance(P[1],Point):
				self.P = P[0]
				P2 = P[1]
				self.vec = P2 - self.P
			else:
				raise TypeError, 'Line expects a list of two Points or a Point a directional Point'
		elif isinstance(P,Point) and vec.abs() > 0:
			self.P = P
			self.vec = vec
		else:
			print P, vec
			raise TypeError, 'Line expects a list of two Points or a Point a directional Point'
		self.P = self.findnearestPoint()
		self.vec = self.vec / self.vec.abs()
	
	def findnearestPoint(self, to=Point(0,0,0)):
		P = self.P
		vec = self.vec
		u = vec.dot(to - P)/vec.dot(vec)
		return (P + u*vec)

	def __repr__(self):
		return "Line(" + repr(self.P) + ", " + repr(self.vec) + ")"
		
	def intersect(self, other):
		if isinstance(other,Plane):
			return other.lineintersect(self)
		else:
			return None
#		elif isinstance(other,Line):
#			return self.lineintersect(other)

class Plane:
	def __init__(self, N, D=0):
		if isinstance(N,Point):
			self.N = N
			self.D = D
		else:
			raise TypeError, 'The arguments passed to Plane must be POINT'

	def __init__(self, P1, P2, P3):
		def chk_type(p_list):
			ret_list = []
			for p in p_list:
				if isinstance(p, Point):
					ret_list.append(p)
				elif type(p) == list and len(p) > 2:
					p = Point(p[0],p[1],p[2])
					ret_list.append(p)
				else:
					ret_list.append(None)
			return ret_list
		
		[P1, P2, P3] = chk_type([P1, P2, P3])
		if None not in [P1, P2, P3]:
			self.N, self.D = self.plane_def(P1, P2, P3)
		else:
			raise TypeError, 'The arguments passed to Plane must be POINT'
	
	def __repr__(self):
		return "Plane(" + repr(self.N) + ", " + repr(self.D) + ")"

	def plane_def(self, p1, p2, p3):
		N = (p2-p1).cross(p3-p1)
		D = (-N).dot(p1)
		return N, D

	def planeintersect(self, other):
		N1 = self.N
		N2 = other.N
		D1 = self.D
		D2 = other.D
		if (N1.cross(N2)).abs() == 0:
			# Planes are parallel
			return None
		else:
			det =  N1.dot(N1) * N2.dot(N2) - (N1.dot(N2))**2
			c1 = D1 * N2.dot(N2) - D2 * N1.dot(N2) / det
			c2 = D2 * N1.dot(N1) - D1 * N1.dot(N2) / det
			return Line((c1 * N1) + (c2 * N2), N1.cross(N2))

	def lineintersect(self,other):
		N = self.N
		D = self.D
		P = other.P
		vec = other.vec
		u1 = N.dot(P) + D
		u2 = N.dot(vec)
		u = u1 / u2
		return P + u * vec

	def intersect(self, other):
		if isinstance(other,Plane):
			return self.planeintersect(other)
		elif isinstance(other,Line):
			return self.lineintersect(other)

if __name__ == '__main__':
	p1 = Point(1,1,0)
	p2 = Point(1,0,0)
	p3 = Point(0,1,0)

	p4 = Point(0,0,1)
	p5 = Point(0,1,1)
	p6 = Point(0,1,0)

	p7 = Point(0,0,1)
	p8 = Point(1,0,1)
	p9 = Point(1,0,0)

	plane1 = Plane(p1,p2,[0,1,0])
	plane2 = Plane(p4,p5,p6)
	plane3 = Plane(p7,p8,p9)

	line1 = plane1.intersect(plane2)
	point1 = plane3.intersect(line1)

	print plane1
	print line1
	print point1
