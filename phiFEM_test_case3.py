from __future__ import print_function
import numpy as np
from dolfin import *
import sympy
import matplotlib.pyplot as plt
parameters['allow_extrapolation'] = True

# Number of iterations
init_Iter = 1
Iter = 7

# parameter of the ghost penalty
sigma = 20
gamma_div = 20
gamma1 = 20
gamma2 = 20

# Polynome Pk
polV = 1
polPhi = 2
parameters["form_compiler"]["quadrature_degree"]=2*(polV+1)+1
parameters["form_compiler"]["representation"] = 'uflacs'

# Ghost penalty
ghost = False

# plot the solution
Plot = True

# Compute the conditioning number
conditioning = False

R = 0.47

def Phi(x,y):
	r = (x**2+y**2)**(0.5)
	if x!=0:
		theta = np.arctan2(y,x)
	if x==0 and y>0:
		theta = 0.5*np.pi
	if x==0 and y<0:
		theta = -0.5*np.pi
	if x==0 and y==0:
		theta = 0 
	R = 0.47
	return r**4*(5.0 + 3.0*np.sin(7.0*theta + 7*np.pi/36.0))/2.0 - R**4

class phi_expression(UserExpression):
	def eval(self, value, x):
		xxx = x[0]
		yyy = x[1]
		rrr = (xxx**2+yyy**2)**(0.5)
		if xxx!=0:
			theta = np.arctan2(yyy,xxx)
		if xxx==0 and yyy>0:
			theta = 0.5*np.pi
		if xxx==0 and yyy<0:
			theta = -0.5*np.pi
		if xxx==0 and yyy==0:
			theta = 0 
		value[0] = rrr**4*(5.0 + 3.0*np.sin(7.0*theta + 7*np.pi/36.0))/2.0 - R**4
	def value_shape(self):
		return (2,)


class g_expression(UserExpression):
	def eval(self, value, x):
		xxx = x[0]
		yyy = x[1]
		rrr = (xxx**2+yyy**2)**(0.5)
		if xxx!=0:
			theta = np.arctan2(yyy,xxx)
		if xxx==0 and yyy>0:
			theta = 0.5*np.pi
		if xxx==0 and yyy<0:
			theta = -0.5*np.pi
		if xxx==0 and yyy==0:
			theta = 0 
		value[0] = (rrr**4*(5.0 + 3.0*np.sin(7.0*theta + 7*np.pi/36.0))/2.0 - R**4)*np.exp(xxx)*np.cos(yyy)
	def value_shape(self):
		return (2,)


# Function used to write in the outputs files
def output_latex(f,A,B):
	for i in range(len(A)):
		f.write('(')
		f.write(str(A[i]))
		f.write(',')
		f.write(str(B[i]))
		f.write(')\n')
	f.write('\n')


# Computation of the Exact solution and exact source term
x, y = sympy.symbols('xx yy')
u_sympy = sympy.exp(y)*sympy.sin(x)
a_sympy = 1.0+x*x+y*y
f_sympy = -sympy.diff(a_sympy*sympy.diff(u_sympy, x),x)-sympy.diff(a_sympy*sympy.diff(u_sympy, y),y)+u_sympy

# Initialistion of the output
size_mesh_vec = np.zeros(Iter)
error_L2_vec = np.zeros(Iter)
error_H1_vec = np.zeros(Iter)
cond_vec = np.zeros(Iter)
for i in range(init_Iter-1,Iter):
	print('##################')
	print('## Iteration ',i+1,'##')
	print('##################')

	# Construction of the mesh
	N = int(10*2**((i)))
	mesh_macro = RectangleMesh(Point(-1.0, -1.0), Point(1.0, 1.0), N, N)
	VP2 = FunctionSpace(mesh_macro, 'CG',polV+2)
	phi = phi_expression(element=VP2.ufl_element())
	phi = interpolate(phi,VP2)
	domains = MeshFunction("size_t", mesh_macro, mesh_macro.topology().dim())
	domains.set_all(0)
	for ind in range(mesh_macro.num_cells()):
		mycell = Cell(mesh_macro,ind)
		v1x,v1y,v2x,v2y,v3x,v3y = mycell.get_vertex_coordinates()
		if Phi(v1x,v1y)<=0 or Phi(v2x,v2y)<=0 or Phi(v3x,v3y)<=0 or phi(v1x,v1y)<=0 or phi(v2x,v2y)<=0 or phi(v3x,v3y)<=0:
			domains[ind] = 1
	mesh = SubMesh(mesh_macro, domains, 1)
	print("num cells mesh :",mesh.num_cells())
	V = FunctionSpace(mesh,'CG',polV)
	VP1 = FunctionSpace(mesh, 'CG',polV)
	VP2 = FunctionSpace(mesh, 'CG',polV+2)
	VP3 = FunctionSpace(mesh, 'CG',polV+7)


	# Computation of the source term
	f_expr = Expression(sympy.ccode(f_sympy).replace('xx', 'x[0]').replace('yy', 'x[1]'),degree=polV+5,domain=mesh)
	u_expr = Expression(sympy.ccode(u_sympy).replace('xx', 'x[0]').replace('yy', 'x[1]'),degree=polV+7,domain=mesh)


	# Construction of phi
	phi = phi_expression(element=VP2.ufl_element())
	phi = interpolate(phi,VP2)
	g = g_expression(element=VP3.ufl_element())
	g = interpolate(g,VP3)
	g = g+u_expr
	a_expr = Expression('1.0+x[0]*x[0]+x[1]*x[1]',degree=polV+5,domain=mesh)


	# Facets and cells where we apply the ghost penalty
	mesh.init(1,2)
	vertex_ghost = MeshFunction("size_t", mesh, mesh.topology().dim()-2)
	facet_ghost = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
	cell_ghost = MeshFunction("size_t", mesh, mesh.topology().dim())
	vertex_ghost.set_all(0)
	facet_ghost.set_all(0)
	cell_ghost.set_all(0)

	for mycell in cells(mesh):
		for myfacet in facets(mycell):
			v1, v2 = vertices(myfacet)
			if Phi(v1.point().x(),v1.point().y())*Phi(v2.point().x(),v2.point().y())<=0 or phi(v1.point().x(),v1.point().y())*phi(v2.point().x(),v2.point().y())<=0:
				cell_ghost[mycell] = 1
				for myfacet2 in facets(mycell):
					facet_ghost[myfacet2] = 1
					v1, v2 = vertices(myfacet2)
					vertex_ghost[v1] = 1
					vertex_ghost[v2] = 1

	for mycell in cells(mesh):
		v1, v2, v3 = vertices(mycell)
		if vertex_ghost[v1]==0 and vertex_ghost[v2]==0 and vertex_ghost[v3]==0:
			cell_ghost[mycell] = 2


	print("Selection of the cells and facets : ok")			
	# Initialize cell function for domains
	dx = Measure("dx")(domain = mesh,subdomain_data = cell_ghost)
	ds = Measure("ds")(domain = mesh)
	dS = Measure("dS")(domain = mesh,subdomain_data = facet_ghost)

	# Initialize cell function for domains
	dx = Measure("dx")(domain = mesh,subdomain_data = cell_ghost)
	ds = Measure("ds")(domain = mesh)
	dS = Measure("dS")(domain = mesh,subdomain_data = facet_ghost)

	# Resolution
	n = FacetNormal(mesh)
	h = CellDiameter(mesh)
	u = TrialFunction(V)
	v = TestFunction(V)
	if ghost == False:
		a = a_expr*inner(grad(phi*u),grad(phi*v))*dx + phi*u*phi*v*dx - a_expr*dot(inner(grad(phi*u),n),phi*v)*ds
		L = f_expr*v*phi*dx+a_expr*inner(grad(-g),grad(phi*v))*dx - g*phi*v*dx - a_expr*dot(inner(grad(-g),n),phi*v)*ds
	if ghost == True:
		a = a_expr*inner(grad(phi*u),grad(phi*v))*dx + phi*u*phi*v*dx - a_expr*dot(inner(grad(phi*u),n),phi*v)*ds + sigma*avg(h)*dot(jump(grad(phi*u),n),jump(grad(phi*v),n))*dS(1)+sigma*h**2*inner(phi*u+div(a_expr*grad(phi*u)),phi*v+div(a_expr*grad(phi*v)))*dx(1)
		L = f_expr*v*phi*dx-sigma*h**2*inner(f_expr,phi*v+div(a_expr*grad(phi*v)))*dx(1)+a_expr*inner(grad(-g),grad(phi*v))*dx  - g*phi*v*dx - a_expr*dot(inner(grad(-g),n),phi*v)*ds + sigma*avg(h)*dot(jump(grad(-g),n),jump(grad(phi*v),n))*dS(1)+sigma*h**2*inner(div(a_expr*grad(-g)),phi*v+div(a_expr*grad(phi*v)))*dx(1)

	# Define solution function
	u_h = Function(V)
	print('ready o solve')
	solve(a == L, u_h)#, solver_parameters={'linear_solver': 'mumps'})
	sol = u_h*phi+g

	# Computation of the error
	approx_L2 = assemble((sol)**2*dx(0)+(sol)**2*dx(2))**0.5
	exact_L2 = assemble((u_expr)**2*dx(0)+(u_expr)**2*dx(2))**0.5
	err_L2 = assemble((sol-u_expr)**2*dx(0)+(sol-u_expr)**2*dx(2))**0.5/assemble((u_expr)**2*dx(0)+(u_expr)**2*dx(2))**0.5
	err_H1 = assemble((grad(sol-u_expr))**2*dx(0)+(grad(sol-u_expr))**2*dx(2))**0.5/assemble((grad(u_expr))**2*dx(0)+(grad(u_expr))**2*dx(2))**0.5
	size_mesh_vec[i] = np.sqrt(2)/N
	error_L2_vec[i] = err_L2
	error_H1_vec[i] = err_H1
	print('h :',np.sqrt(2)/N)
	print('relative L2 error : ',err_L2)
	print('relative H1 error : ',err_H1)
	print('L2 norm of u_h : ',approx_L2)
	print('L2 norm of u : ',exact_L2)	
	if conditioning == True:
		A = np.matrix(assemble(a).array())
		ev, eV = np.linalg.eig(A)
		ev = abs(ev)
		cond = np.max(ev)/np.min(ev)
		cond_vec[i] = cond
		print("conditioning number",cond)
	print('')


# Print the output vectors
print('Vector h :',size_mesh_vec)
print('Vector relative L2 error : ',error_L2_vec)
print('Vector relative H1 error : ',error_H1_vec)
print("conditioning number",cond_vec)


#  Write the output file for latex
if ghost == False:
	f = open('output_no_ghost_case3.txt','w')
if ghost == True:
	f = open('output_ghost2_case3.txt','w')
f.write('relative L2 norm : \n')	
output_latex(f,size_mesh_vec,error_L2_vec)
f.write('relative H1 norm : \n')	
output_latex(f,size_mesh_vec,error_H1_vec)
if conditioning == True:
	f.write('conditioning number  : \n')	
	output_latex(f,size_mesh_vec,cond_vec)
f.close()


# Plot and save
if Plot == True:
	#sol = project(u_h*phi,V)
	plot_sol = plot(sol)
	#file = File('poisson.pvd')
	#file << sol
	plt.savefig('myfig.png')
