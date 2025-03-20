from __future__ import print_function
import numpy as np
from dolfin import *
import matplotlib.pyplot as plt
parameters['allow_extrapolation'] = True
from mshr import *

# Number of rellienement
init_Iter = 1
Iter = 3

# parameter of the ghost penalty
sigma = 100

# Polynome Pk
polV = 1
parameters["form_compiler"]["quadrature_degree"]=2*(polV+1)

# Ghost penalty
ghost = True

# plot the solution
Plot = False

# Compute the conditioning number
conditioning = False



def Omega(x,y):
	return (y<-x/pi+pi and y<pi*x+pi and y>-x/pi-pi and y>pi*x-pi)


# Function used to write in the outputs files
def output_latex(f,A,B):
	for i in range(len(A)):
		f.write('(')
		f.write(str(A[i]))
		f.write(',')
		f.write(str(B[i]))
		f.write(')\n')
	f.write('\n')


#N = 10*2**(Iter+1)
domain_vertices = [Point(0.0, -pi),
                  Point(2.0*pi**2/(pi**2+1),2.0*pi**3/(pi**2+1)-pi),
                  Point(0.0, pi),
                  Point(-2.0*pi**2/(pi**2+1), -2.0*pi**3/(pi**2+1)+pi),
                  Point(0.0, -pi)]
domain = Polygon(domain_vertices)
mesh_exact = generate_mesh(domain,600)
print("num cells mesh exact:",mesh_exact.num_cells())
V_exact = FunctionSpace(mesh_exact,'CG',polV)


# Computation of the source term
f_exact = Expression('1.0',degree=polV,domain=mesh_exact)

class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V_exact, u0, DirichletBoundary())

dx_exact = Measure("dx")(domain = mesh_exact)

# Resolution
u_exact = TrialFunction(V_exact)
v_exact = TestFunction(V_exact)
a_exact = inner(grad(u_exact),grad(v_exact))*dx_exact
L_exact = f_exact*v_exact*dx_exact
u_h_exact = Function(V_exact)
solve(a_exact == L_exact, u_h_exact,bcs=bc,solver_parameters={'linear_solver': 'mumps'})



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
	mesh_macro = RectangleMesh(Point(-4.0, -4.0), Point(4.0, 4.0), 4*N, 8*N)
	domains = MeshFunction("size_t", mesh_macro, mesh_macro.topology().dim())
	domains.set_all(0)
	for ind in range(mesh_macro.num_cells()):
		mycell = Cell(mesh_macro,ind)
		v1x,v1y,v2x,v2y,v3x,v3y = mycell.get_vertex_coordinates()
		if Omega(v1x,v1y) or Omega(v2x,v2y) or Omega(v3x,v3y):
			domains[ind] = 1
	mesh = SubMesh(mesh_macro, domains, 1)
	print("num cells mesh :",mesh.num_cells())
	V = FunctionSpace(mesh,'CG',polV)

	# Construction of phi
	phi = Expression('-(x[1]+x[0]/pi-pi)*(x[1]+x[0]/pi+pi)*(x[1]-pi*x[0]-pi)*(x[1]-pi*x[0]+pi)',degree=4*polV,domain=mesh)
	phi = interpolate(phi,V)

	# Computation of the source term
	f_expr = Expression('1.0',degree=polV,domain=mesh)

	# Facets and cells where we apply the ghost penalty
	mesh.init(1,2)
	facet_ghost = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
	cell_ghost = MeshFunction("size_t", mesh, mesh.topology().dim())
	facet_ghost.set_all(0)
	cell_ghost.set_all(0)
	for mycell in cells(mesh):
		for myfacet in facets(mycell):
			v1, v2 = vertices(myfacet)
			if phi(v1.point().x(),v1.point().y())*phi(v2.point().x(),v2.point().y())<0:
				cell_ghost[mycell] = 1
				for myfacet2 in facets(mycell):
					facet_ghost[myfacet2] = 1

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
		a = inner(grad(phi*u),grad(phi*v))*dx - dot(inner(grad(phi*u),n),phi*v)*ds
		L = f_expr*v*phi*dx
	if ghost == True:
		a = inner(grad(phi*u),grad(phi*v))*dx - dot(inner(grad(phi*u),n),phi*v)*ds + sigma*avg(h)*dot(jump(grad(phi*u),n),jump(grad(phi*v),n))*dS(1)+sigma*h**2*inner(div(grad(phi*u)),div(grad(phi*v)))*dx(1)
		L = f_expr*v*phi*dx-sigma*h**2*inner(f_expr,div(grad(phi*v)))*dx(1)
	u_h = Function(V)
	#problem = LinearVariationalProblem(a, L, u_h)
	#solver = LinearVariationalSolver(problem)
	#solver = KrylovSolver("mumps")
	#solver.parameters["relative_tolerance"] = 5e-6
	#solver.parameters["maximum_iterations"] = 1000
	#solver.parameters["monitor_convergence"] = True
	#solver.solve()
	solve(a == L, u_h,solver_parameters={'linear_solver': 'mumps'})
	sol = u_h*phi

	# Computation of the error
	Iu_h = project(sol,V_exact)
	err_L2 = assemble((Iu_h-u_h_exact)**2*dx_exact)**0.5/assemble((u_h_exact)**2*dx_exact)**0.5
	err_H1 = assemble((grad(Iu_h-u_h_exact))**2*dx_exact)**0.5/assemble((grad(u_h_exact))**2*dx_exact)**0.5
	size_mesh_vec[i] = np.sqrt(2)/N
	error_L2_vec[i] = err_L2
	error_H1_vec[i] = err_H1
	print('h :',np.sqrt(2)/N)
	print('relative L2 error : ',err_L2)
	print('relative H1 error : ',err_H1)	
	#print('L2 norm : ',assemble((u_expr)**2*dx(0))**0.5)	
	if conditioning == True:
		A = np.matrix(assemble(a).array())
		ev, eV = np.linalg.eig(A)
		ev = abs(ev)
		#cond = mesh.hmax()**2*np.max(ev)/np.min(ev)
		cond = np.max(ev)/np.min(ev)
		cond_vec[i] = cond
		print("conditioning number x h^2",cond)
	print("num cells mesh",mesh.num_cells())
	print('')


# Print the output vectors
print('Vector h :',size_mesh_vec)
print('Vector relative L2 error : ',error_L2_vec)
print('Vector relative H1 error : ',error_H1_vec)
print("conditioning number",cond_vec)



#  Write the output file for latex
if ghost == False:
	f = open('output_no_ghost_case2.txt','w')
if ghost == True:
	f = open('output_ghost2_case2.txt','w')
f.write('relative L2 norm : \n')	
output_latex(f,size_mesh_vec,error_L2_vec)
f.write('relative H1 norm : \n')	
output_latex(f,size_mesh_vec,error_H1_vec)
f.write('conditioning number  : \n')	
output_latex(f,size_mesh_vec,cond_vec)
f.close()


# Plot and save
if Plot == True:
	sol = project(u_h*phi,V)
	plot_sol = plot(sol)
	file = File('poisson.pvd')
	file << sol
	plt.savefig('myfig.png')
