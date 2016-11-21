# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:34:16 2016

@author: afalaize
"""


import dolfin as dol
import numpy as np
import matplotlib.pyplot as plt


pars = {'k': 5.,                   # thermal diffusion
        'T_nord': 1.,              # boundary temperature 1
        'T_sud': 0.,               # boundary temperature 2
        'T_est': 0.,               # boundary temperature 1
        'T_ouest': 0.,             # boundary temperature 2
        'Amp_init': 10.,           # Amplitude of inital condition (gaussian)
        'ts': 5e-3,                # sample rate
        'T': 2e-1,                 # simu duration
        't0': 0.,                  # init time
        'res_folder': 'results',   # folder to store results
        'TS_name': 'TS',           # TimeSerie label
        'pvd_name': 'u',           # PVD label
        'Nx': 10,                  # number of elements in first direction
        'Ny': 10,                  # number of elements in second direction
        'N_grid': 1,               # number of grid oscillation periods in the simu duration
        'A_grid': 0.1,             # grid oscillation amplitude
        }

def u_init():
	return dolfin.Expression('A*exp(pow(x[0]-0.5, 2) + pow(x[1]-0.5, 2))', A=dolfin.Constant(pars['Amp_init']))

# init boundary conditions on function spaces
def BC_u(FS_):
    # Base Boundary Conditions
    bc_sud = dolfin.DirichletBC(FS_, pars['T_sud'],
                                "on_boundary && (x[1] < DOLFIN_EPS )")
    bc_nord = dolfin.DirichletBC(FS_,
    							 pars['T_nord'],
                                 "on_boundary && (x[1] > 1-DOLFIN_EPS )")
    bc_ouest = dolfin.DirichletBC(FS_,
    							  pars['T_ouest'],
                                  "on_boundary && (x[0] < DOLFIN_EPS )")
    bc_est = dolfin.DirichletBC(FS_,
    							pars['T_est'],
                                "on_boundary && (x[0] > 1-DOLFIN_EPS )")
    return [bc_nord, bc_sud, bc_est, bc_ouest]

vars_ = {'u': {'type': 'scalar',
               'elements': ("Lagrange", 1),
               'BC': BC_u,
               'init': u_init()
               },
        }


def mesh():
    # Base mesh
    return dolfin.UnitSquareMesh(pars['Nx'],
                                 pars['Ny'])

class Functions(object):
    def __init__(self, fem, mesh):
        self.mesh = mesh
        # Define functions with their functionspaces
        for name in fem.vars_.keys():
            type_ = fem.vars_[name]['type']
            setattr(self, name+'_type', type_)
            elements = fem.vars_[name]['elements']
            setattr(self, name+'_elements', elements)
            init = fem.vars_[name]['init']
            setattr(self, name+'_init', init)
            tempBC = fem.vars_[name]['BC']
            if type_ == 'vector':
                prefix = 'Vector'
            elif type_ == 'scalar':
                prefix = ''
            # Dolfin FunctionSpace
            funcspace = eval('dolfin.' + prefix + 'FunctionSpace')
            def builder(funcspace_, elements_, name_, tempBC_):
                def FS():
                    return funcspace_(self.mesh, *elements_)
                def trial():
                    return dolfin.TrialFunction(getattr(self, name_+'_FS')())
                def test():
                    return dolfin.TestFunction(getattr(self, name_+'_FS')())
                def func():
                    return dolfin.Function(getattr(self, name_+'_FS')())
                def BC():
                    return tempBC_(getattr(self, name_+'_FS')())
                return FS, trial, test, func, BC
            FS, trial, test, func, BC = builder(funcspace, elements, name, tempBC)
            setattr(self, name+'_FS', FS)
            setattr(self, name+'_trial', trial)
            setattr(self, name+'_test', test)
            setattr(self, name, func)
            setattr(self, name+'0', func)
            setattr(self, name+'_BC', BC)

class Fem(object):
    """
Class that contains all Dolfin routines for the FE simulation of PDE \
problems with Finite elements method

Usage
------
fe = Fem(mesh, vars_, pars)

Parameters
----------

mesh :
    dolfin mesh

vars_ :
    dictionary with keys as variables names and values are dictionary with keys

    * 'type': one of 'scalar' or 'vector';
    * 'elements': dolfin element tupple, e.g. ("Lagrange", 1);
    * 'BC': a list of Dolfin.DirichletBC
    * 'init': a Dolfin.Function which is used as initial condition

pars :
    dictionary of physical and simulation parameters

   """
    def __init__(self, mesh, vars_, pars):
        # physical parameters
        self.pars = pars
        self.vars_ = vars_

        # base mesh geometry dimension
        self.mesh = mesh
        # init spatial mesh 'x',  material mesh 'X' and reference mesh 'xi'
        # with base mesh self.mesh()
        self.mesh_names = ['space', 'material', 'reference']
        for name in self.mesh_names:
        	setattr(self, name, Functions(self, mesh))

        self.boundaryConditions = dict()
#        for name in self.mesh_names:
#            self.boundaryConditions.update({name: BCs(self.FunctionSpaces[name])})

#        # Define source term
#        src_const = '0,0'  # float in C syntax
#        self.source = dolfin.Expression(src_const)
#
#        # reference mesh velocity update
#        def mesh_xi_velocity(t):
#            freq = self.sim_pars['N_periods_grid']/self.sim_pars['duration']
#            Amp_velocity = self.phy_pars['A_grid']*numpy.cos(2*numpi.pi*freq*t)
#            velocity  = dolfin.Expression(("A*sin(x[0]*(2*3.14))",
#                                           "A*sin(x[1]*(2*3.14))"),
#                                          A=dolfin.Constant(Amp_velocity))
#            self.mesh_xi_velocity.assign(dolfin.project(velocity,
#                                                        FunctionSpace(self.meshes['xi'])))
#
#        setattr(self, 'xi_velocity', xi_velocity)
#
#        # Define the weak formulation for the heat equation in the ALE reference mesh
#        transient_term = (self.sim_pars['ts']**-1)*dolfin.inner(self.trial_funcs['xi']-self.trial_funcs_initial_guess['xi'], self.test_funcs['xi'])*dolfin.dx
#        diffusion_term = self.phy_pars['k']*dolfin.inner(dolfin.grad(self.trial_func),
#                                                         dolfin.grad(self.test_func))*dolfin.dx
#        ALE_term = self.phy_pars['k']*dolfin.inner(dolfin.grad(self.trial_func),
#                                                         dolfin.grad(self.test_func))*dolfin.dx
#        source_term = -self.source*self.test_func*dolfin.dx
#        self.weak_form = transient_term + diffusion_term + source_term
#
#        A = dolfin.assemble(dolfin.dot(self.trial_func, self.test_func)*dolfin.dx)
#        B = dolfin.assemble(diffusion_term)
#        f = dolfin.assemble(dolfin.dot(self.source, self.test_func)*dolfin.dx)
#
#        self.assembly_matrices = {'A': A.array(),
#                                  'B': B.array(),
#                                  'f': f.array()}
#
#    def run(self, hold_plot=False, rescale=True):
#
#        lhs = dolfin.lhs(self.weak_form)
#        rhs = dolfin.rhs(self.weak_form)
#
#        A = dolfin.assemble(lhs)
#
#        self.trial_func = dolfin.Function(self.FunctionSpace)
#        self.trial_func_initial_guess.assign(self.trial_init)
#
#        #  define file where to store data
#        ufile = dolfin.File(self.sim_pars['res_folder'] + '/' + self.sim_pars['pvd_name'] + '.pvd')
#        # Create empty time series
#        series = dolfin.TimeSeries(self.sim_pars['res_folder'] + '/' + self.sim_pars['TS_name'])
#
#        t = self.sim_pars['t0']
#        self.Nt = 0
#        while t < self.sim_pars['duration'] + dolfin.DOLFIN_EPS:
#            B = dolfin.assemble(rhs)
#            [bc.apply(A, B) for bc in self.boundaryConditions]
#            dolfin.solve(A, self.trial_func.vector(), B)
#
#            # Move to next time step
#            self.trial_func_initial_guess.assign(self.trial_func)
#
#            # Plot solution
#            dolfin.plot(self.trial_func, title="Temperature", rescale=rescale)
#
#            # Save to file
#            ufile << self.trial_func_initial_guess
#
#            series.store(self.mesh, t)
#            series.store(self.trial_func_initial_guess.vector(), t)
#
#            t += self.sim_pars['ts']
#            print "t = " + str(t)
#
#            self.Nt += 1
#
#        # Hold plot
#        if hold_plot:
#            dolfin.interactive()
#
#    def get_coordinates(self, name):
#        # return coordinates of mesh from its name
#        dofmap = self.FunctionSpaces['name'].dofmap()
#        coordinates = dofmap.tabulate_all_coordinates(self.meshes['name'])
#        return coordinates.reshape((-1, self.gdim))
#
#    def data_matrix(self):
#        snapshots = dolfin.TimeSeries(self.sim_pars['res_folder'] +
#                                      '/' +
#                                      self.sim_pars['TS_name'])
#        TS = list()
#        t = self.sim_pars['t0']
#        while t < self.sim_pars['duration'] + dolfin.DOLFIN_EPS:
#            x = dolfin.Vector()
#            snapshots.retrieve(x, t)
#            u = dolfin.Function(self.FunctionSpace)
#            u.vector()[:] = x.array()
#            TS.append(u)
#            t += self.sim_pars['ts']
#        data = np.zeros((self.Nx, self.Nt))
#        for ind_t, d in enumerate(TS):
#            data[:, ind_t] = d.vector().array()
#        return data
#
#
#    def imshow2Ddata(self, data, npoints=100):
#        """
#        plot the data on the point defined by self.mesh
#        return xi, yi, zi, with xi, yi the interpolation grid coordinates and zi the data over that grid
#        """
#        from plot_tools import imshow2Ddata
#        imshow2Ddata(data, self.coordinates, ninterp=20)
#
##
##displacement=dolf.Expression(("0.1*sin(x[0]*(2*3.14))","0.1*sin(x[1]*(2*3.14))"),element=V.ufl_element())
##        u=dolf.project(displacement,V=V)
##        dolf.plot(u)
##        mesh.move(u)
##        dolf.plot(mesh)
##
#
############################################################################
#
#	def gdim(self):
#		'''
#		base mesh geometry dimension
#
#			self.mesh.geometry().dim()
#
#		'''
#		return self.mesh.geometry().dim()

if __name__ == '__main__':
    fem = Fem(mesh(), vars_, pars)


