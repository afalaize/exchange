# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 17:01:29 2016

@author: afalaize
"""
import numpy as np


class POD_SYS(object):

    def __init__(self, DOLFIN_SYS):
        self.DOLFIN_SYS = DOLFIN_SYS
        self.DOLFIN_SYS.run()
        self.D = self.DOLFIN_SYS.data_matrix()


    def split_MEANFLUCT(self, set_mean_to_zero=False):
        # Compute the mean field (can be set to zero)
        if set_mean_to_zero:
            self.D_mean = 0*np.mean(self.D, axis=1)
        else:
            self.D_mean = np.mean(self.D, axis=1)
        self.D_fluc = np.zeros(self.D.shape)
        for i, data in enumerate(self.D.T):
            self.D_fluc[:, i] = data - self.D_mean


    def compute_POD_basis(self, weighting_matrix="Id"):
        """
        Return the POD basis from the snapshots
        """

        self.Nx, self.Nt = self.D.shape

        if weighting_matrix=="Id":
            W = np.eye(self.Nx)
        elif weighting_matrix=="A":
            W = self.DOLFIN_SYS.assembly_matrices['A']

        # Temporal correlation matrix
        self.C = np.dot(self.D_fluc.T, np.dot(W, self.D_fluc))

        # Temporal correlation matrix eigen decomposition
        self.C_eigen_vals, self.C_eigen_vecs = np.linalg.eig(self.C)

        # Remove the imaginary part (which should be numerically close to zero)
        self.C_eigen_vals = np.real(self.C_eigen_vals)
        self.C_eigen_vecs = np.real(self.C_eigen_vecs)

        # Define POD basis
        self.POD_basis = np.dot(self.D_fluc, self.C_eigen_vecs)
        self.Npod = self.POD_basis.shape[1]

    def truncate_POD(self, threshold=1e-6):
        self.threshold = threshold
        self.POD_modes_energy = list()
        for i, val in enumerate(self.C_eigen_vals):
            mode_energy = sum(self.C_eigen_vals[:i+1])/sum(self.C_eigen_vals)
            self.POD_modes_energy.append(mode_energy)
        self.Npod = [me >= 1-self.threshold for me in self.POD_modes_energy].index(True)+1
        self.POD_basis = self.POD_basis[:, :self.Npod]


    def check_orthonormal_basis(self):
        print "Phi.T . Phi = \n", np.dot(self.POD_basis.T, self.POD_basis)


    def normalize_basis(self):
        for i in range(self.Npod):
            scalarprod = np.dot(self.POD_basis[:,i], self.POD_basis[:,i])
            self.POD_basis[:,i] = self.POD_basis[:,i]/np.sqrt(scalarprod)


    def plotPODmodes(self, nmax=9):
        from plot_tools import imshow2Ddata_ax
        import matplotlib.pyplot as plt
        nmax = np.min((nmax, self.Npod))
        n_axe0 = np.ceil(np.sqrt(nmax))
        n_axe1 = np.ceil(nmax/n_axe0)
        fig = plt.figure()
        for n in range(nmax):
            ax = fig.add_subplot(n_axe0, n_axe1, n+1)
            imshow2Ddata_ax(self.POD_basis[:,n], self.DOLFIN_SYS.coordinates, ax)
            ax.set_title('POD mode '+ str(n+1))


    def compute_gradient(self, data, Ninterp=100):

        def interp_2D_data(values, Ninterp=Ninterp):
            """
            values: 1D np.array with values on the mesh points
            mesh: a dolfin mesh object
            Ninterp: int, specifies the number of interpolation grid points
            """
            from scipy.interpolate import griddata
            X = np.array([el[0] for el in self.DOLFIN_SYS.coordinates])
            Y = np.array([el[1] for el in self.DOLFIN_SYS.coordinates])
            xi = np.linspace(X.min(), X.max(), Ninterp)
            yi = np.linspace(Y.min(), Y.max(), Ninterp)
            zi = griddata((X, Y), values, (xi[None, :], yi[:, None]), method='cubic')
            return xi, yi, zi

        def discrete_2D_spatial_gradient(xi, yi, zi):
            """
            """
            def grad2D(i, j):
                if i==0:
                    grad_x = (zi[i+1, j]-zi[i, j])/float(xi[i+1]-xi[i])
                elif i==nx-1:
                    grad_x = (zi[i, j]-zi[i-1, j])/float(xi[i]-xi[i-1])
                else:
                    grad_x = (zi[i+1, j]-zi[i-1, j])/float(xi[i+1]-xi[i-1])

                if j==0:
                    grad_y = (zi[i, j+1]-zi[i, j])/float(yi[j+1]-yi[j])
                elif j==ny-1:
                    grad_y = (zi[i, j]-zi[i, j-1])/float(yi[j]-yi[j-1])
                else:
                    grad_y = (zi[i, j+1]-zi[i, j-1])/float(yi[j+1]-yi[j-1])

                return np.array((grad_x, grad_y))
            nx, ny = zi.shape
            gradient = np.ndarray((nx, ny, 2))
            for i in range(nx):
                for j in range(ny):
                    gradient[i, j, :] = grad2D(i, j)
            return gradient

        def evaluate(xi, yi, zi):
            """
            return the evaluation of zi on the mesh; zi is defined on the grid (xi, yi)
            """
            from scipy.interpolate import bisplrep, bisplev
            nx, ny = zi.shape
            coord = list()
            for i in range(nx):
                for j in range(ny):
                    coord.append((xi[i], yi[j]))

            x_data = [c[0] for c in coord]
            y_data = [c[1] for c in coord]
            z_data = zi.reshape((1, nx*ny), order='F')
            tck = bisplrep(x_data, y_data, z_data)

            X = np.array([el[0] for el in self.DOLFIN_SYS.coordinates])
            Y = np.array([el[1] for el in self.DOLFIN_SYS.coordinates])
            Z = np.zeros(len(X))
            i = 0
            for x, y in zip(X, Y):
                Z[i] = bisplev(x, y, tck)
                i += 1
            return X, Y, Z

        xi, yi, zi = interp_2D_data(data)
        grad = discrete_2D_spatial_gradient(xi, yi, zi)
        X, Y, Grad_x = evaluate(xi, yi, grad[:, :, 0])
        X, Y, Grad_y = evaluate(xi, yi, grad[:, :, 1])
        return Grad_y, Grad_x


    def compute_basis_gradient(self, Ninterp=500):
        self.POD_basis_gradient = {'x1': np.zeros(self.POD_basis.shape),
                                   'x2': np.zeros(self.POD_basis.shape)}

        for n, pod in enumerate(self.POD_basis.T):
            Grad_x, Grad_y = self.compute_gradient(pod, Ninterp=Ninterp)
            self.POD_basis_gradient['x1'][:,n] = Grad_x
            self.POD_basis_gradient['x2'][:,n] = Grad_y



if __name__ == '__main__':
    from DOLFIN_SYS import DOLFIN_SYS
    pod = POD_SYS(DOLFIN_SYS())
    pod.split_MEANFLUCT(set_mean_to_zero=False)
    pod.compute_POD_basis(weighting_matrix="Id")
    pod.truncate_POD(threshold=1e-3)
    pod.check_orthonormal_basis()
    pod.normalize_basis()
    pod.check_orthonormal_basis()
    pod.plotPODmodes()
