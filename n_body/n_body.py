import numpy
from numpy.random import rand
from numpy.linalg import solve
import h5py
from matplotlib import pyplot, gridspec

class Nbody(object):
    def __init__(self, nbodies=100, nt=1000, dt=0.0001, output_file="nbody_data.h5", integrator="dormand-prince"):
        self.G = 5.0 # Gravitational constant
        self.nbodies = nbodies # number of bodies
        self.nt = nt # number of timesteps
        self.dt = dt # size of timestep
        self.time = 0 # current time

        self.box_size = (10.,10.) # size of box, centred on origin

        self.q = rand(nbodies, 2, 2)
        self.q[:,0,:] = 1.0 * self.q[:,0,:] - 0.5 # velocities
        self.q[:,1,:] = 8. * self.q[:,1,:] - 4. # positions

        self.masses = 15. * (0.1 + rand(nbodies)) # masses of bodies

        self.U = self.potential_energy()
        self.T = self.kinetic_energy()
        self.I = self.moment_of_inertia()

        # set up hdf5 file with datasets
        self.output_file = h5py.File(output_file, "w")
        self.output_file.create_dataset("q", (nt+1, nbodies, 2, 2))
        self.output_file.create_dataset("U", (nt+1,))
        self.output_file.create_dataset("T", (nt+1,))
        self.output_file.create_dataset("I", (nt+1,))

        self.write_to_file()

        integrators = {"rk3": self.rk3, "dirk3":self.dirk3, "dormand-prince": self.dormand_prince, "adams-moulton": self.adams_moulton, "adams-bashforth": self.adams_bashforth}

        if integrator in integrators:
            self.integrator = integrators[integrator]
        else:
            print("Could not find integrator {} - defaulting to dormand-prince.".format(integrator))
            self.integrator = self.dormand_prince

    def step(self, q):
        # step simulation through one timestep
        f = numpy.zeros_like(q)

        for n in range(self.nbodies):
            r = q[numpy.newaxis, n, 1,:] - q[:,1,:]
            dist = numpy.sqrt(numpy.sum(r**2, axis=1))
            dist[numpy.absolute(dist) < 1.e-10] = 1.e12 # hard-core potential
            f[:,0,:] += self.G * self.masses[n] * r / dist[:,numpy.newaxis]**3

        f[:,1,:] = q[:,0,:]

        return f

    def rk3(self, f):
        # Third-order Runge-Kutta
        q1 = self.q + f(self.q) * self.dt

        q2 = 0.25 * (3. * self.q + q1 + f(q1) * self.dt)

        return (1./3.) * (self.q + 2. * q2 + 2. * f(q2) * self.dt)

    @staticmethod
    def isinteger(x):
        return numpy.equal(numpy.mod(x, 1), 0)

    def pexprb43(self, f):
        def phi(k, z):
            assert self.isinteger(k), "k = %f is not an integer" % k
            if k == 0:
                return numpy.exp(z)
            elif z < 1.e-12:
                return 1.0
            else:
                return (phi(k-1, z) - phi(k-1, 0)) / z
        Un2 = self.q + 0.5 * self.dt * phi(1,0.0) * f(self.q)
        Un3 = self.q + self.dt * phi(1,0.0) * f(self.q)

        Dn2 = f(Un2) - f(self.q)
        Dn3 = f(Un3) - f(self.q)

        return self.q + self.dt * phi(1, 0.0) * f(self.q) + self.dt * 16. * phi(3, 0.0) * Dn2 + self.dt * (-2. * phi(3, 0.0)) * Dn3

    def dirk3(self, f):
        # Implicit diagonal third order Runge-Kutta
        # decomposed flux into linear part (0, v)^T and non-linear part
        # (sum (Gm r / |r|^3), 0)^T
        mu = 0.5 * (1. - 1./numpy.sqrt(3.))
        nu = 0.5 * (numpy.sqrt(3.))
        gamma = 3. / (2. * (3. + numpy.sqrt(3.)))
        lbda = 3. * (1. + numpy.sqrt(3.)) / (2. * (3. + numpy.sqrt(3.)))

        f1 = f(self.q + self.dt * mu)[:,:,:]
        f1[:,1,:] = 0.
        f2 = f(self.q + self.dt * (nu + 2. * mu))[:,:,:]
        f2[:,1,:] = 0.
        q_new = numpy.zeros_like(self.q)

        for i in range(self.nbodies):
            A = numpy.array([[0., 0.], [1., 0.]])
            M = numpy.eye(2) - self.dt * mu * A
            b1 = self.q[i,:,:] + self.dt * mu * f1[i, :, :]
            y1 = solve(M, b1)

            b2 = y1 + self.dt * nu * (A * y1 + f1[i, :, :]) + self.dt * mu * f2[i, :, :]
            y2 = solve(M, b2)

            q_new[i,:,:] = (1 - lbda) * self.q[i,:,:] + lbda * y2 + self.dt * gamma * (A * y2 + f2[i,:,:])

        return q_new

    def dormand_prince(self, f):
        # see https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
        k = numpy.zeros((7 , self.nbodies, 2, 2))
        b = numpy.array([35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0.])

        k[0] = f(self.q)
        k[1] = f(self.q + self.dt * 1/5 * k[0])
        k[2] = f(self.q + self.dt * (3/40 * k[0] + 9/40 * k[1]))
        k[3] = f(self.q + self.dt * (44/45 * k[0] - 56/15 * k[1] + 32/9 * k[2]))
        k[4] = f(self.q + self.dt * (19372/6561 * k[0] - 25360/2187 * k[1] + 64448/6561 * k[2] - 212/729 * k[3]))
        k[5] = f(self.q + self.dt * (9017/3168 * k[0] - 355/33 * k[1] + 46732/5247 * k[2] + 49/176 * k[3] - 5103/18656 * k[4]))
        k[6] = f(self.q + self.dt * (35/384 * k[0] + 500/1113 * k[2] + 125/192 * k[3] - 2187/6784 * k[4] + 11/84 * k[5]))

        return self.q + self.dt * numpy.sum(k * b[:,numpy.newaxis, numpy.newaxis, numpy.newaxis], axis=0)

    def adams_bashforth(self, f, previous_fs):
        # See http://www.math.iit.edu/~fass/478578_Chapter_2.pdf
        # update previous_fs
        for i in range(4):
            previous_fs[i,:,:,:] = previous_fs[i+1,:,:,:]

        previous_fs[4,:,:,:] = f(self.q)
        return self.q + self.dt / 720. * (1901 * previous_fs[4,:,:,:]
                                          - 2774 * previous_fs[3,:,:,:]
                                          + 2616 * previous_fs[2,:,:,:]
                                          - 1274 * previous_fs[1,:,:,:]
                                          + 251 * previous_fs[0,:,:,:])

    def adams_moulton(self, f, previous_fs):
        # See http://www.math.iit.edu/~fass/478578_Chapter_2.pdf
        # predictor
        q_temp = self.adams_bashforth(f, previous_fs)

        # corrector
        previous_fs[0,:,:,:] = f(q_temp)
        q_new = self.q + self.dt / 720. * (251 * previous_fs[0,:,:,:]
                                          + 646 * previous_fs[4,:,:,:]
                                          - 264 * previous_fs[3,:,:,:]
                                          + 106 * previous_fs[2,:,:,:]
                                          - 19 * previous_fs[1,:,:,:])

        return q_new

    def potential_energy(self):
        U = 0
        for n in range(self.nbodies):
            r = self.q[numpy.newaxis, n, 1,:] - self.q[:,1,:]
            dist = numpy.sqrt(numpy.sum(r**2, axis=1))
            dist[numpy.absolute(dist) < 1.e-10] = 1.e12 # hard-core potential
            U += self.G * numpy.sum(self.masses[n]*self.masses[:] / dist)
        return U

    def kinetic_energy(self):
        p_squared = self.masses**2 * numpy.sum(self.q[:,1,:]**2, axis=1)

        return numpy.sum(p_squared / (2.0 * self.masses))

    def moment_of_inertia(self):
        return numpy.sum(self.masses * numpy.sum(self.q[:,1,:]**2, axis=1))

    def write_to_file(self):
        self.output_file["q"][self.time,:,:,:] = self.q
        self.output_file["U"][self.time] = self.U
        self.output_file["T"][self.time] = self.T
        self.output_file["I"][self.time] = self.I

    def contain_particles(self):
        # Stop particles moving outside of box
        for n in range(self.nbodies):
            # x-direction
            if self.q[n,1,0] < -0.5*self.box_size[0]:
                self.q[n,1,0] = -self.box_size[0] - self.q[n,1,0]
                self.q[n,0,0] *= -1.
            elif self.q[n,1,0] > 0.5*self.box_size[0]:
                self.q[n,1,0] = self.box_size[0] - self.q[n,1,0]
                self.q[n,0,0] *= -1.
                # y-direction
            if self.q[n,1,1] < -0.5*self.box_size[1]:
                self.q[n,1,1] = -self.box_size[1] - self.q[n,1,1]
                self.q[n,0,1] *= -1.
            elif self.q[n,1,1] > 0.5*self.box_size[1]:
                self.q[n,1,1] = self.box_size[1] - self.q[n,1,1]
                self.q[n,0,1] *= -1.

    def evolve(self):
        # run simulation
        fig = pyplot.figure(figsize=(18,8))
        gridspec.GridSpec(2,4)
        ax1 = pyplot.subplot2grid((2,4), (0,0), colspan=2, rowspan=2)
        ax2 = pyplot.subplot2grid((2,4), (0,2), colspan=2)
        ax3 = pyplot.subplot2grid((2,4), (1,2), colspan=2)
        pyplot.ion()

        previous_fs = numpy.zeros((5, self.nbodies, 2, 2))

        for i in range(self.nt):
            self.time += 1

            if self.integrator == self.adams_moulton or self.integrator == self.adams_bashforth:
                if self.time > 5:
                    self.q = self.integrator(self.step, previous_fs)
                else:
                    previous_fs[self.time-1,:,:,:] = self.step(self.q)
                    self.q = self.dormand_prince(self.step)
            else:
                self.q = self.integrator(self.step)

            self.contain_particles()

            self.U = self.potential_energy()
            self.T = self.kinetic_energy()
            self.I = self.moment_of_inertia()

            #self.print_state()
            self.write_to_file()
            self.plot_bodies(ax1)
            if i % 20 == 0:
                self.plot_energies(ax2)
                self.plot_inertia(ax3)

        self.output_file.close()

    def print_state(self):
        # print state to screen
        print("time = {},  U = {},  T = {},  I = {}".format(self.time, self.U, self.T, self.I))

    def plot_bodies(self, ax):
        ax.clear()
        ax.set_xlim(-self.box_size[0]/2., self.box_size[0]/2)
        ax.set_ylim(-self.box_size[1]/2.,self.box_size[1]/2.)
        ax.set_title('Particles')
        ax.scatter(self.q[:,1,0], self.q[:,1,1], s=20.*self.masses, c=numpy.array(range(self.nbodies))/self.nbodies, alpha=0.5)
        pyplot.pause(.00001)

    def plot_energies(self, ax):
        ax.set_xlim(0,self.nt)
        ax.set_ylim(0,5.e6)
        ax.set_title('Energy')
        ax.scatter(self.time, self.T + self.U, c='k', marker='o', label=r'$E$')
        ax.scatter(self.time, self.U, c='r', marker='x', label=r'$U$')
        ax.scatter(self.time, self.T, c='b', marker='+', label=r'$T$')
        if self.time == 1:
            ax.legend()

    def plot_inertia(self, ax):
        ax.set_xlim(0,self.nt)
        ax.set_ylim(0,3.e4)
        ax.set_title('Moment of inertia')
        ax.scatter(self.time, self.I, marker='x')

if __name__ == "__main__":
    nbody = Nbody(integrator="adams-moulton")
    nbody.evolve()
