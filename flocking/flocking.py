"""
Simulation based on flocking behaviour, as described here: https://en.wikipedia.org/wiki/Flocking_(behavior)#Flocking_rules

Behaviour is controlled by three simple rules:
    1. Separation - avoid crowding neighbors (short range repulsion)
    2. Alignment - steer towards average heading of neighbors
    3. Cohesion - steer towards average position of neighbors (long range attraction)
"""

import numpy
from matplotlib import pyplot, gridspec
from numpy.random import rand
from numpy.linalg import norm

class Bird(object):
    """
    Bird class. Bird has position and velocity.
    """
    def __init__(self, position, velocity):
        self.pos = numpy.array(position[:])
        self.velocity = numpy.array(velocity[:])

    def update_bird(self, new_pos, new_vel):
        self.pos[:] = new_pos[:]
        self.velocity[:] = new_vel[:]

class Simulation(object):

    def __init__(self, nbirds=50, nt=500):
        self.birds = set()

        self.n_dims = 2
        self.ax_lims = [0, 10., 0, 10.]
        self.nt = nt
        self.separation_distance = 2.0

        for i in range(nbirds):
            # randomly generate positions and velocities
            position = [rand() * (self.ax_lims[2*j+1] - self.ax_lims[2*j]) + self.ax_lims[2*j] for j in range(self.n_dims)]
            # assume all birds have a speed of 1
            velocity = numpy.array([rand() - 0.5 for j in range(self.n_dims)])
            velocity /= norm(velocity)

            self.birds.add(Bird(position, velocity))

        # initialise plotting
        fig = pyplot.figure(figsize=(16,10))
        gridspec.GridSpec(2,3)
        self.ax1 = pyplot.subplot2grid((2,3), (0,0), colspan=2, rowspan=2)
        self.ax2 = pyplot.subplot2grid((2,3), (0,2), colspan=1)
        self.ax3 = pyplot.subplot2grid((2,3), (1,2), colspan=1)
        fig.show()
        #pyplot.ion()

    def evolve_birds(self):
        nbirds = len(self.birds)
        average_position = numpy.average([b.pos for b in self.birds], axis=0)
        average_alignment = numpy.average([b.velocity for b in self.birds], axis=0)

        for b in self.birds:
            avoid_direction = numpy.zeros(self.n_dims)
            for c in self.birds - {b}:
                distance = numpy.sqrt(numpy.sum((b.pos - c.pos)**2))
                if distance < self.separation_distance:
                    avoid_direction += (c.pos - b.pos) / norm(c.pos - b.pos)

            if norm(avoid_direction) > 0.0:
                avoid_direction /= norm(avoid_direction)
                b.velocity = (average_position - b.pos) / norm(b.pos - average_position) + average_alignment / norm(average_alignment) - avoid_direction
            else:
                b.velocity = (average_position - b.pos) / norm(b.pos - average_position) + average_alignment / norm(average_alignment)

            b.velocity /= norm(b.velocity)

        # update positions
        for b in self.birds:
            b.pos += b.velocity

    def run(self):
        for t in range(self.nt):
            self.evolve_birds()
            self.print_birds(t)

    def print_birds(self, t):
        positions = numpy.array([b.pos for b in self.birds])
        velocities = numpy.array([b.velocity for b in self.birds])

        self.ax1.clear()
        self.ax2.clear()

        self.ax1.quiver(positions[:,0], positions[:,1], velocities[:,0], velocities[:,1])

        average_direction = numpy.average(velocities, axis=0)
        average_direction /= norm(average_direction)

        self.ax2.quiver(0, 0, average_direction[0], average_direction[1] ,units='x', scale=1/0.05)
        self.ax2.set_title('Average direction')

        if t % 10 == 0:
            average_position = numpy.average(positions, axis=0)

            self.ax3.quiver(average_position[0], average_position[1], average_direction[0], average_direction[1])
            self.ax3.set_title('Path')

        pyplot.pause(.00001)

if __name__ == "__main__":
    sim = Simulation()
    sim.run()
