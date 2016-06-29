#!/usr/bin/env python

from larf import drift
import numpy
import matplotlib.pyplot as plt

def plot_steps2d(steps, grad):
    px, py = numpy.meshgrid(grad.domain[0], grad.domain[1], indexing='ij')
    fig,ax = plt.subplots()
    im = ax.pcolor(px, py, grad.scalar)
    fig.colorbar(im)
    ax.quiver(px, py, grad.components[0], grad.components[1])
    ptsX,ptsY = numpy.vstack(steps.points).T
    ax.scatter(ptsX, ptsY, s = 1000*numpy.asarray(steps.distance), c = numpy.asarray(steps.times), alpha=0.5)


def test_stepper_default():
    velo = drift.test_velocity()
    def velocity(notused, r):
        return velo(r)
    stepper = drift.Stepper(velocity, lcar=0.1)
    steps = stepper(0.0, [-0.5,1.0])
    print steps
    plot_steps2d(steps, velo)
    
    
def test_stepper_tweaked():
    velo = drift.test_velocity()
    def velocity(notused, r):
        return velo(r)
    stepper = drift.Stepper(velocity, lcar=0.1)
    visitor=drift.CollectSteps(stuck = drift.StuckDetection(distance = 0.001, nallowed = 2),
                               bounds = drift.BoundPrecision(prec = 1e-6, maxrelval = 2.0))
    steps = stepper(0.0, [-0.5,1.0], visitor)

    print steps
    plot_steps2d(steps, velo)
    
    
