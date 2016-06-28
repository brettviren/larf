#!/usr/bin/env python

from larf import drift, units
import math
import numpy as np
import matplotlib.pyplot as plt

def test_gradient1d():
    x = np.linspace(-1,1,11)
    px, = np.meshgrid(x, indexing='ij') # this is actually same as x but keep notation same as Nd.
    pot = x**2
    g = drift.Gradient(pot, x)
    g(0); g(0.5); g(-.938727)
    try:
        g(1.1)
    except ValueError:
        print 'got value error as expected'
    else:
        raise AssertionError
        

def test_gradient2d():
    xy = np.linspace(-2,2,21), np.linspace(-4,4,21)
    px, py = np.meshgrid(*xy, indexing='ij')
    pot = px * np.exp(-px**2 - py**2)
    gradx, grady = np.gradient(pot)

    fig,ax = plt.subplots()
    im = ax.pcolor(px, py, pot)
    fig.colorbar(im)
    ax.quiver(px, py, gradx, grady)
    return fig,ax

def test_gradient3d():
    xyz = np.linspace(-2,2,21), np.linspace(-4,4,21), np.linspace(-10,10,21)
    px, py, pz = np.meshgrid(*xyz, indexing='ij')
    pot = px * np.exp(-px**2 - py**2)
    gradx, grady, gradz = np.gradient(pot)

    fig,ax = plt.subplots()
    zind = 10
    im = ax.pcolor(px[:,:,zind], py[:,:,zind], pot[:,:,zind])
    fig.colorbar(im)
    ax.quiver(px[:,:,zind], py[:,:,zind], gradx[:,:,zind], grady[:,:,zind])
    return fig,ax

def test_gradient_resample():
    xy = np.linspace(-4,4,21), np.linspace(-2,2,21)
    px, py = np.meshgrid(*xy, indexing='ij')
    pot = px * np.exp(-px**2 - py**2)

    # resample to finer scale
    gfun = drift.Gradient(pot, xy)

    X,Y = np.linspace(-4,4,201), np.linspace(-2,2,201)
    PX, PY = np.meshgrid(X, Y, indexing='ij')

    fine = gfun.regrid((X,Y))
    finemag = np.sqrt(fine[0]**2 + fine[1]**2)

    fig,ax = plt.subplots()
    im = ax.pcolor(PX, PY, finemag)
    fig.colorbar(im)
    ax.quiver(PX, PY, fine[0], fine[1])
    return fig,ax

def test_step():
    xlin, ylin = np.linspace(-4,4,21), np.linspace(-2,2,21)
    px, py = np.meshgrid(xlin, ylin, indexing='ij')
    pot = px * np.exp(-px**2 - py**2)
    velo = drift.Gradient(pot, (xlin, ylin))
    
    def velocity(r, t):
        return velo(r)

    points = list()
    point = np.asarray([-0.5,1.0])
    tnow = 0.0
    dt = 1
    for count in range(1000):
        points.append(point)
        v = velo(point)
        try:
            pnext = drift.step_rk4(point, tnow, tnow+dt, velocity)
        except ValueError:
            print 'last point: %s' % str(pnext)
            break
        #print "%f: %s --> %s @ %s" % (tnow, point, pnext, v)
        tnow += dt
        point = pnext

    ptsX, ptsY = np.vstack(points).T

    fig,ax = plt.subplots()
    im = ax.pcolor(px, py, pot)
    fig.colorbar(im)
    ax.quiver(px, py, velo.components[0], velo.components[1])
    ax.scatter(ptsX, ptsY)
    

def plot_pot(pot, r, figax=None):
    if not figax:
        figax = plt.subplots()
    fig,ax = figax
    u, v = np.gradient(pot)
    #im = ax.pcolormesh(r[0], r[1], pot)
    im = ax.pcolor(r[0], r[1], pot)
    fig.colorbar(im)
    ax.quiver(r[0], r[1], u, v)
    return fig,ax

def get_field():
    #x,y = np.linspace(-100,100,21), np.linspace(-100,100,21)
    xy = np.linspace(-2,2,21), np.linspace(-2,2,21)
    px, py = np.meshgrid(*xy, indexing='ij')
    #pot = np.sqrt(px**2 + py**2)
    pot = px * np.exp(-px**2 - py**2)
    return pot,xy

def get_velo():
    pot,XYlinspaces = get_field()
    Edrift = drift.field(pot)
    Emag = drift.mag(Edrift)
    mu = drift.mobility(Emag)
    return mu*Edrift, XYlinspaces
    

def find_index(xypoint, XYlinspaces):
    indices = list()
    for point,space in zip(xypoint, XYlinspaces):
        ind = np.searchsorted(space, point)
        indices.append(min(ind, len(space)-1))
    return tuple(indices)

def test_rk4():
    velo, XYlinspaces = get_velo()

    from scipy.interpolate import RegularGridInterpolator
    vx = RegularGridInterpolator(XYlinspaces, velo[0])
    vy = RegularGridInterpolator(XYlinspaces, velo[1])

    def velocity(r, t):
        v = np.asarray((vx(r)[0],vy(r)[0]))
        return v

    points = list()
    r = np.asarray([1.0,1.0])
    tnow = 0
    dt = 1e-5
    while True:
        v = velocity(r,tnow)
        rnext = drift.step_rk4(r, tnow, tnow+dt, velocity)
        print "%s -> %s at %f with %s" % (r,rnext,tnow, v)
        rnext2 = r[0]**2+r[1]**2
        if rnext2 > 10000:      # fixme: assume knowledge of domain from get_field()
            break
        points.append(r)
        r = rnext
        tnow += dt
    return np.asarray(points)


def test_grad():
    x, y = np.meshgrid(np.linspace(-2,2,21), np.linspace(-2,2,21), indexing='ij')
    z = x * np.exp(-x**2 - y**2)

    u, v = np.gradient(z, .2, .2)

    fig, ax = plt.subplots()
    im = ax.pcolormesh(x,y,z)
    fig.colorbar(im)
    ax.quiver(x, y, u, v)
    plt.show()    


    
