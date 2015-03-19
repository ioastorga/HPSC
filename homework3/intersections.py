

"""
Show intersection of two curves using Newton  method
Modified by: ** Io Astorga **
"""


from numpy import *
from newton import solve
import matplotlib.pyplot as plt

def g1val(x):
    g1 = x*cos(pi*x)
    g1p = cos(pi*x) - x*pi*sin(pi*x)
    return g1, g1p

def g2val(x):
    g2 = 1. - 0.6*x**2
    g2p = -1.2*x
    return g2,g2p

def f_val(x):
    g1,g1p = g1val(x)
    g2,g2p = g2val(x)
    fval = g1 - g2
    fpval = g1p - g2p
    return fval, fpval


def intersectionG(x0input):
    intersec = []
    xinter = []
    for x0 in x0input:
        x,iters = solve(f_val, x0, debug = False)
        print "solve returns x= %22.15e after %i iterations" % (x,iters )
        fx,fpx = f_val(x)
        g1vx,g1vx_p = g1val(x)
        print "The value of intersection is %22.15e" % fx
        if fx not in intersec:
            intersec.append(g1vx)
            xinter.append(x)

    xplot = linspace(-5.,5.,100.)
    g1 = xplot*cos(pi*xplot)
    g2 = 1.-0.6*xplot**2
    plt.plot(xplot,g1, 'b', xplot, g2,'g',xinter,intersec,'ro')
    plt.show()


