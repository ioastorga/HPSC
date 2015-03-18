
"""
Newton Method Module
Modified by: ** Io Astorga **
"""


import numpy as np
# import matplotlib.pyplot as plt

def solve(fval, x0, debug):
    """

    """
    x=x0
    maxiter=20
    toler=1.e-14
    if debug:
        print "Initial guess: x= %22.15e " % x
        print "Max iterations = %d" % maxiter

    iternum = 1

    while iternum <= maxiter:
        fx,fpx = fval(x)
        
        if abs(fx) < toler:
            break

        deltax = fx/fpx
        x = x - deltax

        iternum = iternum+1

        if debug:
            print "After %d iterations, x = %22.15e " % (iternum,x)

    """
    May not have converged
    """

    if iternum > maxiter:
        fx,fpx = fval(x)
        war_message = "Warning: has not yet converged"
        assert abs(fx) > tol, war_message
        
    iterations = iternum-1

    return x, iterations



#######

def fvals_sqrt(x):
    """
    Return f(x) and f'(x) for applying Newton to find a square root.
    """
    f = x**2 - 4.
    fp = 2.*x
    return f, fp

def test1(debug_solve=False):
    """
    Test Newton iteration for the square root with different initial
    conditions.
    """
    from numpy import sqrt
    for x0 in [1., 2., 100.]:
        print " "  # blank line
        x,iters = solve(fvals_sqrt, x0, debug=debug_solve)
        print "solve returns x = %22.15e after %i iterations " % (x,iters)
        fx,fpx = fvals_sqrt(x)
        print "the value of f(x) is %22.15e" % fx
        assert abs(x-2.) < 1e-14, "*** Unexpected result: x = %22.15e"  % x


if __name__=="__main__":
    # "main program"
    
    # the code below is executed only if the module is executed at the command line,
    #    $ python demo2.py
    # or run from within Python, e.g. in IPython with
    #    In[ ]:  run demo2
    # not if the module is imported.
    print "Running test..."
    test1()

