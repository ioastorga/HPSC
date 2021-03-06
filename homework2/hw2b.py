
"""
Demonstration module for quadratic interpolation.
Update this docstring to describe your code.
Modified by: ** Io Astorga **
"""


import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve

def quad_interp(xi,yi):
    """
    Quadratic interpolation.  Compute the coefficients of the polynomial
    interpolating the points (xi[i],yi[i]) for i = 0,1,2.
    Returns c, an array containing the coefficients of
      p(x) = c[0] + c[1]*x + c[2]*x**2.
    """

    # check inputs and print error message if not valid:

    error_message = "xi and yi should have type numpy.ndarray"
    assert (type(xi) is np.ndarray) and (type(yi) is np.ndarray), error_message

    error_message = "xi and yi should have length 3"
    assert len(xi)==3 and len(yi)==3, error_message

    # Set up linear system to interpolate through data points:
    A = np.vstack([np.ones(3), xi, xi**2]).T
    """
    A = np.array([[1.,xi[0],xi[0]**2 ],[1.,xi[1],xi[1]**2 ],[1.,xi[2],xi[2]**2]])
    "The np.vstack function takes as its argument a list of numpy arrays and 
    stacks them vertically into a matrix."
    """
    b = yi

    ### Fill in this part to compute c ###
    c=solve(A,b)

    return c


def test_quad1():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([-1.,  0.,  2.])
    yi = np.array([ 1., -1.,  7.])
    c = quad_interp(xi,yi)
    c_true = np.array([-1.,  0.,  2.])
    print "c =      ", c
    print "c_true = ", c_true
    # test that all elements have small error:
    assert np.allclose(c, c_true), \
        "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)

def test_quad2():
    c_true = np.array([-1., 1.5, 0.5])
    xi = np.array([1., 2., 3.])
    yi =np.array([1., 4., 8.])
    c = quad_interp(xi,yi)
    print "c= ", c
    print "c_true_ ", c_true
    assert np.allclose(c,c_true), \
            "Incorrect result, c = %s, but expected: c = %s" % (c, c_true)


def plot_quad(xi,yi):
    """
    Print quad and points
    """
    c = quad_interp(xi, yi)
    print "c =  " , c
    x = np.linspace(xi.min() - 1, xi.max() +1, 1000)
    y = c[0] +c[1]*x + c[2]*x**2

    plt.figure(1)
    plt.clf()
    plt.plot(x,y,'g-')

    plt.plot(xi, yi, 'bo')
    plt.ylim(-5, 10)

    plt.title("Data points and interpolating polynomial")
    plt.savefig('hwb2.png')


    """
    Ini cubic
    """

def cubic_interp(xi,yi):
    """
    Cubic interpolation.  Compute the coefficients of the polynomial
    interpolating the points (xi[i],yi[i]) for i = 0,1,2, 3.
    Returns c, an array containing the coefficients of
      p(x) = c[0] + c[1]*x + c[2]*x**2 +c[3]*x**3.
    """

    # check inputs and print error message if not valid:

    error_message = "xi and yi should have type numpy.ndarray"
    assert (type(xi) is np.ndarray) and (type(yi) is np.ndarray), error_message

    error_message = "xi and yi should have length 3"
    assert len(xi)==4 and len(yi)==4, error_message

    # Set up linear system to interpolate through data points:
    A = np.vstack([np.ones(4), xi, xi**2, xi**3]).T
    """
    A = np.array([[1.,xi[0],xi[0]**2 ],[1.,xi[1],xi[1]**2 ],[1.,xi[2],xi[2]**2]])
    "The np.vstack function takes as its argument a list of numpy arrays and 
    stacks them vertically into a matrix."
    """
    b = yi

    ### Fill in this part to compute c ###
    c=solve(A,b)

    return c


def test_cubic1():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([0.,  1.,  2., 3.])
    yi = np.array([ 1., 0.,  1., -2.])
    c = cubic_interp(xi,yi)
    c_true = np.array([1.,  -4.,  4., -1.])
    print "c =      ", c
    print "c_true = ", c_true
    # test that all elements have small error:
    assert np.allclose(c, c_true), \
        "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)

def plot_cubic(xi,yi):
    """
    Print quad and points
    """
    c = cubic_interp(xi, yi)
    print "c =  " , c
    x = np.linspace(xi.min() - 1, xi.max() +1, 1000)
    y = c[0] +c[1]*x + c[2]*x**2 +c[3]*x**3

    plt.figure(1)
    plt.clf()
    plt.plot(x,y,'g-')

    plt.plot(xi, yi, 'bo')
    plt.ylim(-5, 10)

    plt.title("Data points and interpolating polynomial")
    plt.savefig('hwb2cubic.png')


    """
    Ini poly
    """

def poly_interp(xi,yi):
    """
    Cubic interpolation.  Compute the coefficients of the polynomial
    interpolating the points (xi[i],yi[i]) for i = 0,1,2,3. . .
    Returns c, an array containing the coefficients of
      p(x) = c[0] + c[1]*x + c[2]*x**2 +c[3]*x**3 . . .
    """

    # check inputs and print error message if not valid:

    error_message = "xi and yi should have type numpy.ndarray"
    assert (type(xi) is np.ndarray) and (type(yi) is np.ndarray), error_message

    error_message = "xi and yi should have length 3"
    assert len(xi)==len(yi), error_message
    print "No. eq. =    ", len(xi)
    print "Pol. degree =    ", len(xi)-1
    
    # Set up linear system to interpolate through data points:
    A= np.ones(len(xi))
    for j in range(1,len(xi)):
        A=np.vstack([A,xi**j])
    #print "A before T =   ", A
    A=A.T
   # print "A after T =  ", A
    """
    A = np.vstack([np.ones(4), xi, xi**2, xi**3]).T
    A = np.array([[1.,xi[0],xi[0]**2 ],[1.,xi[1],xi[1]**2 ],[1.,xi[2],xi[2]**2]])
    "The np.vstack function takes as its argument a list of numpy arrays and 
    stacks them vertically into a matrix."
    """
    b = yi

    ### Fill in this part to compute c ###
    c=solve(A,b)
    return c


def test_poly1():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([0.,  1.,  2., 3., 4.])
    yi = np.array([ 1., 0.,  1., -2., -1])
    c = poly_interp(xi,yi)
    c_true = np.array([0.58333333, -4.5, 10.41666667, -7.5, 1.])
    print "c =      ", c
    print "c_true = ", c_true
    # test that all elements have small error:
    assert np.allclose(c, c_true), \
        "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)

def test_poly2():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([0.,  1.,  2., 3., 4.,5.])
    yi = np.array([ 1., 0.,  1., -2., -1, 0.])
    c = poly_interp(xi,yi)
    c_true = np.array([ -0.21666667, 2.75 ,-12.08333333, 21.25, -12.7, 1.])
    print "c =      ", c
    print "c_true = ", c_true
    # test that all elements have small error:
    assert np.allclose(c, c_true), \
        "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)

def plot_poly(xi,yi):
    """
    Print quad and points
    """
    c = poly_interp(xi, yi)
    print "c =  " , c
    x = np.linspace(xi.min() - 1, xi.max() +1, 1000)
    n=len(xi)
    y = c[n-1]
    for j in range(n-1, 0, -1):
        y = y*x + c[j-1]

    # y = c[0] +c[1]*x + c[2]*x**2 +c[3]*x**3

    plt.figure(1)
    plt.clf()
    plt.plot(x,y,'g-')

    plt.plot(xi, yi, 'bo')
    plt.ylim(-5, 10)

    plt.title("Data points and interpolating polynomial")
    plt.savefig('hwb2poly.png')

if __name__=="__main__":
    # "main program"
    # the code below is executed only if the module is executed at the command line,
    #    $ python demo2.py
    # or run from within Python, e.g. in IPython with
    #    In[ ]:  run demo2
    # not if the module is imported.
    print "Running test..."
    test_quad1()

