"""
Some details here
"""

def sqrt2(x, debug=False):
    """
    by default debug is false
    more details
    """
    s=1.
    kmax = 100
    tol = 1.e-14
    for k in range(kmax):
        if debug:
            print "Before iteration %s, s = %s" % (k,s) 
        s0 = s
        s=0.5*(s+x/s)
        delta_s = s - s0
        if abs(delta_s / x) < tol:
            break
    if debug:
        print "After %s iterations, s= %s" % (k+1,s)
    return s
