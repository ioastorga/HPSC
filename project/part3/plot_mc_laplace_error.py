
from pylab import *

# read in three columns from file and unpack into 3 arrays:
n,pde_approx,error = loadtxt('mc_Laplace_error.txt',unpack=True)

figure(1)
clf()
loglog(n,error,'-o',label='Monte-Carlo')
loglog([1,1e7],[1,sqrt(1e-7)],'k',label='1 / sqrt(N)')
legend()
xlabel('number of MC points used')
ylabel('abs(error)')
title('Log-log plot of relative error in MC Laplace')
savefig('mc_laplace_error.png')
