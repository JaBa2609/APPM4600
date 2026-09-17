import scipy.special as sci
import numpy as np
import matplotlib.pyplot as plt

#constants
Ti = 20
Ts = -15
alpha = 0.138*10**(-6)
t = 5.184*10**(6)

#lambda funcs
#f(x) = lambda x: Ts + (Ti - Ts)*sci.erf(x/2*np.sqrt(alpha*t))
#fp(x) = lambda x: (Ti - Ts) 2/np.sqrt(np.pi) * np.exp(-(x**2/4*alpha*t))

#Goal:
#   find x after 60 days
#   let t = 60 days, t = 5.184e+6 seconds

#xval = 2*np.sqrt(alpha*t)*sci.erfinv(-Ts/(Ti-Ts))
#print('depth of freezing after 60 days =', xval)

def driver():

# use routines    
    f = lambda x: Ts + (Ti - Ts)*sci.erf(x/(2*np.sqrt(alpha*t)))
    fp = lambda x: (Ti - Ts)*2/np.sqrt(np.pi) * np.exp(-(x**2/4*alpha*t))
    a = 0
    b = 3
    p0 = 0.5

    Nmax = 100
    tol = 1.e-13

    (p,pstar,info,it) = newton(f,fp,p0,tol, Nmax)
    print('the approximate root is', '%16.16e' % pstar)
    print('the error message reads:', '%d' % info)
    print('Number of iterations:', '%d' % it)

    #x_range = np.arange(0,15,.1)
    #y_range = f(x_range)
    #plt.plot(x_range, y_range)
    #plt.show()
    

# # define routines
# def bisection(f,a,b,tol):
    
# #    Inputs:
# #     f,a,b       - function and endpoints of initial interval
# #      tol  - bisection stops when interval length < tol

# #    Returns:
# #      astar - approximation of root
# #      ier   - error message
# #            - ier = 1 => Failed
# #            - ier = 0 == success

# #     first verify there is a root we can find in the interval 

#     fa = f(a)
#     fb = f(b);
#     if (fa*fb>0):
#        ier = 1
#        astar = a
#        return [astar, ier]

# #   verify end points are not a root 
#     if (fa == 0):
#       astar = a
#       ier =0
#       return [astar, ier]

#     if (fb ==0):
#       astar = b
#       ier = 0
#       return [astar, ier]

#     count = 0
#     d = 0.5*(a+b)
#     while (abs(d-a)> tol):
#       fd = f(d)
#       if (fd ==0):
#         astar = d
#         ier = 0
#         return [astar, ier]
#       if (fa*fd<0):
#          b = d
#       else: 
#         a = d
#         fa = fd
#       d = 0.5*(a+b)
#       count = count +1
# #      print('abs(d-a) = ', abs(d-a))
      
#     astar = d
#     ier = 0
#     print('count = ', count)
#     return [astar, ier]


def newton(f,fp,p0,tol,Nmax):
  """
  Newton iteration.
  
  Inputs:
    f,fp - function and derivative
    p0   - initial guess for root
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
     
  """
  p = np.zeros(Nmax+1);
  p[0] = p0
  for it in range(Nmax):
      p1 = p0-f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p,pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]

driver()

