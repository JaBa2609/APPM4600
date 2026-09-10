import numpy as np

def order_of_convergence(p_true, p_hat):
    """
    If the method converges with order alpha:
        e_{n+1} ~ C * e_n^alpha        with e_n = |p_hat[n] - p_true|

    To get alpha:
        alpha ~ ln(e_{n+1}/e_n) / ln(e_n/e_{n-1})
    """
    p_hat = np.asarray(p_hat, dtype=float)
    #print(p_hat)
    #print('shape of p_hat', np.shape(p_hat))
    #print(errors)
    errors = np.abs(p_hat - p_true)

    alphas = []
    for n in range(1, len(errors) - 2):
        e_prev, e_curr, e_next = errors[n - 1], errors[n], errors[n + 1]
        
        #print(e_next,e_curr,e_prev)
        alpha = np.log(e_next / e_curr) / np.log(e_curr / e_prev)
        alphas.append(alpha)
    lamb = abs(e_next)/abs(e_curr)**alphas[-1]
    #print(np.array(alphas))
    return [np.array(alphas[-1]), lamb]


"""
 This script explores the use of the fixed point method.  
 Two functions are considered that have different properties.
 I like to use this code before I talk about convergence analysis
 for the fixed point method as motivation.
"""

############################################# 
"""
Copyright (C) 2025  Adrianna M. Gillman

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
############################################# 

    
def driver():

# test functions 
    f1 = lambda x: np.sqrt(10/(x+4))

    Nmax = 100
    tol = (1/2)*(10**(-10)) #adjusted to have 10 digits of precision?

# test f1 '''
    x0 = [1.5]
    for x in x0:
        [xstar,ier,count,vec_approx] = fixedpt(f1,x,tol,Nmax)
        print('the approximate fixed point is:',xstar)
        print('f1(xstar):',f1(xstar))
        print('Error message reads:',ier)

    return vec_approx


# define routines
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''
    ''' vec_approx = vector of approximations '''
    vec_approx = [x0]

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       #print(x1, x0)
       vec_approx.append(x1)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          #print('The count is',count)
          return [xstar,ier,count,np.asarray(vec_approx)]
       x0 = x1
       #print(vec_approx)

    xstar = x1
    ier = 1
    return [xstar,ier,count,np.asarray(vec_approx)]
    
test=driver()
#print(test)
#print(order_of_convergence(1.3652300134140976,test))


# Section 3.1 Aitkens Del squared acceleration technique

def aitkens(p_input):
    p_n = lambda i: p_input[i-2] - (p_input[i-1]-p_input[i])**2/(p_input[i]-2*p_input[i-1]+p_input[i-2])
    p_approx = [p_input[0]]
    p_approx.append(p_n(1))
    print(p_approx)
    nMax = 100
    count = 0
    tol = 10**(-10)

    while count < nMax:
        count += 1
        p_approx.append(p_n(count))
        if abs( p_approx[count+1] - p_approx[count] ) < tol:
            return p_approx

    print("Aitken's output is", np.asarray(p_approx))
    return p_approx

aitkens(test)
order_of_convergence(1.3652300134140976, aitkens(test))