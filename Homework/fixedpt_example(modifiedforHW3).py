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



# import libraries
import numpy as np
from decimal import Decimal
    
def driver():

# test functions 
     f1 = lambda x: x*(1+(7-x**5)/x**2)**3
# fixed point is alpha1 = 7**1/5

     f2 = lambda x: x-(x**5-7)/x**2
#fixed point is alpha2 = 3.09... 
     f3 = lambda x: x-(x**5-7)/(5*x**4)

     f4 = lambda x: x-(x**5-7)/12

     Nmax = 1000
     tol = 1e-10

# # test f1 '''
#      x0 = 1
#      [xstar,ier] = fixedpt(f1,x0,tol,Nmax)
#      print('the approximate fixed point is:',xstar)
#      print('f1(xstar):',f1(xstar))
#      print('Error message reads:',ier)
    
# #test f2 '''
#      x0 = 1
#      [xstar,ier] = fixedpt(f2,x0,tol,Nmax)
#      print('the approximate fixed point is:',xstar)
#      print('f2(xstar):',f2(xstar))
#      print('Error message reads:',ier)

#test f3 '''
     x0 = 1
     [xstar,ier,count] = fixedpt(f3,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar)
     print('f2(xstar):',f3(xstar))
     print('Error message reads:',ier)
     print('Number of Required Iterations:', count)

#test f4 '''
     x0 = 1
     [xstar,ier,count] = fixedpt(f4,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar)
     print('f2(xstar):',f4(xstar))
     print('Error message reads:',ier)
     print('Number of Required Iterations:', count)

# define routines
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       #print('On iteration:', count)
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier,count]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier,count]
    

driver()
