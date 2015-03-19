# -*- coding: iso-8859-1 -*-

from numpy import *
from scipy import integrate
import matplotlib.pyplot as pl
from scipy.integrate import simps

import time

#################################################################
#1D finite element method on the interval I = [0,1] for the equation
#-u''(x) = f(x), u(0) = u(1) = 0, 
#################################################################

#################################################################
#19032015, Wilhelm Braun, wilhelm.braun@cantab.net
#################################################################

#define RHS of equation 
def f(x):
  #first two cases are giving rise to exact solutions below
  #return pi**(2)*sin(pi*x)
  #return 4.*pi**(2)*cos(pi*x)*sin(pi*x)
  return sin(10.*pi*x)*cos(25.*pi*x)
  #return 1.
#uniform partition of the interval I
M = 1000

partition = linspace(0, 1, num = M+2, endpoint = True)

#print partition, len(partition)

#construction of the M basis functions of basis functions
#These are piecewise linear functions, phi(j,xj) = 1

hj = partition[1]-partition[0]

def phi(j,x):

#evaluation of basis function at a point

###############################
##definition of support points
  x_minus = partition[j-1]
  xj = partition[j]
  x_plus = partition[j+1]
###############################

  if(x == xj):
    phi_point =1.
  elif(x > x_minus and x < xj):
    phi_point = (1./hj)*x
  elif (x > xj and x < x_plus):
    phi_point = (-1/hj) *x
  else:
    phi_point = 0
  return phi_point

def basis_function(j,x):
#evaluation over a whole range
 #print x, len(x)
 bf = zeros(len(x))
 for index in arange(0,len(x)):
  bf[index] = phi(j,x[index])
  #print index, out[index]

 return bf

#print len(vectorized_phi(1,partition))
#pl.plot(partition, basis_function(3, partition))
#pl.plot(partition, basis_function(10, partition))
#pl.show()

def load_function(j,x):
 return f(x) * basis_function(j,x)


#testing
for plot_index in arange(0, 500):
  pl.plot(partition, load_function(plot_index, partition))

pl.show()
#print integrate.quad(, 0.0, 1)[0]


def assemble_load_vector(interval):
 print "assembling load vector"
 start = time.clock()
 load_vector = zeros(M)
 for index in arange(0, M):
  # integration procedure
  load_vector[index] = simps(load_function(index+1, partition), partition)
  #print index, load_vector[index] 
 
 end = time.clock()
 print "done in ", end-start, "seconds"
 return load_vector

#assemble the load vector
b = assemble_load_vector(partition)
print M == len(b)

pl.plot(b, 'x')
pl.show()


#assemble stiffness matrix
start = time.clock()

print "assembling stiffness matrix for the problem"
A = zeros((M,M), dtype = 'float')
for i in range(M):
  A[i,i] = 2.
  if(i != 0):
    A[i,i-1] = -1
  if( i!= M-1):
    A[i, i+1] = -1

end = time.clock()

print "done in ", end-start, "seconds"
#divide by discretisation step
A = A/hj

#solve system of equations with numpy linear algebra package
solution_vector = linalg.solve(A, b)

print len(solution_vector)
print len(partition[1:M])
#def exact_solution(x):
 ##return sin(pi*x)
 #return sin(pi*x)*cos(pi*x)

pl.plot(partition[1:M+1], solution_vector)
#pl.plot(partition, exact_solution(partition), '-x')
pl.show()

#pl.plot(partition[1:M+1], (solution_vector - exact_solution(partition[1:M+1])))
#pl.show()

