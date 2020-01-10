from numpy import pi
from sympy import *
r = Symbol('r')
h = Symbol('h')

# f = 15/(pi*h**6)*(h-r)**3

f = simplify(15/(2*pi*h**3)*(-r**3/(2*h**3)+r**2/h**2+h/(2*r)-1))

print(f)

h_diff = f.diff(h)
r_diff = f.diff(r)

print(h_diff+r_diff)

f2 = h_diff+r_diff

h_diff = f2.diff(h)
r_diff = f2.diff(r)

print(simplify(h_diff+r_diff))