#!/usr/bin/python

from __future__ import division
from mpmath import *

mp.dps = 100

def prettyprint(x, length=40):
    s = str(x)
    return s[:length]


def lnLambda(l1,l2,m):
    return log(sqrt( (2*l1+1)*(2*l2+1)*factorial(l1-m)*factorial(l2-m) /( factorial(l1+m)*factorial(l2+m)*l1*(l1+1)*l2*(l2+1) ) ))



print prettyprint(lnLambda(1,1,0))
print prettyprint(lnLambda(1,1,1))

print prettyprint(lnLambda(2,1,1))
print prettyprint(lnLambda(4,5,3))
print prettyprint(lnLambda(5,7,3))
print prettyprint(lnLambda(16,6,4))

print prettyprint(lnLambda(50,50,0))
print prettyprint(lnLambda(50,50,1))
print prettyprint(lnLambda(50,50,50))
print prettyprint(lnLambda(50,20,10))

print prettyprint(lnLambda(100,1,0))
print prettyprint(lnLambda(100,1,1))
print prettyprint(lnLambda(100,6,4))
print prettyprint(lnLambda(100,33,17))
print prettyprint(lnLambda(100,100,0))
print prettyprint(lnLambda(100,100,50))
print prettyprint(lnLambda(100,201,10))

print prettyprint(lnLambda(200,200,0))
print prettyprint(lnLambda(200,100,70))
print prettyprint(lnLambda(500,500,0))

print prettyprint(lnLambda(1000,1000,0))
print prettyprint(lnLambda(1000,1000,1))
print prettyprint(lnLambda(1000,1000,2))
print prettyprint(lnLambda(1000,1000,3))
print prettyprint(lnLambda(1000,1000,10))
print prettyprint(lnLambda(1000,1000,20))
print prettyprint(lnLambda(1000,1000,50))
print prettyprint(lnLambda(1000,1000,100))
print prettyprint(lnLambda(1000,1000,499))
print prettyprint(lnLambda(1000,1000,500))
print prettyprint(lnLambda(1000,1000,501))
print prettyprint(lnLambda(1000,1000,999))
print prettyprint(lnLambda(1000,1000,1000))
