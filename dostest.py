#! /bin/bash
from numpy.linalg import inv
from numpy import zeros
import dos 

# Now you can call defined function that module as follows
t1=1.0
t2=0.0
a=1.0
Nk=200
mass=0.0
dim=0.0*t1
totomega=800

dos.band(t1,t2,a,Nk,mass,dim,totomega,graphene=False)

