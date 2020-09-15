#! /bin/bash/python3
import numpy as np # 
from numpy import tile
import matplotlib.pyplot as plt
#from numpy import linalg as LA
from scipy import linalg as LA
import math,cmath
pi = math.pi
cos = math.cos
sin = math.sin
#exp = cmath.exp
exp=np.exp
sqrt = math.sqrt
#log=math.log
log=np.log
from numpy.linalg import inv
from numpy import zeros
import lieb

# Now you can call defined function that module as follows
t1=1.0
t2=0.0
a=1.0
Nk=20
mass=0.0
dim=0.0
totomega=400
energy=-.01
#Ef=25.0
#temp=0.5
lieb.band(t1,t2,a,Nk,mass,dim,totomega,square=False)
#for temp in range(1,100):
#	g=square.fermi(energy , 0.0 , temp)
#	print temp, g


