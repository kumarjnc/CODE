#!/usr/bin/python

import numpy as np
from numpy import tile
import matplotlib.pyplot as plt
from numpy import linalg as LA
import math,cmath
pi = math.pi
cos = math.cos
sin = math.sin
exp = cmath.exp
sqrt = math.sqrt

#----------------------------------------------------------------------------------------------------------
def band(t1,t2,a,Nk,mass,dim,totomega,graphene=False):
    print t1,t2,a,Nk
    if graphene:
	assert (t1 > 0), 'Only positive hopping are allowed!'
    	print('Band structure plot of the Graphene with nearest and next nearest neighbour hopping......')
	sigma_z = np.array([[1,0],[0,-1]])
	gapplus=np.zeros((Nk,Nk))
	gapminus=np.zeros((Nk,Nk))
	f = open("graphene.dat", 'w')


	for nkx in range(1, Nk):
        	for nky in range(1,Nk):
    			kx = 3*pi/Nk*(nkx-Nk-1)
    			ky = 3*pi/Nk*(nky-Nk-1)
 			fk = -t1*exp(-1j*kx*a)*(1+2*exp(1j*1.5*kx*a)*cos(sqrt(3)/2*ky*a))
 			gk = -t2/2*(2*cos(sqrt(3)*ky*a)+4*cos(0.5*sqrt(3)*ky*a)*cos(1.5*kx*a))
 
 			Ek = np.array([[0.0, fk],[np.conjugate(fk), 0.0]])+ mass*sigma_z+gk*sigma_z

       			vk=LA.eigvals(Ek)
    			gapplus[nkx][nky]=vk[0].real
    			gapminus[nkx][nky]=vk[1].real
    			f.write("%5.2f %5.2f %5.2f %5.2f\n" % (kx,ky,gapplus[nkx][nky],gapminus[nkx][nky]))
		f.write("\n")                  
	f.close()
    else:
    
	print('Band structure plot of the Lieb lattice (a decorated lattice) with nearest and next nearest neighbour hopping......')

	assert ( t1 > 0), 'Only positive hoppin are allowed!'	
    	gapplus=np.zeros((Nk,Nk))
	gapflat=np.zeros((Nk,Nk))
	gapminus=np.zeros((Nk,Nk))

	f  = open('dispersion_lieb.dat','w');

	    for nkx in range(1, Nk):
		for nky in range(1, Nk):
			kx = 2*pi/Nk*nkx
			ky = 2*pi/Nk*nky
			fk = cos(0.5*a*kx)
			gk = cos(0.5*a*ky)
			fkk = 2*cos(0.5*a*kx)*cos(0.5*a*ky)
			gkk = 2*cos(0.5*a*kx)*cos(0.5*a*ky)
			fdim = 2*(1j)*sin(0.5*a*kx)
			gdim = 2*(1j)*sin(0.5*a*ky)
			Ek = -2.0*t1*np.array([[0.0, fk, gk],[fk, 0.0, 0.0],[gk, 0.0, 0.0]])-4.0*t2*np.array([[0.0, 0.0, 0.0],[0.0, 0.0, fkk],[0.0, gkk, 0.0]])+dim*np.array([[0.0, -fdim, -gdim],[fdim, 0.0, 0.0],[ gdim, 0.0, 0.0]])
	 
			vk=LA.eigvals(Ek)
			gapplus[nkx][nky]=vk[0].real
			gapflat[nkx][nky]=vk[1].real
			gapminus[nkx][nky]=vk[2].real

			f.write("%5.2f %5.2f %5.2f %5.2f %5.2f\n" % (kx,ky,gapplus[nkx][nky],gapflat[nkx][nky],gapminus[nkx][nky]))
		f.write("\n")
	    f.close()
	    f = open('dos.dat','w')
	    for Nomega in range(1,totomega):
		dosplus=0.0
		dosflat =0.0
		dosminus=0.0
		for nkx in range (1,Nk):
			for nky in range(1,Nk):
				omega=10.0/totomega*(Nomega-totomega/2)
				dosplus= dosplus+0.01/((omega-gapplus[nkx][nky])**2+(0.01)**2)
				dosflat= dosflat+0.01/((omega-gapflat[nkx][nky])**2+(0.01)**2)
				dosminus= dosminus+0.01/((omega-gapminus[nkx][nky])**2+(0.01)**2)
		dostot=dosflat+dosminus
		f.write("%6.2f %6.2f %6.2f %6.2f\n" %(omega,dosplus/(Nk*Nk),dosflat/(Nk*Nk),dostot/(Nk*Nk)))
    
