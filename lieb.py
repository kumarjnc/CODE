#!/usr/bin/python3
import numpy as np
from numpy import tile
#import matplotlib.pyplot as plt
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

#----------------------------------------------------------------------------------------------------------
def band(t1,t2,a,Nk,mass,dim,totomega,square=True):
	print (t1,t2,a,Nk)
	if square:
		assert (t1 > 0), "Only positive hopping are allowed!"
		print('Band structure plot of the square lattice with nearest and next nearest neighbour hopping......')
		gap1=np.zeros((Nk,Nk))
		gap2=np.zeros((Nk,Nk))
		gap3=np.zeros((Nk,Nk))
		gap4=np.zeros((Nk,Nk))
		  
		f = open("square.dat", 'w')
		
		alpha=0.0
		for nkx in range(0, Nk):
				for nky in range(0,Nk):
						kx = (2*pi)/(Nk)*(nkx)
						ky = (2*pi)/(Nk)*(nky)
						fk = 1.0+exp(-1j*kx*a)
						gk = 1.0+exp(-1j*ky*a)
						
						Ek = np.array([[0.0,(1+alpha)*t1*fk,(1+alpha)*t1*gk,0],[(1+alpha)*t1*np.conjugate(fk),0.0,0,(1-alpha)*t1*gk],[(1+alpha)*t1*np.conjugate(gk),0,0.0,(1-alpha)*t1*fk],[0,(1-alpha)*t1*np.conjugate(gk),(1-alpha)*t1*np.conjugate(fk),0.0]])


						
						vk=LA.eigvals(Ek)
						vk=sorted(vk)
						

						gap1[nkx][nky]=vk[0].real
						gap2[nkx][nky]=vk[1].real
						gap3[nkx][nky]=vk[2].real
						gap4[nkx][nky]=vk[3].real
						
						f.write("%5.2f %5.2f %5.2f %5.2f %5.2f\n" % (kx, vk[0].real, vk[1].real, vk[2].real, vk[3].real))
		f.close()
		dos1=np.zeros(totomega)
		dos2=np.zeros(totomega)
		dos3=np.zeros(totomega)
		dos4=np.zeros(totomega)
		f = open('dos.dat','w')
		for Nomega in range(1,totomega):
				omega=15.0/totomega*(Nomega-totomega/2)
				dossum1 =0.0
				dossum2 =0.0
				dossum3 =0.0
				dossum4 =0.0

				for nkx in range (0,Nk):
						for nky in range(0,Nk):
								z=omega+0.04j
								A1=-(1.0/pi)*(1.0/(z-gap1[nkx][nky])).imag
								dossum1=dossum1+A1
								A2=-(1.0/pi)*(1.0/(z-gap2[nkx][nky])).imag
								dossum2=dossum2+A2
								A3=-(1.0/pi)*(1.0/(z-gap3[nkx][nky])).imag
								dossum3=dossum3+A3
								A4=-(1.0/pi)*(1.0/(z-gap4[nkx][nky])).imag
								dossum4=dossum4+A4
				dos1[Nomega]=dossum1	
				dos2[Nomega]=dossum2
				dos3[Nomega]=dossum3
				dos4[Nomega]=dossum4
				f.write("%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n" %(omega,dos1[Nomega]/(Nk*Nk),dos2[Nomega]/(Nk*Nk),dos3[Nomega]/(Nk*Nk),dos4[Nomega]/(Nk*Nk), dos1[Nomega]/(Nk*Nk)+dos2[Nomega]/(Nk*Nk)+dos3[Nomega]/(Nk*Nk)+dos4[Nomega]/(Nk*Nk)))
#-----------------------------------------------------------norm of the dos-------------------------------------------------------------------
		norm=0.0
		for Nomega in range(2,totomega):
				omega=15.0/totomega*(Nomega-totomega/2)
				norm=norm+dos1[Nomega]/(Nk*Nk)+dos2[Nomega]/(Nk*Nk)+dos3[Nomega]/(Nk*Nk)+dos4[Nomega]/(Nk*Nk)
		totnorm=norm*(15.0/totomega)
		print (totnorm)
#----------------------------------------------------------------------------------------------------------------------------------------------	
		entropy = open('entropy.dat','w')
		for i in range(1,400):
				temp=(2.0/100)*i
				sum_entropy=0.0	
				for Nomega in range(1,totomega):
						omega=15.0/totomega*(Nomega-totomega/2)
						s=fermi(omega,0.0,temp)
						#sum_entropy=sum_entropy-2.0*(s*log(s)+(1.0001-s)*log(1.0001-s))*exp(-0.5*(omega/0.01)**2).real/(sqrt(2.0*pi)*0.01) #gaussian 
						sum_entropy=sum_entropy-2.0*(s*log(s))*(dos1[Nomega]+dos2[Nomega]+dos3[Nomega]+dos4[Nomega])*1.0/(Nk*Nk)
						#print omega, 2.0*s*log(s)
				sum_entropy_tot=sum_entropy*(15.0/totomega)				  
				entropy.write("%6.2f %6.2f\n" % (temp, sum_entropy_tot))
#------------------------------------------------------------------------------------------------------------------------------------------------
	else:			
		print('Band structure plot of the Lieb lattice (a decorated lattice) with nearest and next nearest neighbour hopping......')

		assert ( t1 > 0), 'Only positive hoppin are allowed!'	
		gapplus=np.zeros((Nk,Nk))
		gapflat=np.zeros((Nk,Nk))
		gapminus=np.zeros((Nk,Nk))

		f = open('dispersion_lieb.dat','w');

		for nkx in range(0, Nk):
				for nky in range(0, Nk):
						kx = 2*pi/Nk*nkx
						ky = 2*pi/Nk*nky
						fk = cos(0.5*a*kx)
						gk = cos(0.5*a*ky)
						fkk = 2*cos(0.5*a*kx)*cos(0.5*a*ky)
						gkk = 2*cos(0.5*a*kx)*cos(0.5*a*ky)
						fdim = 2*(1j)*sin(0.5*a*kx)
						gdim = 2*(1j)*sin(0.5*a*ky)
						#Ek = -2.0*t1*np.array([[0.0, fk, gk],[fk, 0.0, 0.0],[gk, 0.0, 0.0]])-4.0*t2*np.array([[0.0, 0.0, 0.0],[0.0, 0.0, fkk],[0.0, gkk, 0.0]])+dim*np.array([[0.0, -fdim, -gdim],[fdim, 0.0, 0.0],[ gdim, 0.0, 0.0]])
 
						#Ek = -2.0*t1*np.array([[0.0, fk, 0.0],[fk, 0.0, 0.0],[0.0, 0.0, 0.0]])-4.0*t2*np.array([[0.0, 0.0, 0.0],[0.0, 0.0, fkk],[0.0, gkk, 0.0]])+dim*np.array([[0.0, -fdim, -gdim],[fdim, 0.0, 0.0],[ gdim, 0.0, 0.0]])
						Ek = -2.0*t1*np.array([[0.0, fk, gk],[fk, 0.0, 0.0],[gk, 0.0, 0.0]])-4.0*t2*np.array([[0.0, 0.0, 0.0],[0.0, 0.0, fkk],[0.0, gkk, 0.0]])+dim*np.array([[0.0, -fdim, 0.0],[fdim, 0.0, 0.0],[0.0, 0.0, 0.0]])
						vk=LA.eigvals(Ek)
						gapplus[nkx][nky]=vk[0].real
						gapflat[nkx][nky]=vk[1].real
						gapminus[nkx][nky]=vk[2].real

						f.write("%5.2f %5.2f %5.2f %5.2f %5.2f\n" % (kx,ky,gapplus[nkx][nky],gapflat[nkx][nky],gapminus[nkx][nky]))
				f.write("\n")
		f.close()

		f = open('dos.dat','w')
		n=0.0
		n1=0.0
		for Nomega in range(1,totomega):
				omega=10.0/totomega*(Nomega-totomega/2)
				dosplus=0.0
				dosflat =0.0
				dosminus=0.0
				for nkx in range (0,Nk):
						for nky in range(0,Nk):	
								z=omega+0.02j
								A1=-(1.0/pi)*(1.0/(z-gapplus[nkx][nky])).imag
								dosplus=dosplus+A1
								A2=-(1.0/pi)*(1.0/(z-gapflat[nkx][nky])).imag
								dosflat=dosflat+A2
								A3=-(1.0/pi)*(1.0/(z-gapminus[nkx][nky])).imag
								dosminus=dosminus+A3

				dostot=dosflat+dosminus+dosplus
				n=n+(dosplus/(Nk*Nk))
				n1=n1+(dostot/(Nk*Nk))
				f.write("%6.2f %6.2f %6.2f %6.2f %6.2f\n" %(omega,dosplus/(Nk*Nk),dosflat/(Nk*Nk),dosminus/(Nk*Nk),dosplus/(Nk*Nk)+dosflat/(Nk*Nk)+dosminus/(Nk*Nk)))
		print  (n1*(10.0/totomega))


		entropy = open('entropy.dat','w')
		for i in range(1,60):
				temp=(1.0/10)*i
				sum_entropy=0.0
				for Nomega in range(1,totomega):
						omega=12.0/totomega*(Nomega-totomega/2)
						dosplus=0.0
						dosflat =0.0
						dosminus=0.0
						for nkx in range (0,Nk):
								for nky in range(0,Nk):
										z=omega+0.02j
										A1=-(1.0/pi)*(1.0/(z-gapplus[nkx][nky])).imag
										dosplus=dosplus+A1
										A2=-(1.0/pi)*(1.0/(z-gapflat[nkx][nky])).imag
										dosflat=dosflat+A2
										A3=-(1.0/pi)*(1.0/(z-gapminus[nkx][nky])).imag
										dosminus=dosminus+A3
						s=fermi(omega,0.0,temp)
						sum_entropy=sum_entropy-2.0*(s*log(s))*(dosplus/(Nk*Nk)+dosflat/(Nk*Nk)+dosminus/(Nk*Nk))
				sum_entropy_tot=sum_entropy*(12.0/totomega)
				entropy.write("%6.2f %6.2f\n" % (temp, sum_entropy_tot))
		entropy.close()

def fermi(energy , Ef, temp):
		fermi=1.0/(exp((energy-Ef)/temp)+1.0)
		
		return float(fermi.real)
