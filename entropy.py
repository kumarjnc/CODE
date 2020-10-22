# python-env/2.7.13 should be loaded 

import numpy as np
#-----------------------------------
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

def fit(input_dmft,temp_dmft,mu_min,mu_max,temp_min,temp_max,Nfit_mu,Nfit_temp):
	from scipy import interpolate
	print (temp_dmft)
	N=len(temp_dmft)
	M=int(len(input_dmft[:,0])/N)
	print (M,N)

	
	mu_dmft=input_dmft[0:M,0]
	densA=np.absolute(input_dmft[:,1]);
	densB=np.absolute(input_dmft[:,2]);

	dbleA=np.absolute(input_dmft[:,4]);
	dbleB=np.absolute(input_dmft[:,5]);

	densA=densA.reshape(N,M);
	densB=densB.reshape(N,M);

	dbleA=dbleA.reshape(N,M);
	dbleB=dbleB.reshape(N,M);
#------------------Fitting-----------------------------------------------------
	f = interpolate.interp2d(mu_dmft, temp_dmft, densA, kind='linear')
	g = interpolate.interp2d(mu_dmft, temp_dmft, densB, kind='linear')


	h = interpolate.interp2d(mu_dmft, temp_dmft, dbleA, kind='linear')
	k = interpolate.interp2d(mu_dmft, temp_dmft, dbleB, kind='linear')

#-----------------------Updated Grid-------------------------------------------


	mu_new=np.linspace(mu_min,mu_max,num=Nfit_mu,endpoint=True);
	temp_new=np.linspace(temp_min, temp_max,num=Nfit_temp,endpoint=True);

#-----------------------------------New density -----------------------
	densAnew=f(mu_new,temp_new); # New Fitted data for A site
	densBnew=g(mu_new,temp_new); # New Fitted data for B site
#------------------------------------New double------------------------


	dbleAnew=h(mu_new,temp_new); # New Fitted data for A site; double occupancy
	dbleBnew=k(mu_new,temp_new); # New Fitted data for B site; double occpancy
#-----------------------------------------------------------------------------------


	f=open('Fitted_data.dat','w')

	for i in range(Nfit_temp):
		for j in range(Nfit_mu):
			f.write("%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n" %(mu_new[j], temp_new[i],densAnew[i,j],densBnew[i,j],dbleAnew[i,j],dbleBnew[i,j]))
	f.close

	return densAnew, densBnew, dbleAnew, dbleBnew
#---------------------------------------------------------------------------------------------------------


def difference(Nfit_temp,Nfit_mu,temp_min,temp_max, densA,densB, doublA, doubleB):

	temp_h=(temp_max-temp_min)/Nfit_temp;     # temperature grid

	dfdtA=np.zeros((Nfit_temp,Nfit_mu))       # df(temp,mu)_A/dt 
	dfdtB=np.zeros((Nfit_temp,Nfit_mu))       

	print (temp_h)
	for i in range(1,Nfit_temp-1):
		for j in range(Nfit_mu):
			dfdtA[i,j]=(densA[i+1,j]-densB[i-1,j])/(2*temp_h)
			dfdtB[i,j]=(densB[i+1,j]-densB[i-1,j])/(2*temp_h)

	for j in range(Nfit_mu):
        	dfdtA[0,j]=(densA[1,j]-densB[0,j])/temp_h
        	dfdtA[-1,j]=(densA[-1,j]-densA[-2,j])/temp_h
        	dfdtB[0,j]=(densB[1,j]-densB[0,j])/temp_h
        	dfdtB[-1,j]=(densB[-1,j]-densB[-2,j])/temp_h


	#g=open('Derivative.dat','w')
	#for j in range(Nfit_mu):
        #	g.write("%6.4f %6.4f\n" %(dfdtA[i,j],dfdtB[i,j]))
	#g.close
	
	return dfdtA, dfdtB
#------------------------------------------------------------------------------------------------------------------------------------
def ent(temp_index,mu_min,mu_max,Nfit_mu,Nfit_temp,densA,densB,dfdtA,dfdtB):

	mu_h=(mu_min-(-mu_max))/Nfit_mu;
	print (mu_h)
	print ('This is temperature')
	i=temp_index

	j=1
	while True:
        	j=j+1
        	if (densA[i,j] >= 0.008):
                	llim=j
                	break

	print (llim,zzAnew[i,llim])

	sumA=np.zeros(Nfit_mu)
	sumB=np.zeros(Nfit_mu)


	s=open("Entropy.dat",'w')
	for mu in range(llim+2,Nfit_mu):
        	sumdA=(densA[i,llim]+densA[i,mu])*0.5 #density sum
        	sumdB=(densB[i,llim]+densB[i,mu])*0.5
        	sA=(dfdtA[i,llim]+dfdtA[i,mu])*0.5      #Entropy
        	sB=(dfdtB[i,llim]+dfdtB[i,mu])*0.5
	for j in range(llim+1,mu-1):
		sA=sA+dfdtA[i,j]
		sB=sB+dfdtB[i,j]
		sumdA=sumdA+densA[i,j]
		sumdB=sumdB+densB[i,j]
	sumA[mu]=hmu*sA
	sumB[mu]=hmu*sB
	sumdtotA=hmu*sumdA
	sumdtotB=hmu*sumdB

	print (sumtotA, sumtotB)

	#for j in range(Nfit_mu):
        #	s.write("%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n"  %(zzAnew[i,j],zzBnew[i,j],dAnew[i,j],dBnew[i,j],sumA[j]+2.0*sumB[j]))
	#s.close()

	return sumA, sumB
