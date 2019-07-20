# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 12:41:26 2018

@author: Hangjian Zhao
"""

###############################################################################
#This script is based on a solution proposed by Bower, J. W., Motz, L. H., & Durden, D. W. (1999)(Reference below). 
####Using scipy package, this graphic solution could be used to determine the critical condtion
###(i.e critical rise raio and critical pumping rate) of saltwater upconing in a leaky aquifer system.

##Script consists of two components: 
#Base graph (a series of solutions based on different sets of parameters) 
#Solution (Solution from "Parameters.txt" user input)


from scipy import special
import numpy as np
import matplotlib.pyplot as plt

#Define the constants
sigma=40           #Specific weight of freshwater/(Specific weight of saltwater-Specific weight of freshwater)
pi=np.pi           #Define the pi

#Define the container for storing the results
class Results():
    def __init__(self, y1, y2,delta_ratio):
        self.y1 = y1
        self.y2 = y2
        self.delta_ratio = delta_ratio
        
#Define the aquifer class container
class Aquifer():
    def __init__(self,h,m,l,b,d,r,Kz,Kr,k1,b1,T,n_all):
        self.h=h
        self.m=m
        self.l=l
        self.b=b
        self.d=d
        self.r=r
        self.Kz=Kz
        self.Kr=Kr
        self.k1=k1
        self.b1=b1
        self.T=T
        self.n_all=n_all
        
#Define the function (EqualPressure function and EqualPressureGradient function)
def Calculation(h,m,l,b,d,r,Kz,Kr,k1,b1,T,n_all):
    delta_ratio=np.arange(0.1,1.01,0.01)  
    y1=[]
    y2=[]
    for i in range(delta_ratio.size):
        fsn=np.zeros((n_all))
        fsn_deriv=np.zeros((n_all))
        
        for n in range(n_all+1)[1:]:
            fsn[n-1]=1/n*(np.sin(n*pi*l/b)-np.sin(n*pi*d/b))*np.cos(n*pi*(1-(delta_ratio[i])*(1-l/b)))*special.kn(0,n*pi*m)  #
            fsn_deriv[n-1]=(np.sin(n*pi*l/b)-np.sin(n*pi*d/b))*np.sin((n*pi)*(1-(delta_ratio[i])*(1-l/b)))*special.kn(0,n*pi*m)
        
        fscum=np.sum(fsn)
        fs=4*(pi*(l/b-d/b))**(-1.0)*fscum 
        suby1=2*pi*(1-l/b)*(delta_ratio[i])/(special.kn(0,h)+fs/2.0)
        y1.append(suby1)
      
        fsderivcum=np.sum(fsn_deriv)       
        fs_devz=(-4/b)*((l/b-d/b)**(-1.0))*fsderivcum
        suby2=-4*pi*((b*fs_devz)**(-1.0))
        y2.append(suby2)
        
    return Results(y1,y2,delta_ratio)
        
#Define the function for reading the input parameters for txt.file
def ReadParams():
    param = []
    with open('/Parameters.txt','r') as input_file:
        for line in input_file:
            input_line=line.split()
            param.append(input_line[1])
    #Reading input parameters
    l=float(param[0])                 #Distance from the top of the aquifer to the bottom of the well screen (m)
    b=float(param[1])                 #Thickness of the aquifer from the bottom of the confining unit to the initial condition of the saltwater-freshwater interface (m)
    d=float(param[2])                 #Distance from the top of the aquifer to the top of the screen (m)
    r=float(param[3])                 #Radial distance from the pumping well (rw) (m)
    Kz=float(param[4])                #Vertical hydraulic conductivity of the aquifer (m/day)
    Kr=float(param[5])                #Horizontal hydraulic conductivity of the aquifer (m/day)
    k1=float(param[6])                #Vertical hydraulic conductivity of the overlying confining unit (m/day)
    b1=float(param[7])                #Thickness of the overlying confining unit (m)
    T=float(param[8])                 #Transmissivity of the aquifer (m2/day)
    n_all=int(param[9])               #Summation factor (Define the accuracy of the solution, with a greater summation factor, accuracy improves and runtime increases) (dimensionless)
    
    #Calculating secondary parameters
    h=r*((k1/(b1*T))**(1/2))
    m=(Kz/Kr)**(1/2)*r/b
    aquifer=Aquifer(h,m,l,b,d,r,Kz,Kr,k1,b1,T,n_all)
    return aquifer

def axstyle(ax): 
    ax.set_yscale('log')
    ax.set_xlim(0.1,1.0)
    ax.set_xlabel('Interface rise ratio',fontsize=8.0)
    ax.set_ylabel('Nondimensional pumping rate 40Q/(Tb)',fontsize=8.0)
    ax.set_yscale('log')
    ax.tick_params(labelsize=8,axis='both',colors='grey')


################################Execution######################################
#%Construct the base graph based on a sereis of parameter sets (Solution is equivalent to Fig 4 in the paper referenced)
mlist=[1e-2,1e-3,1e-4,1e-5]  
Q1all=[]
Q2all=[]
for x in mlist:
    results=Calculation(1e-5,x,6,10,0,0.3,0.2,2,0.0001,2,20,1000)
    Q1all.append(results.y1)
    Q2all.append(results.y2)
    x=results.delta_ratio
    

############Execution the function based on input parameters###################
aquifer=ReadParams()
results=Calculation(aquifer.h,aquifer.m,aquifer.l,aquifer.b,aquifer.d,aquifer.r,aquifer.Kz,aquifer.Kr,aquifer.k1,aquifer.b1,aquifer.T,aquifer.n_all)
calc1=np.array(results.y1)
calc2=np.array(results.y2)
sloc=np.argmin((abs(calc1-calc2)),axis=0)       #Find the index where the smallest difference exists
Critical_rise_ratio=x[sloc]                     #Find the corresponding critical rise ratio
Critical_rise=Critical_rise_ratio*(aquifer.b-aquifer.l)
Q=(calc1[sloc]*aquifer.T*aquifer.b)/sigma       #Backcalculate the critical pumping rate Q

#%%Graphing

Q1col=['lawngreen','greenyellow','forestgreen','darkolivegreen']
Q2col=['skyblue','darkturquoise','deepskyblue','darkcyan']

fig,(ax1,ax2)=plt.subplots(2,figsize=(6.0,9.0))
for i in range(len(Q1col)):
    ax1.plot(x,Q1all[i],'--',c=Q1col[i],alpha=0.8)
    ax1.plot(x,Q2all[i],'--',c=Q2col[i],alpha=0.8)

axstyle(ax1)

#Plot the solution from input parameters
ax2.plot(x,calc1,'-',c='red',alpha=1.0)
ax2.plot(x,calc2,'-',c='blue',alpha=1.0)
ax2.plot(Critical_rise_ratio,calc1[sloc],'*',c='k',ms=10.0)             #Find the corresponding critical rise ratio)
ax2.text(0.12,0.0015,'Critical rise: %.1f m\nCritical pumping rate: %.1f $m^3$/day'%(Critical_rise,Q))

axstyle(ax2)

#Reference:
#Bower, J. W., Motz, L. H., & Durden, D. W. (1999). 
#Analytical solution for determining the critical condition of saltwater upconing in a leaky artesian aquifer. Journal of Hydrology, 221(1-2), 43-54.





















    
    
    

