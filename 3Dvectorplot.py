import numpy as np
import matplotlib.pyplot as plt
from numba import njit
import time
start = time.time()


@njit()
def threeDvecplotter(b, l, p, Vec, Drate, Dpr, G, Omega, h, logTheta):
    
    
    hKb = h*Kb
    
    
    
    i = 0
    
    for b_now in b:
        i+=1
        j = 0
        for l_now in l:
            k = 0
            j +=1
            one_over_ExpB_plus_ExpL_plus_Kb = 1/(np.exp(b_now)+np.exp(l_now)+Kb)
    
            for p_now in p:
                
                Vec[i,j,k,0] =  G - np.exp(p_now)*one_over_ExpB_plus_ExpL_plus_Kb
                Vec[i,j,k,1] = np.exp(p_now+l_now-b_now)*one_over_ExpB_plus_ExpL_plus_Kb - Omega
                Vec[i,j,k,2] = -h+hKb*one_over_ExpB_plus_ExpL_plus_Kb + np.exp(l_now - p_now+logTheta) 
                
                k+=1
    return Vec




# here we make a vectorplot of the log values of the B L & P:

    
    # #--------------parameters---------------------- 
    
DBmax   = 5.25 * 10**5  # um^2/h
Kv      = .001          # um^-2
#Kc      = .02           # um^-2
#alphac  = 2             # 1
Dp      = 1             # um^2/h
#Dn      = 4.5*10**6     # um^2/h
eta     = 8 * 10**4     # um^2/h
beta    = 80            # 1              
kl      = 2#00             # um^2/h
gmax    = 6             # 1/h
Kn      =  .1           # um^-2
#Lambda  = .2            # 1 
Kb      =  .1           # um^-2
#B0      = 10**7         # 1
#P0      = 1.5 * 10**8   # 1
n0      = 1             # um^-2
Length_petri = 1000     # um   
B_0 = 1                 # um^-2
L_0 = 1                 # um^-2
P_0 = 1                 # um^-2



#--------------Condensed constants---------------------- 

Drate = DBmax * (n0/(n0 +Kv))**2                # um^2/h
Dpr = Dp  / Drate                               # 1
G = gmax * n0 / (n0 + Kn) /Drate                # um^-2
Omega = kl * n0 / Drate                         # um^-2
h = eta * Kb / Drate                            # um^-2
Theta  = beta * kl * n0 * eta * Kb / Drate**2   # um^-4
logTheta = np.log(Theta)

# ------------ user set parameters --------------------

b_mean = 0; b_dist = 2; b_num = 21
l_mean = 0; l_dist = 2; l_num = 21
p_mean = 0; p_dist = 2; p_num = 21

# ------------ computations --------------------------

b = np.linspace(b_mean-b_mean, b_mean+b_dist, num=b_num,dtype = float)
l = np.linspace(l_mean-l_dist, l_mean+l_dist, num=l_num,dtype = float)
p = np.linspace(p_mean-p_dist, p_mean+p_dist, num=p_num,dtype = float) 

Vec = np.zeros((len(b),len(l),len(p),3))

    
Vec = threeDvecplotter(b,l,p,Vec, Drate, Dpr, G, Omega, h, logTheta)
print("The simulation took " + str(time.time() - start) + " seconds")