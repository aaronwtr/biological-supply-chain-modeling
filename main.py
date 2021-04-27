import numpy as np


def petri_rod(b,l,p,Vc, Vv):
    
    
    DBmax   = 5.25 * 10**5 * um_to_m**2 /h_to_s   # um^2/h to m^2/s
    Kv      = .001 * (10**(-6))**(-2)             # um^-2 to m^-2
    Kc      = .02 * (10**(-6))**(-2)              # um^-2 to m^-2
    alphac  = 2                                   #
    Dp      = 1 * ((10**(-6)))**2 /3600           # um^2/h to m^2/s
    Dn      = 4.5*10**6 * um_to_m**2 /h_to_s      #
    eta     = 8 * 10**4 * um_to_m**2 /h_to_s      # um^2/h to m^2/s
    beta    = 80                                  #
    kl      = 2 * um_to_m**2 /h_to_s              # um^2/h to m^2/s
    gmax    = 6 / h_to_s                          # 1/h to 1/s
    Kn      =  .1 * um_to_m**(-2)                 # um^-2 to m^-2
    Lambda  = .2                                  # 
    Kb      =  .1 * um_to_m**(-2)                 # um^-2 to m^-2
    B0      = 10**7                               #
    P0      = 1.5 * 10**8                         #
    n0      = 1 * um_to_m**(-2)                   # um^-2 to m^-2
    r       = .25 * cm_to_m                       # cm to m
    R       = 5.5 * cm_to_m                       # cm to m
    
    B = np.log(b); L = np.log(l); P = np.log(p); 
    T = np.exp(B)+np.exp(L)
    F = np.divide(np.exp(P),T+Kb)
    
    # Bydx = B[2:N] ; By0 = B[1:N-1]; By_dx = B[0:N-2]
    
    dx = 2*R/(N+1)
    
    Left_BOperator = ((B[1]-B[0]+2)**2-4)/(4*dx**2)
    Mid_BOperator = ((B[2:N] - B[0:N-2] + 2)**2 + 8*(B[0:N-2] -  B[1:N-1]) -4)/(4*dx**2)
    Right_BOperator = ((B[N-2]-B[N-1]+2)**2-4)/(4*dx**2)
    
    BOperator = np.concatenate((Left_BOperator,Mid_BOperator,Right_BOperator))
    
    Left_LOperator = ((L[1]-L[0]+2)**2-4)/(4*dx**2)
    Mid_LOperator = ((L[2:N] - L[0:N-2] + 2)**2 + 8*(L[0:N-2] -  L[1:N-1]) -4)/(4*dx**2)
    Right_LOperator = ((L[N-2]-L[N-1]+2)**2-4)/(4*dx**2)
    
    LOperator = np.concatenate((Left_LOperator,Mid_LOperator,Right_LOperator))
    
    Left_POperator = ((P[1]-P[0]+2)**2-4)/(4*dx**2)
    Mid_POperator = ((P[2:N] - P[0:N-2] + 2)**2 + 8*(P[0:N-2] -  P[1:N-1]) -4)/(4*dx**2)
    Right_POperator = ((P[N-2]-P[N-1]+2)**2-4)/(4*dx**2)
    
    POperator = np.concatenate((Left_POperator,Mid_POperator,Right_POperator))
    
    Drate = DBmax * (n0/(n0 +Kv))**2
    G = gmax * n0 / (n0 + Kn) /Drate
    Omega = kl * n0 / Drate
    Dpr = Dp  / Drate
    h = eta * Kb / Drate
    Theta  = beta * kl * n0 * eta * Kb / Drate**2
    
    
    B += G - F + BOperator
    L += np.exp(B-L)*F - Omega + LOperator
    P += np.divide(T*h,T+Kb) + np.exp(L-P) * Theta + Dpr * POperator
    
    
    
    
    
    
    return np.exp(B) , np.exp(L) , np.exp(P)

# we do all units in SI units 

um_to_m = 10**(-6)
cm_to_m = 10**(-3)
h_to_s  = 3600


#--------------parameters---------------------- 

DBmax   = 5.25 * 10**5 * um_to_m**2 /h_to_s   # um^2/h to m^2/s
Kv      = .001 * (10**(-6))**(-2)             # um^-2 to m^-2
Kc      = .02 * (10**(-6))**(-2)              # um^-2 to m^-2
alphac  = 2                                   #
Dp      = 1 * ((10**(-6)))**2 /3600           # um^2/h to m^2/s
Dn      = 4.5*10**6 * um_to_m**2 /h_to_s      #
eta     = 8 * 10**4 * um_to_m**2 /h_to_s      # um^2/h to m^2/s
beta    = 80                                  #
kl      = 2 * um_to_m**2 /h_to_s              # um^2/h to m^2/s
gmax    = 6 / h_to_s                          # 1/h to 1/s
Kn      =  .1 * um_to_m**(-2)                 # um^-2 to m^-2
Lambda  = .2                                  # 
Kb      =  .1 * um_to_m**(-2)                 # um^-2 to m^-2
B0      = 10**7                               #
P0      = 1.5 * 10**8                         #
n0      = 1 * um_to_m**(-2)                   # um^-2 to m^-2
r       = .25 * cm_to_m                       # cm to m
R       = 5.5 * cm_to_m                       # cm to m
#vector of constants
Vc = [DBmax,Kv,Kc,alphac,Dp,Dn,eta,beta,kl,gmax,Kn,Lambda,Kb,B0,P0,n0,r,R]
#--------------Condensed constants---------------------- 

# Drate = DBmax * (n0/(n0 +Kv))**2
# G = gmax * n0 / (n0 + Kn) /Drate
# Omega = kl * n0 / Drate
# Dpr = Dp  / Drate
# h = eta * Kb / Drate
# Theta  = beta * kl * n0 * eta * Kb / Drate**2



#------------------Simulation Variables---------------------

N = 10 ; M = 1 #discretisaion in X (N) and Y (M)
dt = 10**(-3)  #time step
dur = 10 # duration
samples = 10 #how often one wants data to be returned
Steps = np.round(dur/dt)
# parameters for bacteria, infected bacteria and phages
min_val = np.finfo(float).eps
b = min_val*np.ones((N)) ; l = min_val*np.ones((N)) ; p = min_val*np.ones((N))

#vector of variables
Vv = [N,M,dt,dur]

b,l,p = petri_rod(b,l,p,Vc, Vv)

print(b)

#F = np.divide(p,b+l+Kb)









