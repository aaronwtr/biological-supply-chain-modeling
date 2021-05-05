import numpy as np
import matplotlib.pyplot as plt

# This code is exclusively for a 1D simulation




def petri_rod(b,l,p,Vc, Vv):
    
    
           
    
    B = np.log(b)
    L = np.log(l)
    P = np.log(p)
    [Drate,Dpr,G,Omega,h,Theta,Length_petri] = Vc  
    
    [N,dt,dur, b_min, l_min, p_min] = Vv

    
    
    
    # Left_BOperator = [((B[1]-B[0]+2)**2-4)/(4*dx**2)]
    # Mid_BOperator = ((B[2:N] - B[0:N-2] + 2)**2 + 8*(B[0:N-2] -  B[1:N-1]) -4)/(4*dx**2)
    # Right_BOperator = [((B[N-2]-B[N-1]+2)**2-4)/(4*dx**2)]
    
    # BOperator = np.concatenate((Left_BOperator,Mid_BOperator,Right_BOperator))
    
    # Left_LOperator = [((L[1]-L[0]+2)**2-4)/(4*dx**2)]
    # Mid_LOperator = ((L[2:N] - L[0:N-2] + 2)**2 + 8*(L[0:N-2] -  L[1:N-1]) -4)/(4*dx**2)
    # Right_LOperator = [((L[N-2]-L[N-1]+2)**2-4)/(4*dx**2)]
    
    # LOperator = np.concatenate((Left_LOperator,Mid_LOperator,Right_LOperator))
    
    # Left_POperator = [((P[1]-P[0]+2)**2-4)/(4*dx**2)]
    # Mid_POperator = ((P[2:N] - P[0:N-2] + 2)**2 + 8*(P[0:N-2] -  P[1:N-1]) -4)/(4*dx**2)
    # Right_POperator = [((P[N-2]-P[N-1]+2)**2-4)/(4*dx**2)]
    
    # POperator = np.concatenate((Left_POperator,Mid_POperator,Right_POperator))
    
    dx = Length_petri/N
    DX = 1/(4*dx**2)
    
    #print(dx)
    
    Steps = np.round(dur/(dt*samples))
    #print(Steps)
    Dtau = Drate * dt # um^2
    #print(Dtau)
    #Steps = 1
    i = 0 
    
    # DtauG = Dtau * G
    # DtauDX = DX * Dtau
    # DtauDXDpr = DtauDX * Dpr
    # Dtauh = Dtau*h
    # DtauTheta = Dtau*Theta
    minh = -h
    while i <  Steps: 
        
        # B1 = B + G - np.divide(np.exp(P),np.exp(B)+np.exp(L)+Kb) + np.concatenate(([((B[1]-B[0]+2)**2-4)],((B[2:N] - B[0:N-2] + 2)**2 + 8*(B[0:N-2] -  B[1:N-1]) -4),[((B[N-2]-B[N-1]+2)**2-4)]))/(4*dx**2)
        # L1 = L + np.exp(B-L)*np.divide(np.exp(P),np.exp(B)+np.exp(L)+Kb) - Omega + np.concatenate(([((L[1]-L[0]+2)**2-4)],((L[2:N] - L[0:N-2] + 2)**2 + 8*(L[0:N-2] -  L[1:N-1]) -4),[((L[N-2]-L[N-1]+2)**2-4)]))/(4*dx**2)
        # P1 = P - np.divide(T*h,T+Kb) + np.exp(L-P) * Theta + Dpr * np.concatenate(([((P[1]-P[0]+2)**2-4)],((P[2:N] - P[0:N-2] + 2)**2 + 8*(P[0:N-2] -  P[1:N-1]) -4),[((P[N-2]-P[N-1]+2)**2-4)]))/(4*dx**2)
        
        # B1 = B + G - F + np.concatenate(([((B[1]-B[0]+2)**2-4)],((B[2:N] - B[0:N-2] + 2)**2 + 8*(B[0:N-2] -  B[1:N-1]) -4),[((B[N-2]-B[N-1]+2)**2-4)]))/(4*dx**2)
        # L1 = L + np.exp(B-L)*F - Omega + np.concatenate(([((L[1]-L[0]+2)**2-4)],((L[2:N] - L[0:N-2] + 2)**2 + 8*(L[0:N-2] -  L[1:N-1]) -4),[((L[N-2]-L[N-1]+2)**2-4)]))/(4*dx**2)
        # P1 = P - np.divide((np.exp(B)+np.exp(L))*h,np.exp(B)+np.exp(L)+Kb) + np.exp(L-P) * Theta + Dpr * np.concatenate(([((P[1]-P[0]+2)**2-4)],((P[2:N] - P[0:N-2] + 2)**2 + 8*(P[0:N-2] -  P[1:N-1]) -4),[((P[N-2]-P[N-1]+2)**2-4)]))/(4*dx**2)
        
        LminP = L - P
        B1 = B + Dtau *(G - np.divide(1,np.exp(B-P)+np.exp(LminP)+Kb/np.exp(P)) + np.concatenate(([((B[1]-B[0]+2)**2-4)],((B[2:N] - B[0:N-2] + 2)**2 + 8*(B[0:N-2] -  B[1:N-1]) -4),[((B[N-2]-B[N-1]+2)**2-4)]))*DX )
        L1 = L + Dtau * (np.divide(1,np.exp(LminP)+np.exp(L-B+LminP)+Kb/np.exp(B-LminP)) - Omega + np.concatenate(([((L[1]-L[0]+2)**2-4)],((L[2:N] - L[0:N-2] + 2)**2 + 8*(L[0:N-2] -  L[1:N-1]) -4),[((L[N-2]-L[N-1]+2)**2-4)]))*DX )
        P = P + Dtau *(np.divide(minh,1+Kb/((np.exp(B)+np.exp(L)))) + np.exp(LminP) * Theta + Dpr * np.concatenate(([((P[1]-P[0]+2)**2-4)],((P[2:N] - P[0:N-2] + 2)**2 + 8*(P[0:N-2] -  P[1:N-1]) -4),[((P[N-2]-P[N-1]+2)**2-4)]))*DX )
        
        
        
        
        # LminP = L - P
        # B1 = B +  DtauG - np.divide(Dtau,np.exp(B-P)+np.exp(LminP)+Kb/np.exp(P)) + np.concatenate(([((B[1]-B[0]+2)**2-4)],((B[2:N] - B[0:N-2] + 2)**2 + 8*(B[0:N-2] -  B[1:N-1]) -4),[((B[N-2]-B[N-1]+2)**2-4)]))*DtauDX 
        # L1 = L + np.divide(1,np.exp(LminP)+np.exp(L-B+LminP)+Kb/np.exp(B-LminP)) - Omega + np.concatenate(([((L[1]-L[0]+2)**2-4)],((L[2:N] - L[0:N-2] + 2)**2 + 8*(L[0:N-2] -  L[1:N-1]) -4),[((L[N-2]-L[N-1]+2)**2-4)]))*DtauDX 
        # P = P - np.divide(Dtauh,1+Kb/((np.exp(B)+np.exp(L)))) + np.exp(LminP) * DtauTheta + np.concatenate(([((P[1]-P[0]+2)**2-4)],((P[2:N] - P[0:N-2] + 2)**2 + 8*(P[0:N-2] -  P[1:N-1]) -4),[((P[N-2]-P[N-1]+2)**2-4)]))*DtauDXDpr
        
        
        
        B = B1
        L = L1
        
        B = np.where(B < b_min, b_min, B)
        L = np.where(L < l_min, l_min, L)
        P = np.where(P < p_min, p_min, P)

        
        i += 1
    
    
    
    
    return np.exp(B) , np.exp(L) , np.exp(P)



# #--------------parameters---------------------- 
    
DBmax   = 5.25 * 10**5  # um^2/h
Kv      = .001          # um^-2
Kc      = .02           # um^-2
alphac  = 2             # 1
Dp      = 1             # um^2/h
Dn      = 4.5*10**6     # um^2/h
eta     = 8 * 10**4     # um^2/h
beta    = 80            # 1              
kl      = 2#00             # um^2/h
gmax    = 6             # 1/h
Kn      =  .1           # um^-2
Lambda  = .2            # 1 
Kb      =  .1           # um^-2
B0      = 10**7         # 1
P0      = 1.5 * 10**8   # 1
n0      = 1             # um^-2
Length_petri = 1000     # um   
B_0 = 1                 # um^-2
L_0 = 1                 # um^-2
P_0 = 1                 # um^-2
B_min = 10**(-10)       # um^-2 = 1 cm^-2
L_min = 10**(-10)       # um^-2 = 1 cm^-2
P_min = 10**(-10)       # um^-2 = 1 cm^-2


#--------------Condensed constants---------------------- 

Drate = DBmax * (n0/(n0 +Kv))**2                # um^2/h
Dpr = Dp  / Drate                               # 1
G = gmax * n0 / (n0 + Kn) /Drate                # um^-2
Omega = kl * n0 / Drate                         # um^-2
h = eta * Kb / Drate                            # um^-2
Theta  = beta * kl * n0 * eta * Kb / Drate**2   # um^-4

#vector of constants
Vc = [Drate,Dpr,G,Omega,h,Theta,Length_petri]


#------------------Simulation Variables---------------------

N = 10 #discretisaion in X (N)
dt = 10**(-8)  #time steps in h 
dur = .0030 # duration of total simulation in h
samples = 30 #how often one wants data to be returned

Bresults = np.zeros((N,samples+1))
Lresults = np.zeros((N,samples+1))
Presults = np.zeros((N,samples+1))

# parameters for bacteria, infected bacteria and phages in # um^-2
min_val = np.finfo(float).eps
B = min_val*np.ones((N)) ; L = min_val*np.ones((N)) ; P = min_val*np.ones((N))
B = np.ones((N)) ; L = np.ones((N))*10**(-6) ; P = np.ones((N))*10**(-2)
B[0] = 100
# convert BLP to blp

b = B/B_0 
b_min = np.log(B_min/B_0)
l = L/L_0
l_min = np.log(L_min/L_0)
p = P*h/P_0
p_min = np.log(P_min*h/P_0)



#vector of variables
Vv = [N,dt,dur, b_min, l_min, p_min]

Bresults[:,0] , Lresults[:,0]  , Presults[:,0]  = B , L , P


for i in range(samples):
    b,l,p = petri_rod(b,l,p,Vc, Vv)
    Bresults[:,i+1] , Lresults[:,i+1]  , Presults[:,i+1]  = b*B_0 , l*L_0 , p*P_0/h



fig, axs = plt.subplots(3)
fig.suptitle('results over time')
axs[0].imshow(np.log(Bresults), cmap='viridis')
axs[0].set_title('Bacteria')
axs[0].set(ylabel='Distance')
axs[1].imshow(np.log(Lresults), cmap='viridis')
axs[1].set_title('Infected bacteria')
axs[1].set(ylabel='Distance')
axs[2].imshow(np.log(Presults), cmap='viridis')
axs[2].set_title('Phages')
axs[2].set(xlabel='Time', ylabel='Distance')

#fig.colorbar()
plt.show()

print('the end results')
print(b*B_0)
print(l*L_0)
print(p*P_0/h)

#F = np.divide(p,b+l+Kb)









