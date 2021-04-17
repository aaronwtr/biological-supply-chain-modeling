import numpy as np

# we do all units in SI units 

um_to_m = 10**(-6)
cm_to_m = 10**(-3)
h_to_s  = 3600

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
