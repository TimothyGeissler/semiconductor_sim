import sys
import math

Q = 1.602177*math.pow(10, -19)
ni = 1.07*math.pow(10, 10) #Intrinsic carrier concentration
kB = 1.380649*math.pow(10, -23)#0.00008617333262
T = 300
epsilon_r = 11.68 #rel permitivity of Si
epsilon_o = 8.854*math.pow(10, -14) #permitivity of vacuum

## Room temp, non-compensated n-type Si with Nd dopants ##

# Majority carrier concentration
def nN(Nd):
    return (0.5)*(Nd + math.sqrt((Nd)*(Nd) + 4*(ni)*(ni)))

# Minority carrier concentration
def pN(Na):
    return (2*(ni)*(ni))/(Nd + math.sqrt((Nd)*(Nd) + 4*(ni)*(ni)))

# Majority electron low-field bulk mobility (cm^2/Vs) 92 + 1268/(1 + (Nd(1.3×10^17))^0.91)
def mu_nN(Na, Nd):
    return 92 + 1268/(1 + math.pow(((Na + Nd)/(1.3*math.pow(10, 17))), 0.91))

#  Minority holes low-field bulk mobility
def mu_pN(Na, Nd):
    return 130 + 370/(1 + math.pow(((Na + Nd)/(8*math.pow(10, 17))), 1.25));


## Room temp,  siliconnon-compensated p-type Si with Na dopants ##

# Majority carrier concentration
def pP(Na):
    return (0.5)*(Na + math.sqrt((Na)*(Na) + 4*(ni)*(ni)))

# Minority carrier concentration
def nP(Na):
    return (2*(ni)*(ni))/(Na + math.sqrt((Na)*(Na) + 4*(ni)*(ni)))

# Majority hole low-field bulk mobility (cm^2/Vs)
def mu_pP(Na, Nd):
    return 49.7 + 418.3/(1 + math.pow(((Na + Nd)/(1.6*math.pow(10, 17))), 0.7))

#  Minority electron low-field bulk mobility
def mu_nP(Na, Nd):
    return 232 + 1180/(1 + math.pow(((Na + Nd)/(8*math.pow(10, 16))), 0.9))


## Compensated Silicon Regions ##

# Hole concentration in p-type
def compensated_pP(Na, Nd):
    return (0.5)*((Na - Nd) + math.sqrt(math.pow((Na - Nd), 2) + 4*ni*ni))

# Electron concentration in p-type
def compensated_nP(Na, Nd):
    return (ni*ni)/compensated_pP(Na, Nd)

# Electron concentration in n-type
def compensated_nN(Na, Nd):
    return (0.5)*((Nd - Na) + math.sqrt(math.pow((Nd - Na), 2) + 4*ni*ni))

# Hole concentration in n-type
def compensated_pN(Na, Nd):
    return (ni*ni)/compensated_nN(Na, Nd)

# Resistivity of an n-type Si (Ohm.cm)
def rho_N(Nd, n, p):
    return 1/(Q * (mu_nN(Na, Nd) * n + mu_pN(Na, Nd) * p))

# Resistivity of a p-type Si (Ohm.cm)
def rho_P(Na, n, p):
    return 1/(Q * (mu_nP(Na, Nd) * n + mu_pP(Na, Nd) * p))

# Electron drift velocity
def v_ndrift(E, mu_n):
    return -mu_n*E

# Hole drift velocity
def v_pdrift(E, mu_p):
    return mu_p*E

# Static electron drift current density (A/cm^2)
def j_nxdrift(n, E, mu_n):
    return -Q*n*v_ndrift(E, mu_n)

# Static hole drift current density
def j_pxdrift(p, E, mu_p):
    return Q*p*v_pdrift(E, mu_p)

def diff_coeff(T, mu):
    return (kB * T * mu)/Q

def tau_rec_p(Na):
    try:
        return 1/((3.45*10**-12) * Na + (9.5*10**-32)*(Na*Na))
    except:
        return 0
    

def tau_rec_n(Nd):
    try:
        return 1/((7.8*10**-13) * Nd + (1.8*10**-31)*(Nd*Nd))
    except:
        return 0
    
## PN JUNCTION DIODES ##

# Builtin voltage - Th Eq (Volts)
def v_bi(Na, Nd):
    return ((kB * T)/Q) * math.log(Na * Nd / (ni**2))

# Builtin voltage - p side
def V_biP(Na, Nd):
    return ((Q * Na) / (2 *  epsilon_o * epsilon_r)) * W_dP(Na, Nd)**2

# Builtin voltage - n side
def V_biN(Na, Nd):
    return ((Q * Nd) / (2 *  epsilon_o * epsilon_r)) * W_dN(Na, Nd)*W_dN(Na, Nd)

# P-side depletion width (cm)
def W_dP(Na, Nd):
    return math.sqrt((2*epsilon_o*epsilon_r / Q) * (Nd/(Na*(Na + Nd))) * v_bi(Na, Nd)) 

# N-side depletion width (cm)
def W_dN(Na, Nd):
    return math.sqrt((2*epsilon_o*epsilon_r / Q) * (Na/(Nd*(Na + Nd))) * v_bi(Na, Nd)) 

# Total depletiuon width
def W_d(Na, Nd):
    return W_dN(Na, Nd) + W_dP(Na, Nd)

# Max electric field
def E_max(Na, Nd):
    #return ((-Q * Na)/(epsilon_o * epsilon_r)) * W_dP(Na, Nd)
    return -math.sqrt(((2 * Q)/(epsilon_o * epsilon_r)) * ((Na * Nd)/(Na + Nd)) * (v_bi(Na, Nd)))

# Electron distribution p-side depletion region
def n_p_dep(Na, Nd, x):
    return ((ni**2)/Na) * math.exp((((Q**2)*Na)/(2*epsilon_o*epsilon_r*kB*T))*(W_dP(Na, Nd) + x)**2)

# Hole distribution p-side depletion region
def p_p_dep(Na, Nd, x):
    return Na*math.exp(-(((Q**2)*Na)/(2*epsilon_o*epsilon_r*kB*T))*(W_dP(Na, Nd) + x)**2)

# Electron distribution n-side depletion region
def n_n_dep(Na, Nd, x):
    return Nd*math.exp(-(((Q**2)*Nd)/(2*epsilon_o*epsilon_r*kB*T))*(W_dN(Na, Nd) - x)**2)

# Hole distribution n-side depletion region
def p_n_dep(Na, Nd, x):
    return ((ni**2)/Nd) * math.exp((((Q**2)*Nd)/(2*epsilon_o*epsilon_r*kB*T))*(W_dN(Na, Nd) - x)**2)

## PRINT STATEMENTS ##
def comp_concentration(Na, Nd):
    if (Nd > Na):
        # Compensated n-type  
        print("\n## For compensated n-type (cm^-3) ##\nMajority nN=" + "{:E}".format(compensated_nN(Na, Nd)) + "\tMinority pN=" + "{:E}".format(compensated_pN(Na, Nd)))
        mu_ntype(Nd)
        print("\n## Resistivity (n-type) (Ohm.cm) ##\nrho_n=" + str(rho_N(Nd, compensated_nN(Na, Nd), compensated_pN(Na, Nd))))
        print("\n## Drift velocity (n-type)\nElectron drift v_n,drift=" + "{:E}".format(v_ndrift(E, mu_nN(Na, Nd))) + "\tHole drift v_p,drift=" + "{:E}".format(v_pdrift(E, mu_pN(Na, Nd))))
        print("\n## Diffusion Coefficient (n-type) (cm^2/s) ##\nD_p,N=" + str(diff_coeff(T, mu_pN(Na, Nd))) + "\tD_n,N=" + str(diff_coeff(T, mu_nN(Na, Nd))))
        print("\n## Recombination Lifetime (Sec) ##\ntau_rec,P= " + "{:E}".format(tau_rec_p(Na)) + "\ttau_rec,N=" + "{:E}".format(tau_rec_n(Nd)))
    else:
        # Compensated p-type
        print("\n## For compensated p-type (cm^-3) ##\nMajority pP=" + "{:E}".format(compensated_pP(Na, Nd)) + "\tMinority nP=" + "{:E}".format(compensated_nP(Na, Nd)))
        mu_ptype(Na)
        print("\n## Resistivity (p-type) (Ohm.cm) ##\nrho_p=" + str(rho_P(Na, compensated_nP(Na, Nd), compensated_pP(Na, Nd))))
        print("\n## Drift velocity (p-type)\nElectron drift v_n,drift=" + "{:E}".format(v_ndrift(E, mu_nP(Na, Nd))) + "\tHole drift v_p,drift=" + "{:E}".format(v_pdrift(E, mu_pP(Na, Nd))))
        print("\n## Diffusion Coefficient (p-type) (cm^2/s) ##\nD_n=" + str(diff_coeff(T, mu_nP(Na, Nd))) + str(diff_coeff(T, mu_pP(Na, Nd))))
        print("\n## Recombination Lifetime (sec) ##\ntau_rec,P= " + "{:E}".format(tau_rec_p(Na)) + "\ttau_rec,N=" + "{:E}".format(tau_rec_n(Nd)))

def print_pn_junc(Na, Nd):
    # Print PN Junction data
    print("\n## PN JUNCTION ##\n")
    print("\n## Builtin voltage (V) ##\nV_bi,P=" + str(V_biP(Na, Nd)) + "\tV_bi,N=" + str(V_biN(Na, Nd)) + "\tV_bi=" + str(v_bi(Na, Nd)))
    print("\n## Depletion Region Width (cm) ##\nW_d,p=" + str(W_dP(Na, Nd)) + "\tW_d,n=" + str(W_dN(Na, Nd)) + "\tW_d=" + str(W_d(Na, Nd)))
    print("\n## Max Elextric field (x=0) (V/cm) ##\nE_max=" + str(E_max(Na, Nd)))
    print("\n## Carrier concentration (Depletion Region) x=" + str(x) + " ##\nn_p(x)=" + str(n_p_dep(Na, Nd, x)) + "\tp_p(x)=" + str(p_p_dep(Na, Nd, x))
          + "\tn_n(x)=" + str(n_n_dep(Na, Nd, x)) + "\tp_n(x)=" + str(p_n_dep(Na, Nd, x)))

def mu_ptype(Na):
    print("\n## Low-field bulk mobility (p-type) (cm^2/V.s) ## \nMajority mu_p,P=" + str(mu_pP(Na, Nd)) + "\tMinority mu_n,P=" + str(mu_nP(Na, Nd)))

def mu_ntype(Nd):
    print("\n## Low-field bulk mobility (n-type) (cm^2/V.s) ## \nMajority mu_n,N=" + str(mu_nN(Na, Nd)) + "\tMinority mu_p,N=" + str(mu_pN(Na, Nd)))


## USAGE ##

Na = 0#3.2*math.pow(10, 16)
Nd = 2*math.pow(10, 16)

# PN Junction variables
x = 0

E = 2.5*math.pow(10, 3) #for drift velocity


comp_concentration(Na, Nd)
print_pn_junc(Na, Nd)
#print(j_pxdrift(pN(2*math.pow(10,16)), 2.5*math.pow(10, 3), mu_pN(Nd = 2*math.pow(10, 16))))
#print(rho_N(1*math.pow(10, 18), math.pow(10, 18), math.pow((1.07*math.pow(10,10)), 2)/ math.pow(10, 18)))