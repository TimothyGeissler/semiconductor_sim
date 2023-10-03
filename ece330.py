import sys
import math

Q = 1.602177*math.pow(10, -19)
ni = 1.07*math.pow(10, 10) #Intrinsic carrier concentration
kB = 0.00008617333262
T = 300
epsilon_r = 11.68 #rel permitivity of Si
epsilon_o = 8.854*math.pow(10, -12) #permitivity of vacuum

## Room temp, non-compensated n-type Si with Nd dopants ##

# Majority carrier concentration
def nN(Nd):
    return (0.5)*(Nd + math.sqrt((Nd)*(Nd) + 4*(ni)*(ni)))

# Minority carrier concentration
def pN(Na):
    return (2*(ni)*(ni))/(Nd + math.sqrt((Nd)*(Nd) + 4*(ni)*(ni)))

# Majority electron low-field bulk mobility (cm^2/Vs) 92 + 1268/(1 + (Nd(1.3×10^17))^0.91)
def mu_nN(Nd):
    return 92 + 1268/(1 + math.pow((Nd/(1.3*math.pow(10, 17))), 0.91))

#  Minority holes low-field bulk mobility
def mu_pN(Nd):
    return 130 + 370/(1 + math.pow((Nd/(8*math.pow(10, 17))), 1.25));


## Room temp, non-compensated p-type Si with Na dopants ##

# Majority carrier concentration
def pP(Na):
    return (0.5)*(Na + math.sqrt((Na)*(Na) + 4*(ni)*(ni)))

# Minority carrier concentration
def nP(Na):
    return (2*(ni)*(ni))/(Na + math.sqrt((Na)*(Na) + 4*(ni)*(ni)))

# Majority hole low-field bulk mobility (cm^2/Vs)
def mu_pP(Na):
    return 49.7 + 418.3/(1 + math.pow((Na/(1.6*math.pow(10, 17))), 0.7))

#  Minority electron low-field bulk mobility
def mu_nP(Na):
    return 232 + 1180/(1 + math.pow((Na/(8*math.pow(10, 16))), 0.9))


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
    return 1/(Q * (mu_nN(Nd) * n + mu_pN(Nd) * p))

# Resistivity of a p-type Si (Ohm.cm)
def rho_P(Na, n, p):
    return 1/(Q * (mu_nP(Na) * n + mu_pP(Na) * p))

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
    return ((kB * T)) * math.log(Na * Nd / (ni**2))

# Builtin voltage - p side
def V_biP(Na, Nd):
    return ((Q * Na) / (2 *  epsilon_o * epsilon_r)) * W_dP(Na, Nd)*W_dP(Na, Nd)

# Builtin voltage - n side
def V_biN(Na, Nd):
    return ((Q * Nd) / (2 *  epsilon_o * epsilon_r)) * W_dN(Na, Nd)*W_dN(Na, Nd)

# P-side depletion width (um)
def W_dP(Na, Nd):
    return math.sqrt((2*epsilon_o*epsilon_r / Q) * (Nd/(Na*(Na + Nd))) * v_bi(Na, Nd)) * 10**4

# N-side depletion width (um)
def W_dN(Na, Nd):
    return math.sqrt((2*epsilon_o*epsilon_r / Q) * (Na/(Nd*(Na + Nd))) * v_bi(Na, Nd)) * 10**4

# Total depletiuon width
def W_d(Na, Nd):
    return W_dN(Na, Nd) + W_dP(Na, Nd)

## PRINT STATEMENTS ##
def comp_concentration(Na, Nd):
    if (Nd > Na):
        # Compensated n-type  
        print("\n## For compensated n-type (cm^-3) ##\nMajority nN=" + "{:E}".format(compensated_nN(Na, Nd)) + "\tMinority pN=" + "{:E}".format(compensated_pN(Na, Nd)))
        mu_ntype(Nd)
        print("\n## Resistivity (n-type) (Ohm.cm) ##\nrho_n=" + str(rho_N(Nd, compensated_nN(Na, Nd), compensated_pN(Na, Nd))))
        print("\n## Drift velocity (n-type)\nElectron drift v_n,drift=" + "{:E}".format(v_ndrift(E, mu_nN(Nd))) + "\tHole drift v_p,drift=" + "{:E}".format(v_pdrift(E, mu_pN(Nd))))
        print("\n## Hole Diffusion Coefficient (n-type) ##\nD_n=" + str(diff_coeff(T, mu_pN(Nd))))
        print("\n## Recombination Lifetime (Sec) ##\ntau_rec,P= " + "{:E}".format(tau_rec_p(Na)) + "\ttau_rec,N=" + "{:E}".format(tau_rec_n(Nd)))
    else:
        # Compensated p-type
        print("\n## For compensated p-type (cm^-3) ##\nMajority pP=" + "{:E}".format(compensated_pP(Na, Nd)) + "\tMinority nP=" + "{:E}".format(compensated_nP(Na, Nd)))
        mu_ptype(Na)
        print("\n## Resistivity (p-type) (Ohm.cm) ##\nrho_p=" + str(rho_P(Na, compensated_nP(Na, Nd), compensated_pP(Na, Nd))))
        print("\n## Drift velocity (p-type)\nElectron drift v_n,drift=" + "{:E}".format(v_ndrift(E, mu_nP(Na))) + "\tHole drift v_p,drift=" + "{:E}".format(v_pdrift(E, mu_pP(Na))))
        print("\n## Electron Diffusion Coefficient (p-type) ##\nD_n=" + str(diff_coeff(T, mu_nP(Na))))
        print("\n## Recombination Lifetime (sec) ##\ntau_rec,P= " + "{:E}".format(tau_rec_p(Na)) + "\ttau_rec,N=" + "{:E}".format(tau_rec_n(Nd)))

def print_pn_junc(Na, Nd):
    # Print PN Junction data
    print("\n## PN JUNCTION ##\n")
    print("\n## Builtin voltage (V) ##\nV_bi,P=" + str(V_biP(Na, Nd)) + "\tV_bi,N=" + str(V_biN(Na, Nd)) + "\tV_bi=" + str(v_bi(Na, Nd)))
    print("\n## Depletion Region Width (um) ##\nW_d,p=" + str(W_dP(Na, Nd)) + "\tW_d,n=" + str(W_dN(Na, Nd)) + "\tW_d=" + str(W_d(Na, Nd)))

def mu_ptype(Na):
    print("\n## Low-field bulk mobility (p-type) (cm^2/V.s) ## \nMajority mu_p,P=" + str(mu_pP(Na)) + "\tMinority mu_n,P=" + str(mu_nP(Na)))

def mu_ntype(Nd):
    print("\n## Low-field bulk mobility (n-type) (cm^2/V.s) ## \nMajority mu_n,N=" + str(mu_nN(Nd)) + "\tMinority mu_p,N=" + str(mu_pN(Nd)))


## USAGE ##

Na = 7.5*math.pow(10, 16)
Nd = 4.2*math.pow(10, 17)

E = 2.5*math.pow(10, 3) #for drift velocity


comp_concentration(Na, Nd)
print_pn_junc(Na, Nd)
#print(j_pxdrift(pN(2*math.pow(10,16)), 2.5*math.pow(10, 3), mu_pN(Nd = 2*math.pow(10, 16))))
#print(rho_N(1*math.pow(10, 18), math.pow(10, 18), math.pow((1.07*math.pow(10,10)), 2)/ math.pow(10, 18)))