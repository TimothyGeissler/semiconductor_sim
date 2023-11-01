import sys
import math
import mpmath

Q = 1.602177*math.pow(10, -19)
ni = 1.07*math.pow(10, 10) #Intrinsic carrier concentration
kB = 1.380649*math.pow(10, -23)#0.00008617333262
T = 300
epsilon_r = 11.68 #rel permitivity of Si
epsilon_o = 8.854*math.pow(10, -14) #permitivity of vacuum

V_PN = 0.65# Bias Voltage, PN Junction DiodeSS

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
    #return ((Q * Na) / (2 *  epsilon_o * epsilon_r)) * W_dP(Na, Nd)**2
    return (Nd/(Na + Nd))*(v_bi(Na, Nd) - V_PN)

# Builtin voltage - n side
def V_biN(Na, Nd):
    #return ((Q * Nd) / (2 *  epsilon_o * epsilon_r)) * W_dN(Na, Nd)*W_dN(Na, Nd)
    return (Na/(Na + Nd)) * (v_bi(Na, Nd) - V_PN)

# P-side depletion width (cm)
def W_dP(Na, Nd):
    return math.sqrt((2*epsilon_o*epsilon_r / Q) * (Nd/(Na*(Na + Nd))) * (v_bi(Na, Nd) - V_PN))

# N-side depletion width (cm)
def W_dN(Na, Nd):
    return math.sqrt((2*epsilon_o*epsilon_r / Q) * (Na/(Nd*(Na + Nd))) * (v_bi(Na, Nd) - V_PN)) 

# Total depletiuon width
def W_d(Na, Nd):
    return W_dN(Na, Nd) + W_dP(Na, Nd)

# Max electric field
def E_max(Na, Nd):
    #return ((-Q * Na)/(epsilon_o * epsilon_r)) * W_dP(Na, Nd)
    return -math.sqrt(((2 * Q)/(epsilon_o * epsilon_r)) * ((Na * Nd)/(Na + Nd)) * (v_bi(Na, Nd)))

# Electron distribution p-side depletion region
def n_p_dep(Na, Nd, x):
    try:
        return ((ni**2)/Na) * math.exp((Q * V_PN) / (kB * T)) * math.exp((((Q**2)*Na)/(2*epsilon_o*epsilon_r*kB*T))*(W_dP(Na, Nd) + x)**2)
    except:
        return 0 #valid for 0 < x < w_dP

# Hole distribution p-side depletion region
def p_p_dep(Na, Nd, x):
    return Na*math.exp(-(((Q**2)*Na)/(2*epsilon_o*epsilon_r*kB*T))*(W_dP(Na, Nd) + x)**2)

# Electron distribution n-side depletion region
def n_n_dep(Na, Nd, x):
    return Nd*math.exp(-(((Q**2)*Nd)/(2*epsilon_o*epsilon_r*kB*T))*(W_dN(Na, Nd) - x)**2)

# Hole distribution n-side depletion region
def p_n_dep(Na, Nd, x):
    return ((ni**2)/Nd) * math.exp((Q * V_PN) / (kB * T)) * math.exp((((Q**2)*Nd)/(2*epsilon_o*epsilon_r*kB*T))*(W_dN(Na, Nd) - x)**2)

# Electron concentration at -W_d,P (edge of p-side depletion region)
def nP_depletion_edge(Na, Nd):
    return Nd*math.exp(-((Q * (v_bi(Na, Nd) - V_PN)) / (kB * T)))

# Hole concentration at W_d,N (edge of n-side depletion region)
def pN_depletion_edge(Na, Nd):
    return Na * math.exp(-((Q * (v_bi(Na, Nd) - V_PN)) / (kB * T)))

# Minority Carrier (n) recombination Diffusion length (p side)
def L_rec_n(Na, Nd):
    return math.sqrt(diff_coeff(T, mu_nP(Na, Nd)) * tau_rec_n(Nd))

# Minority Carrier (p) recombination Diffusion length ((n side))
def L_rec_p(Na, Nd):
    return math.sqrt(diff_coeff(T, mu_pN(Na, Nd)) * tau_rec_p(Na))

# Minority carrier (n) generation diff length (p side)
def L_gen_n(Na, Nd):
    return math.sqrt(diff_coeff(T, mu_nP(Na, Nd)) * tau_rec_n(Nd) * 75)

# Minority carrier (p) generation diff length (n side)
def L_gen_p(Na, Nd):
    return math.sqrt(diff_coeff(T, mu_pN(Na, Nd)) * tau_rec_p(Na) + 75)

# Hole diffusion current density, n-side, at W_p,N - long Base
def J_pdiffN_LB(Na, Nd):
    return ((Q * diff_coeff(T, mu_pN(Na, Nd)) * ni**2) / (L_rec_p(Na, Nd) * Nd)) * (math.exp((Q * V_PN) / (kB * T)) - 1)

# electron diffusion saturation current - Short Base
#def Is_diff_p_short(Na, Nd, A):
    #return ((Q * diff_coeff(T, mu_nP(Na, Nd)) * ni**2 * A) / (L_gen_p(Na, Nd, mu_nP(Na, Nd)) * Na)) * mpmath.coth((Wn - W_dN(Na, Nd)) / L_gen_p(Na, Nd))

# the hole diffusion saturation curret (n sidem fwd bias)
#def Is_diff_n_gen(Na, Nd, A):
    #return ((Q * diff_coeff(T, mu_pN(Na, Nd)) * ni**2 * A) / ((math.pow(10, -5) - W_dP(Na, Nd)) * Na))

# hole diffusion saturation current
def Is_diff_N(Na, Nd, A):
    return (Q * diff_coeff(T, mu_pN(Na, Nd)) * (ni**2) * A) / (L_rec_p(Na, Nd) * Nd)

# electron diffusion current at the edge of the p-side depletion region (short base)
def In_x_diff_p(Na, Nd, A):
    return Is_diff_p_short(Na, Nd, A) * (math.exp((Q * V_PN) / (kB * T)) - 1)

# hole diffusion current at edge of n-side depletion region (long base, fwd bias)
def Ip_x_diff_N(Na, Nd, A):
    return ((Q * diff_coeff(T, mu_pN(Na, Nd)) * (ni**2) * A) / (L_rec_p(Na, Nd) * Nd)) * (math.exp((Q * V_PN) / (kB * T)) - 1)

# depletion capacitance
def C_pn_dep(Na, Nd):
    return math.sqrt(((Q * epsilon_o * epsilon_r) / 2) * ((Na * Nd) / (Na + Nd)) * (1/(v_bi(Na, Nd) - V_PN)))

# quasi-neutral diffusion capacitance p-side (long base)
def C_pn_diff_P_long(Na, Nd, A):
    return (((Q**2) * A * (ni**2) * L_rec_n(Na, Nd)) / (kB * T * Na)) * math.exp((Q * V_PN) / (kB * T))

# quasi-neutral diffusion capacitance n-side (long base)
def C_pn_diff_N_long(Na, Nd, A):
    return (((Q**2) * A * (ni**2) * L_rec_p(Na, Nd)) / (kB * T * Nd)) * math.exp((Q * V_PN) / (kB * T))

def C_pn_diff_long(Na, Nd, A):
    C_diff_N_long = (((Q**2) * A * (ni**2) * L_rec_p(Na, Nd)) / (kB * T * Nd)) * math.exp((Q * 0) / (kB * T))
    C_diff_P_long = (((Q**2) * A * (ni**2) * L_rec_n(Na, Nd)) / (kB * T * Na)) * math.exp((Q * 0) / (kB * T))
    return (C_diff_N_long + C_diff_P_long) * math.exp((Q * V_PN) / (kB * T))

# Space-charge recombination time
def tau_rec_SCR(Na, Nd):
    return (tau_rec_n(Nd) + tau_rec_p(Na)) / 2

# Space-charge recombination Current
def I_S_SCR(Na, Nd, A):
    return (Q * ni * W_d(Na, Nd) * A) / (2 * tau_rec_SCR(Na, Nd))

def diff_current_density(Na, Nd, VPN):
    mu_pN(0, Na)
    D_N = diff_coeff(300, mu_pN(0, na))

    tau_N = tau_rec_n(Nd)
    if (VPN < 0): #REV BIAS
        tau_N = tau_N * 75

    L_N = math.sqrt(tau_N * D_N)

    j_px_diff_N =  ((Q * D_N * ni**2) / (L_N * Nd)) * math.exp((Q * VPN) / (kB * T) - 1)

    


## PRINT STATEMENTS ##
def comp_concentration(Na, Nd):
    # Compensated n-type  
    print("### N - TYPE SEMICONDUCTOR ###")
    print("\n## For compensated n-type (cm^-3) ##\nMajority nN=" + "{:E}".format(compensated_nN(Na, Nd)) + "\tMinority pN=" + "{:E}".format(compensated_pN(Na, Nd)))
    mu_ntype(Nd)
    print("\n## Resistivity (n-type) (Ohm.cm) ##\nrho_n=" + str(rho_N(Nd, compensated_nN(Na, Nd), compensated_pN(Na, Nd))))
    print("\n## Drift velocity (n-type)\nElectron drift v_n,drift=" + "{:E}".format(v_ndrift(E, mu_nN(Na, Nd))) + "\tHole drift v_p,drift=" + "{:E}".format(v_pdrift(E, mu_pN(Na, Nd))))
    print("\n## Diffusion Coefficient Compensated (cm^2/s) ##\nD_p,N=" + str(diff_coeff(T, mu_pN(Na, Nd))) + "\tD_n,N=" + str(diff_coeff(T, mu_nN(Na, Nd))))
    print("\n## Diffusion Coefficient (cm^2/s) ##\nD_p,N=" + str(diff_coeff(T, mu_pN(0, Nd))) + "\tD_n,N=" + str(diff_coeff(T, mu_nN(0, Nd))))
    print("\n## Recombination Lifetime (Sec) ##\ntau_rec,P= " + "{:E}".format(tau_rec_p(Na)) + "\ttau_rec,N=" + "{:E}".format(tau_rec_n(Nd)))
    
    
    # Compensated p-type
    print("### P - TYPE SEMICONDUCTOR ###")
    print("\n## For compensated p-type (cm^-3) ##\nMajority pP=" + "{:E}".format(compensated_pP(Na, Nd)) + "\tMinority nP=" + "{:E}".format(compensated_nP(Na, Nd)))
    mu_ptype(Na)
    print("\n## Resistivity (p-type) (Ohm.cm) ##\nrho_p=" + str(rho_P(Na, compensated_nP(Na, Nd), compensated_pP(Na, Nd))))
    print("\n## Drift velocity (p-type)\nElectron drift v_n,drift=" + "{:E}".format(v_ndrift(E, mu_nP(Na, Nd))) + "\tHole drift v_p,drift=" + "{:E}".format(v_pdrift(E, mu_pP(Na, Nd))))
    print("\n## Diffusion Coefficient Compensated (cm^2/s) ##\nD_p,N=" + str(diff_coeff(T, mu_pN(Na, Nd))) + "\tD_n,N=" + str(diff_coeff(T, mu_nN(Na, Nd))))
    print("\n## Diffusion Coefficient (cm^2/s) ##\nD_p,N=" + str(diff_coeff(T, mu_pN(Na, 0))) + "\tD_n,N=" + str(diff_coeff(T, mu_nN(Na, 0))))
    print("\n## Recombination Lifetime (sec) ##\ntau_rec,P= " + "{:E}".format(tau_rec_p(Na)) + "\ttau_rec,N=" + "{:E}".format(tau_rec_n(Nd)))

def print_pn_junc(Na, Nd, A):
    # Print PN Junction data
    print("\n## PN JUNCTION ##\n")
    print("\n## Builtin voltage (V) ##\nV_bi,P=" + str(V_biP(Na, Nd)) + "\tV_bi,N=" + str(V_biN(Na, Nd)) + "\tV_bi=" + str(v_bi(Na, Nd)))
    print("\n## Depletion Region Width (cm) ##\nW_d,p=" + str(W_dP(Na, Nd)) + "\tW_d,n=" + str(W_dN(Na, Nd)) + "\tW_d=" + str(W_d(Na, Nd)))
    print("\n## Carrier Concentration at depletion edge (cm^-3) ##\nnP(-W_d,P, V_PN)=" + "{:E}".format(nP_depletion_edge(Na, Nd)) + "\tpN(W_d,N, V_PN)=" + "{:E}".format(pN_depletion_edge(Na, Nd)))
    print("\n## Max Electric field (x=" + str(x) + ") (V/cm) ##\nE_max=" + str(E_max(Na, Nd)))
    print("\n## Carrier concentration (Depletion Region) x=" + str(x) + " ##\nn_p(x)=" + str(n_p_dep(Na, Nd, x)) + "\tp_p(x)=" + str(p_p_dep(Na, Nd, x))
          + "\tn_n(x)=" + str(n_n_dep(Na, Nd, x)) + "\tp_n(x)=" + str(p_n_dep(Na, Nd, x)))
    print("\n## Hole Diffusion Current Density @ edge of n-side depletion region (Long Base) J_pxdiff,N(W,dN , V_PN)=" + str(J_pdiffN_LB(Na, Nd)))
    print("\n## Hole diffusion current @ edge of n-side depletion region - Long base (A) ##\nI_p.x.diff,N=" + str(Ip_x_diff_N(Na, Nd, A)))
    print("\n## Depletion Capacitance (F) ##\nC_pn.dep=" + str(C_pn_dep(Na, Nd)))
    print("\n## Long Base quasi-neutral diffusion capacitance (F) ##\nC_pn.diff,P(Vpn)=" + str(C_pn_diff_P_long(Na, Nd, A)) + "\tC_pn.diff,N(Vpn)=" + str(C_pn_diff_N_long(Na, Nd, A)) + "\tC_pn.diff(Vpn)=" + str(C_pn_diff_long(Na, Nd, A)))
    print("\n## Carrier diffusion saturation current, long base (A) ##\nI_S.diff,N=" + str(Is_diff_N(Na, Nd, A)))
    print("\n## Space-charge recombination current (A) (Fwd Bias) ##\nI_s.SCR=" + str(I_S_SCR(Na, Nd, A)))

def mu_ptype(Na):
    print("\n## Low-field bulk mobility (p-type compensated) (cm^2/V.s) ## \nMajority mu_p,P=" + str(mu_pP(Na, Nd)) + "\tMinority mu_n,P=" + str(mu_nP(Na, Nd)))

def mu_ntype(Nd):
    print("\n## Low-field bulk mobility (n-type compensated) (cm^2/V.s) ## \nMajority mu_n,N=" + str(mu_nN(Na, 0)) + "\tMinority mu_p,N=" + str(mu_pN(Na, 0)))
    print("\n## Low-field bulk mobility (n-type) (cm^2/V.s) ## \nMajority mu_n,N=" + str(mu_nN(0, Nd)) + "\tMinority mu_p,N=" + str(mu_pN(0, Nd)))


## USAGE ##

Na = 3*math.pow(10, 17)
Nd = 5*math.pow(10, 16)

# PN Junction variables
x = W_dN(Na, Nd)
A = 7*math.pow(10, -5) #diode area A_d

E = 2.5*math.pow(10, 3) #for drift velocity


comp_concentration(Na, Nd)
print_pn_junc(Na, Nd, A)
#print(j_pxdrift(pN(2*math.pow(10,16)), 2.5*math.pow(10, 3), mu_pN(Nd = 2*math.pow(10, 16))))
#print(rho_N(1*math.pow(10, 18), math.pow(10, 18), math.pow((1.07*math.pow(10,10)), 2)/ math.pow(10, 18))

print(str(L_rec_n(Na, Nd)) + ", " + str(L_rec_p(Na, Nd)))
#print(str(Is_diff_p_short(Na, Nd, 7*math.pow(10, -5))))