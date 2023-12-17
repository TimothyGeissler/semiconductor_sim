import math

Q = 1.602177*math.pow(10, -19)
ni = 1.07*math.pow(10, 10) #Intrinsic carrier concentration
kB = 1.380649*math.pow(10, -23)#0.00008617333262
T = 300
epsilon_r = 11.68 #rel permitivity of Si
epsilon_o = 8.854*math.pow(10, -14) #permitivity of vacuum

class PN:
    def __init__(self, Na, Nd, V_PN, W, A):
        self.Nd = Nd
        self.Na = Na
        self.V_PN = V_PN
        self.W = W
        self.A = A
    
    def nP(self):
        np = (ni**2) / self.pP()
        print("n_P=" + str(np) + " cm^-3")
        return np
    
    def pN(self):
        pn = (ni**2) / self.nN()
        print("p_N=" + str(pn) + " cm^-3")
        return pn
    
    def pP(self):
        pp = self.Na
        print("p_P=" + str(pp) + " cm^-3")
        return pp
    
    def nN(self):
        nn = self.Nd
        print("n_N=" + str(nn) + " cm^-3")
        return nn
    
    def v_bi(self):
        v = ((kB * T)/Q) * math.log(self.Na * self.Nd / (ni**2))
        print("V_bi=" + str(v) + " V")
        return v
    
    # Builtin voltage - p side
    def V_biP(self):
        #return ((Q * Na) / (2 *  epsilon_o * epsilon_r)) * W_dP(Na, Nd)**2
        vbp = (self.Nd/(self.Na + self.Nd))*(self.v_bi() - self.V_PN)
        print("V_bi,P=" + str(vbp) + " V")
        return vbp

    # Builtin voltage - n side
    def V_biN(self):
        #return ((Q * Nd) / (2 *  epsilon_o * epsilon_r)) * W_dN(Na, Nd)*W_dN(Na, Nd)
        vbn = (self.Na/(self.Na + self.Nd)) * (self.v_bi() - self.V_PN)
        print("V_bi,N=" + str(vbn) + " V")
        return vbn
    
    # P-side depletion width (cm)
    def W_dP(self):
        wdp = math.sqrt((2*epsilon_o*epsilon_r / Q) * (self.Nd/(self.Na*(self.Na + self.Nd))) * (self.v_bi() - self.V_PN))
        print("W_d,P=" + str(wdp) + " cm")
        return wdp
    
    # N-side depletion width (cm)
    def W_dN(self):
        wdn = math.sqrt((2*epsilon_o*epsilon_r / Q) * (self.Na/(self.Nd*(self.Na + self.Nd))) * (self.v_bi() - self.V_PN)) 
        print("W_d,N=" + str(wdn) + " cm")
        return wdn
    
    # Total depletiuon width
    def W_d(self):
        wd = self.W_dN() + self.W_dP()
        print("W_d=" + str(wd) + " cm")
        return wd
    
        # Diffusion coefficient
    def D_pP(self):
        D = ((kB * T) / Q) * self.mu_pP()
        print("D_pP=" + str(D) + " cm^2/s")
        return D
    
    # Diffusion coefficient
    def D_nP(self):
        D = ((kB * T) / Q) * self.mu_nP()
        print("D_nP=" + str(D) + " cm^2/s")
        return D
    
    # Diffusion coefficient
    def D_pN(self):
        D = ((kB * T) / Q) * self.mu_pN()
        print("D_pN=" + str(D) + " cm^2/s")
        return D
    
    # Diffusion coefficient
    def D_nN(self):
        D = ((kB * T) / Q) * self.mu_nN()
        print("D_nN=" + str(D) + " cm^2/s")
        return D
    
        # Minority Carrier (n) recombination Diffusion length (p side)
    def L_rec_nP(self):
        L =  math.sqrt(self.D_nP() * self.tau_rec_p())
        print("L_rec,N=" + str(L) + " cm")
        return L

    # Minority Carrier (p) recombination Diffusion length ((n side))
    def L_rec_pN(self):
        L =  math.sqrt(self.D_pN() * self.tau_rec_n())
        print("L_rec,P=" + str(L) + " cm")
        return L

    # Minority carrier (n) generation diff length (p side)
    def L_gen_nP(self):
        L =  math.sqrt(self.nP() * self.tau_rec_p() * 75)
        print("L_gen,N=" + str(L) + " cm")
        return L

    # Minority carrier (p) generation diff length (n side)
    def L_gen_pN(self):
        L =  math.sqrt(self.D_pN() * self.tau_rec_n() * 75)
        print("L_gen,P=" + str(L) + " cm")
        return L
    
    def tau_rec_n(self):
        tau = 1/((7.8*math.pow(10,-13)) * self.Nd + (1.8*math.pow(10,-31))*(self.Nd*self.Nd))
        print("tau_rec,N=" + str(tau) + " s")
        return tau
    
    def tau_rec_p(self):
        tau = 1/((3.45*10**-12) * self.Na + (9.5*10**-32)*(self.Na*self.Na))
        print("tau_rec,P=" + str(tau) + " s")
        return tau
    
    def mu_pP(self):
        mu = 49.7 + 418.3/(1 + math.pow(((self.Na)/(1.6*math.pow(10, 17))), 0.7))
        print("mu_p,P=" + str(mu) + " cm^2/V.s")
        return mu

    #  Minority electron low-field bulk mobility
    def mu_nP(self):
        mu = 232 + 1180/(1 + math.pow(((self.Na)/(8*math.pow(10, 16))), 0.9))
        print("mu_n,P=" + str(mu) + " cm^2/V.s")
        return mu
    
    def mu_nN(self):
        mu = 92 + 1268/(1 + math.pow(((self.Nd)/(1.3*math.pow(10, 17))), 0.91))
        print("mu_n,N=" + str(mu) + " cm^2/V.s")
        return mu
    
    #  Minority holes low-field bulk mobility
    def mu_pN(self):
        mu = 130 + 370/(1 + math.pow(((self.Nd)/(8*math.pow(10, 17))), 1.25));
        print("mu_p,N=" + str(mu) + " cm^2/V.s")
        return mu
    
    # hole diffusion current density for a long-base n-side quasi-neutral region
    def J_pxdiff_N_LB(self):
        J = ((Q * self.D_pN() * (ni**2)) / (self.L_rec_pN() * self.Nd)) * (math.exp((Q * self.V_PN) / (kB * T)) - 1)
        print("J_p.x,diff,N=" + str(J) + " A/cm^2")
        return J

    # electron diffusion current density for a long-base n-side quasi-neutral region
    def J_nxdiff_P_LB(self):
        j = ((Q * self.D_nP() * (ni**2)) / (self.L_rec_nP() * self.Na)) * (math.exp((Q * self.V_PN) / (kB * T)) - 1)
        print("J_n.x,diff,P=" + str(j) + " A/cm^2")
        return j
    
    #def V_B(self):
        #E_g_si = 1.124
        #vb = -60 * ((E_g_si/ 1.1) ** (3/2)) * ((10**16 / self.Nd) **(3/4))

    # If base << 1 then S-B else L-B
    def short_base_P(self, WP):
        base = (WP - self.W_dP()) / self.L_rec_nP()
        print("Short-base ratio (p)=" + str(base))

    def short_base_N(self, WN):
        base = (WN - self.W_dN()) / self.L_rec_pN()
        print("Short-base ratio (n)=" + str(base))

Na = 2*math.pow(10, 18)
Nd = 8*math.pow(10, 17)
W = 0.1 * math.pow(10, -4)
A = 7 * math.pow(10, -5)

VPN = 0.35

pn = PN(Na, Nd, VPN, W, A)
pn.short_base_P(0.1*math.pow(10, -4))
