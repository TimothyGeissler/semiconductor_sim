import math
Q = 1.602177*math.pow(10, -19)
ni = 1.07*math.pow(10, 10) #Intrinsic carrier concentration
kB = 1.380649*math.pow(10, -23)#0.00008617333262
T = 300
epsilon_si = 11.68 #rel permitivity of Si
epsilon_ox = 3.9
epsilon_o = 8.854*math.pow(10, -14) #permitivity of vacuum

class MOSCAP:


    def __init__(self, Na, X_ox, Q_it, Q_f):
        self.Na = Na
        self.X_ox = X_ox
        self.Q_it = Q_it
        self.Q_f = Q_f

    # Oxide capacitance
    def C_ox(self):
        c = (epsilon_ox * epsilon_o) / self.X_ox
        print("C_ox=" + str(c) + " F/cm^2")
        return c
    
    # Fermi potential
    def phi_fb(self):
        fb = ((kB * T) / Q) * math.log(self.Na / ni)
        print("phi_FB=" + str(fb) + " V")
        return fb
    
    # Contact potential difference - Al <-> p-type
    def phi_pm(self):
        pm = -0.51165 - ((kB * T) / Q) * math.log(self.Na / ni)
        print("phi_PM=" + str(pm) + " V")
        return pm
    
    #Space-charge voltage @ Vtn
    def V_sc_tn(self):
        print("V_sc(V_TN)=" + str(self.phi_fb()*2) + " V (invert if p-type)")
        return self.phi_fb()*2

    # Depletion width @ Vtn
    def W_d_vtn(self):
        w = math.sqrt(((2 * epsilon_o * epsilon_si) / (Q * self.Na)) * self.V_sc_tn())
        print("W_d,B(Vtn)=" + str(w) + " cm")
        return w
    
    # Max depletion width ( = W_d,inv,HF)
    def W_d_max(self):
        w = math.sqrt((4 * epsilon_o * epsilon_si * self.phi_fb()) / (Q * self.Na))
        print("W_d,max=" + str(w) + " cm")
        return w

    # Charge Density in SCR @ Vtn
    def Q_sc_tn(self):
        q = -Q * self.Na * self.W_d_vtn()
        print("Q_sc(V_tn)=" + str(q) + " C/cm^2")
        return q
    
    def Q_sc_vgb(self, vgb):
        beta = Q / (kB * T)
        ldb = math.sqrt((epsilon_o * epsilon_si * kB * T) / (Q**2 * self.Na))
        qsc = ((math.sqrt(2) * epsilon_o * epsilon_si * kB * T) / (Q * ldb)) * ((math.exp(-beta * vgb) + beta * vgb - 1) + (self.nP() / self.Na) * (math.exp(beta * vgb) - beta * vgb - 1))**(0.5)
        print("Q_sc(V_GB)=" + str(qsc) + " C/cm^2")
        return qsc
    
    # substrate inversion charge (inversion bias range)
    def Q_inv_b(self, vgb):
        qinv = -self.C_ox() * (vgb - self.V_tn())
        print("Q_inv,b=" + str(qinv) + " C/cm^2")
        return qinv
    
    def v_sc_inv(self, vgb):
        ldb = math.sqrt((epsilon_o * epsilon_si * kB * T) / (Q**2 * self.Na))
        #v = 2 * self.phi_fb() + ((2 * kB * T) / Q) * math.log(1 - ((Q * ldb) / (math.sqrt(2) * epsilon_si * epsilon_o * kB * T)) * self.Q_inv_b(vgb))
        v = 2 * self.phi_fb()
        print("V_sc_inv(Vgb)=" + str(v) + " V")
        return v

    # Threshold voltage
    def V_tn(self):
        # If Q_it = Q_f = 0
        # = self.V_fb() + 2 * self.phi_fb() + (math.sqrt(4 * Q * epsilon_si * epsilon_o * self.Na * self.phi_fb()) / self.C_ox())
        # Generalized:
        v = self.phi_pm() - ((self.Q_f + self.Q_it) / self.C_ox()) + 2 * self.phi_fb() + (math.sqrt(4 * Q * epsilon_si * epsilon_o * self.Na * self.phi_fb()) / self.C_ox())
        print("V_TN=" + str(v) + " V")
        return v
    
    # Oxide Voltage @ flatband
    def V_ox(self):
        v = -((self.Q_f + self.Q_it) / self.C_ox())
        print("V_ox(V_FB)=" + str(v) + " V")
        return v

    # Flatband voltage
    def V_fb(self):
        v = self.phi_pm() + self.V_ox()
        print("V_fb=" + str(v) + " V")
        return v
    
    # Transverse E field at Vtn
    def E_ox(self):
        trans_field = (math.sqrt(4 * Q * epsilon_si * epsilon_o * self.Na * self.phi_fb()) - self.Q_f - self.Q_it) / (epsilon_o * epsilon_ox)
        print("E_ox=" + "{:E}".format(trans_field) + "V/cm")
        return trans_field

    # Voltage across gate-oxide layer at vtn
    def V_ox(self):
        v = self.X_ox * self.E_ox()
        print("V_ox=" + str(v) + "V")
        return v
    
    # Transverse E field in substrate @ gate-oxide/substract interface at vtn
    def E_b(self):
        eb = -(self.Q_sc_tn() / (epsilon_si * epsilon_o))
        print("E_bulk=" + "{:E}".format(eb) + " V/cm")
        return eb
 
    # Mobility of holes
    def mu_pP(self):
        mu = 49.7 + 418.3/(1 + math.pow(((self.Na)/(1.6*math.pow(10, 17))), 0.7))
        print("mu_pP=" + str(mu) + " cm^2/Vs")
        return mu
    
    # Mobility of electrons
    def mu_nP(self):
        mu = 232 + 1180/(1 + math.pow(((self.Na)/(8*math.pow(10, 16))), 0.9))
        print("mu_nP=" + str(mu) + " cm^2/Vs")
        return mu
    
    # Majority carrier concentration
    def pP(self):
        p = (0.5)*(self.Na + math.sqrt((self.Na**2) + 4*(ni)*(ni)))
        print("p_P=" + str(p) + " cm^-1")
        return p

    # Minority carrier concentration
    def nP(self):
        n = (2*(ni)*(ni))/(self.Na + math.sqrt((self.Na**2) + 4*(ni)*(ni)))
        print("n_P=" + str(n) + " cm^-1")
        return n

    # Diffusion coefficient
    def D_pP(self):
        D = ((kB * T) / Q) * self.mu_pP()
        print("D_pP=" + str(D) + " cm^2/s")
        return D
    
    # Condutivity
    def sigma(self):
        sig = Q * (self.mu_nP() * self.nP() + self.mu_pP() * self.pP())
        print("sigma=" + str(sig) + " (ohm.cm)^1")
        return sig

    # Dielectric relaxation time
    def tau_diel(self):
        tau = (epsilon_o * epsilon_si) / self.sigma()
        print("tau_diel=" + str(tau) + " s")
        return tau
    
    # Debye Length
    def L_D(self):
        L = math.sqrt(self.D_pP() * self.tau_diel())
        print("L_D,P=" + str(L) + "")
        return L
    
    # SCR capacitance @ Low Freq, flatband
    def C_sc_LF(self):
        c = (epsilon_o * epsilon_si) / self.L_D()
        print("C_sc,LF(V_FB)=" + str(c) + " F/cm^2")
        return c
    
    # SCR capacitance @ High freq, Vtn
    def C_sc_HF(self):
        c = (epsilon_o * epsilon_si) / self.W_d_vtn()
        print("C_sc,HF(V_tn)=" + str(c) + " F/cm^2")
        return c

    # Gate-to-bulk capacitance @ low Freq, Flatband
    def C_gb_LF(self):
        c = (self.C_ox() * (self.C_sc_LF() + self.Q_it)) / (self.C_ox() + self.C_sc_LF() + self.Q_it)
        print("C_gb,LF(V_FB)=" + str(c) + " F/cm^2")
        return c
    
    def C_gb_LF_vgb(self, vgb):
        c = (self.C_ox() * self.C_sc_LF_vgb(vgb)) / (self.C_ox()  + self.C_sc_LF_vgb(vgb))
        print("C_gb_LF(Vgb)=" + str(c) + " C/cm^2")
        return c

    def C_sc_LF_vgb(self, vgb):
        beta = Q / (kB * T)
        ldb = math.sqrt((epsilon_o * epsilon_si * kB * T) / (Q**2 * self.Na))
        csc = ((epsilon_o * epsilon_si) / (math.sqrt(2) * ldb)) * (((1 - math.exp(-beta*self.v_sc_inv(vgb))) + (self.nP() / self.Na) * (math.exp(beta*self.v_sc_inv(vgb)) - 1)) / (math.sqrt((math.exp(-beta*self.v_sc_inv(vgb)) + beta*self.v_sc_inv(vgb) - 1) + (self.nP() / self.Na) * (math.exp(beta*self.v_sc_inv(vgb)) - beta*self.v_sc_inv(vgb) - 1))))
        print("C_sc_LF(Vgb)=" + str(csc) + " C/cm^2")
        return csc

    # Gate-to-bulk capacitance @ high Freq, Vtn
    def C_gb_HF(self):
        c = (self.C_ox() * (self.C_sc_HF() + self.Q_it)) / (self.C_ox() + self.C_sc_HF() + self.Q_it)
        print("C_gb,HF(V_Tn)=" + str(c) + " F/cm^2")
        return c
    
    # Gate-to-bulk capacitance @ high freq, inversion
    def C_gb_HF_inv(self):
        c = math.pow((self.X_ox/(epsilon_o * epsilon_ox)) + (self.W_d_max() / (epsilon_o * epsilon_si)), -1)
        print("C_gb,HF_inv(V_GB)=" + str(c) + " F/cm^2")
        return c

Na = 3.5*math.pow(10, 17)
X_ox = (35)*math.pow(10, -8) # Angstrom * 10^-8 = cm
N_f = 2.8 * math.pow(10, 10)
Q_it = 0
Q_f = Q * N_f


mos = MOSCAP(Na, X_ox, Q_it, Q_f)
mos.phi_pm()