import math
import numpy as np


Q = 1.602177 * math.pow(10, -19)
ni = 1.07 * math.pow(10, 10)  # Intrinsic carrier concentration
kB = 1.380649 * math.pow(10, -23)  # 0.00008617333262
T = 300
epsilon_si = 11.68  # rel permitivity of Si
epsilon_ox = 3.9
epsilon_o = 8.854 * math.pow(10, -14)  # permitivity of vacuum


class BJT:

    def __init__(self, Nd_E, Na_B, Nd_C, A, W_B, W_C, W_E, V_BC, V_BE):
        self.Nd_E = Nd_E
        self.Nd_C = Nd_C
        self.Na_B = Na_B
        self.A = A

        self.W_B = W_B
        self.W_C = W_C
        self.W_E = W_E

        self.V_BC = V_BC

    def p_E(self):
        pe = ni**2/self.Nd_E
        #print("p_E (Min. conc. @ th. eq.)=" + str(pe) + " cm^-3")
        return pe
    
    def n_E(self):
        pe = self.Nd_E
        #print("n_E (Maj. conc. @ th. eq.)=" + str(pe) + " cm^-3")
        return pe

    def n_B(self):
        nb = ni**2/self.Na_B
        #print("n_B (Min. conc. @ th. eq.)=" + str(nb) + " cm^-3")
        return nb
    
    def p_B(self):
        nb = self.Na_B
        #print("p_B (Maj. conc. @ th. eq.)=" + str(nb) + " cm^-3")
        return nb
    
    def p_C(self):
        pc = ni**2/self.Nd_C
        #print("p_C (Min. conc. @ th. eq.)=" + str(pc) + " cm^-3")
        return pc
    
    def n_C(self):
        pc = self.Nd_C
        #print("n_C (Maj. conc. @ th. eq.)=" + str(pc) + " cm^-3")
        return pc
    
    #  Minority holes low-field bulk mobility - Emitter
    def mu_pE(self):
        mu = 130 + 370/(1 + math.pow(((self.n_E() + self.p_E())/(8*math.pow(10, 17))), 1.25));
        print("mu_p,E=" + str(mu) + " cm^2/V.s")
        return mu
    #  Minority holes low-field bulk mobility - Collector
    def mu_pC(self):
        mu = 130 + 370/(1 + math.pow(((self.n_C() + self.p_C())/(8*math.pow(10, 17))), 1.25));
        print("mu_p,C=" + str(mu) + " cm^2/V.s")
        return mu

    #  Minority electron low-field bulk mobility - Base
    def mu_nB(self):
        mu = 232 + 1180/(1 + math.pow(((self.n_B() + self.p_B())/(8*math.pow(10, 16))), 0.9))
        print("mu_n,B=" + str(mu) + " cm^2/V.s")
        return mu

    # Minority carrier diffusivity
    def D_pE(self):
        d = (kB * T * self.mu_pE())/Q
        print("D_p,E=" + str(d))
        return d
    
    def D_pC(self):
        d = (kB * T * self.mu_pC())/Q
        print("D_p,C=" + str(d))
        return d
    
    def D_nB(self):
        d = (kB * T * self.mu_nB())/Q
        print("D_n,B=" + str(d))
        return d
    
    # Builtin voltage Emitter/Base
    def Vbi_BE(self):
        v = ((kB * T) / Q) * math.log((self.p_B() * self.n_E()) / ni**2)
        print("V_bi,BE=" + str(v) + " V")
        return v
    
    # Builtin voltage Base/Collector
    def Vbi_BC(self):
        v = ((kB * T) / Q) * math.log((self.Na_B * self.Nd_C) / ni**2)
        print("V_bi,BC=" + str(v) + " V")
        return v
    
    # Depletion region in E
    def Wd_E(self):
        wd_e =  math.sqrt(((2 * epsilon_o * epsilon_si) / (Q * self.Nd_E)) * (self.Na_B / (self.Na_B + self.Nd_E)) * (self.Vbi_BE() - V_BE))
        print("W_d,E=" + str(wd_e) + " cm (Dep. Width in Emitter)")
        return wd_e
    
    # Depletion region in B on E side
    def Wd_B_E(self):
        wd_b_e = math.sqrt(((2 * epsilon_o * epsilon_si) / (Q * self.Na_B)) * (self.Nd_E / (self.Na_B + self.Nd_E)) * (self.Vbi_BE() - V_BE))
        print("W_d,B/E=" + str(wd_b_e) + " cm (Dep. Width on Emitter side in Base)")
        return wd_b_e

    # Dep. width on Base/Emitter junction (both sides) - TH Eq (pg 9.26)
    def Wd_BE_total(self):
        wd_be = self.Wd_B_E() + self.Wd_E()
        print("W_d,BE=" + str(wd_be) + " cm")
        return wd_be
    
    # Depletion region in Collector
    def Wd_C(self):
        wd_c = math.sqrt(((2 * epsilon_o * epsilon_si) / (Q * self.Nd_C)) * (self.Na_B / (self.Na_B + self.Nd_C)) * (self.Vbi_BC() - V_BC))
        print("W_d,C=" + str(wd_c) + " cm")
        return wd_c
    
    # Depletion region in B on C side
    def Wd_B_C(self):
        wd_b_c = math.sqrt(((2 * epsilon_o * epsilon_si) / (Q * self.Na_B)) * (self.Nd_C / (self.Na_B + self.Nd_C)) * (self.Vbi_BC() - V_BC))
        print("W_d,B/C=" + str(wd_b_c) + " cm (Dep. Width on Collector side in Base)")
        return wd_b_c

    # Dep. width on Base/Collector junction (both sides) - TH Eq (pg 9.27)
    def Wd_BC_total(self):
        wd_bc = math.sqrt(((2 * epsilon_o * epsilon_si * kB * T) / (Q**2)) * ((Na_B + Nd_C) / (Na_B * Nd_C)) * math.log((Na_B * Nd_C) / ni**2))
        print("W_d,BC=" + str(wd_bc) + " cm")

    def tau_rec_n_B(self):
        tau = 1/((3.45*math.pow(10, -12) * self.Na_B) + (9.5*math.pow(10, -32))*(self.Na_B**2))
        print("tau_rec.n,B=" + str(tau))
        return tau
    
    def tau_rec_p_E(self):
        tau = 1/((7.8*10**-13) * self.Nd_E + (1.8*10**-31)*(self.Nd_E**2))
        print("tau_rec.p,E=" + str(tau))
        return tau
    
    def tau_rec_p_C(self):
        tau = 1/((7.8*10**-13) * self.Nd_C + (1.8*10**-31)*(self.Nd_C**2))
        print("tau_rec.p,C=" + str(tau))
        return tau
    
    def tau_gen_n_B(self):
        tau = 75 * self.tau_rec_n_B()
        print("tau_gen.n,B=" + str(tau))
        return tau
    
    def tau_gen_p_E(self):
        tau = 75 * self.tau_rec_p_E()
        print("tau_gen.p,E=" + str(tau))
        return tau
    
    def tau_gen_p_C(self):
        tau = 75 * self.tau_rec_p_C()
        print("tau_gen.p,C=" + str(tau))
        return tau
    
    def L_pE(self):
        L = 0
        if (V_BE < 0 and V_BC < 0):
            # Cut-off bias
            L = math.sqrt(self.D_pE() * self.tau_gen_p_E())
            print("L_p,E=" + "{:E}".format(L) + " cm (Cut-Off)")
        elif (V_BE > 0 and V_BC < 0):
            # Fwd-active
            L = math.sqrt(self.D_pE() * self.tau_rec_p_E())
            print("L_p,E=" + "{:E}".format(L) + " cm (Fwd. Active)")
        elif (V_BE > 0 and V_BC > 0):
            # Saturation
            L = math.sqrt(self.D_pE() * self.tau_rec_p_E())
            print("L_p,E=" + "{:E}".format(L) + " cm (Saturation)")
        elif (V_BE < 0 and V_BC > 0):
            # Rev active
            L = math.sqrt(self.D_pE() * self.tau_gen_p_E())
            print("L_p,E=" + "{:E}".format(L) + " cm (Rev. Active)")
        return L
    
    def L_pC(self):
        L = 0
        if (V_BE < 0 and V_BC < 0):
            # Cut-off bias
            L = math.sqrt(self.D_pC() * self.tau_gen_p_C())
            print("L_p,C=" + "{:E}".format(L) + " cm (Cut-Off)")
        elif (V_BE > 0 and V_BC < 0):
            # Fwd-active
            L = math.sqrt(self.D_pC() * self.tau_rec_p_C())
            print("L_p,C=" + "{:E}".format(L) + " cm (Fwd. Active)")
        elif (V_BE > 0 and V_BC > 0):
            # Saturation
            L = math.sqrt(self.D_pC() * self.tau_rec_p_C())
            print("L_p,C=" + "{:E}".format(L) + " cm (Saturation)")
        elif (V_BE < 0 and V_BC > 0):
            # Rev active
            L = math.sqrt(self.D_pC() * self.tau_gen_p_C())
            print("L_p,C=" + "{:E}".format(L) + " cm (Rev. Active)")
        return L
        
    def L_nB(self):
        L = 0
        if (V_BE < 0 and V_BC < 0):
            # Cut-off bias
            L = math.sqrt(self.D_nB() * self.tau_gen_n_B())
            print("L_n,B=" + "{:E}".format(L) + " cm (Cut-Off)")
        elif (V_BE > 0 and V_BC < 0):
            # Fwd-active
            L = math.sqrt(self.D_nB() * self.tau_rec_n_B())
            print("L_n,B=" + "{:E}".format(L) + " cm (Fwd. Active)")
        elif (V_BE > 0 and V_BC > 0):
            # Saturation
            L = math.sqrt(self.D_nB() * self.tau_rec_n_B())
            print("L_n,B=" + "{:E}".format(L) + " cm (Saturation)")
        elif (V_BE < 0 and V_BC > 0):
            # Rev active
            L = math.sqrt(self.D_nB() * self.tau_gen_n_B())
            print("L_n,B=" + "{:E}".format(L) + " cm (Rev. Active)")
        return L
    
    def coth(self, x):
        return 1 / np.tanh(x)
    
    def tau_SCR_BE(self):
        tau_rec = (self.tau_rec_p_E() + self.tau_rec_n_B()) / 2 
        tau = 0
        if (V_BE < 0 and V_BC < 0):
            # Cut-off bias
            tau = 75 * tau_rec
            print("tau_SCR,BE=" + "{:E}".format(tau) + " s (Cut-Off)")
        elif (V_BE > 0 and V_BC < 0):
            # Fwd-active
            tau = tau_rec
            print("tau_SCR,BE=" + "{:E}".format(tau) + " s (Fwd. Active)")
        elif (V_BE > 0 and V_BC > 0):
            # Saturation
            tau = tau_rec
            print("tau_SCR,BE=" + "{:E}".format(tau) + " s (Saturation)")
        elif (V_BE < 0 and V_BC > 0):
            # Rev active
            tau = 75 * tau_rec
            print("tau_SCR,BE=" + "{:E}".format(tau) + " s (Rev. Active)")
        return tau

    
    # Common-emitter current gain
    def beta_F(self):
        b = (np.cosh((self.W_B - self.Wd_B_E() - self.Wd_B_C()) / self.L_nB()) + ((self.L_nB() * self.D_pE() * self.Na_B) / (self.L_pE() * self.D_nB() * self.Nd_E)) * self.coth((self.W_E - self.Wd_E()) / self.L_pE()) * np.sinh((self.W_B - self.Wd_B_E() - self.Wd_B_C()) / self.L_nB()) + ((self.L_nB() * self.Na_B * self.Wd_BE_total()) / (2 * ni * self.tau_SCR_BE() * self.D_nB())) * np.sinh((self.W_B - self.Wd_B_E() - self.Wd_B_C()) / self.L_nB()) * math.exp(-(Q * V_BE) / (2 * kB * T)) - 1)**(-1)
        print("beta_F=" + str(b))
        return b

    # common-base current gain
    def alpha_F(self):
        a = self.beta_F() / (1 + self.beta_F())
        print("alpha_F=" + str(a))
        return a
    
    # emitter injection efficiency
    def gamma_E(self):
        g = (1 + ((self.L_nB() * self.D_pE() * self.Na_B) / (self.L_pE() * self.D_nB() * self.Nd_E)) * (self.coth((self.W_E - self.Wd_E()) / self.L_pE()) / self.coth((self.W_B - self.Wd_B_E() - self.Wd_B_C()) / self.L_nB())) + ((self.L_nB() * self.Na_B * self.Wd_BE_total()) / (2 * ni * self.tau_SCR_BE() * self.D_nB())) * np.tanh((self.W_B - self.Wd_B_E() - self.Wd_B_C()) / self.L_nB()) * math.exp(-(Q * V_BE) / (2 * kB * T)))**(-1)
        print("gamma_E=" + str(g))
        return g
    
    # Base transport factor
    def alpha_T(self):
        alpha = 1/np.cosh((self.W_B - self.Wd_B_E() - self.Wd_B_C()) / self.L_nB())
        print("alpha_T=" + str(alpha))
        return alpha
    
    # minority-carrier transit time (base QNR, fwd active)
    def tau_tr_B(self):
        tau = math.log(np.cosh((self.W_B - self.Wd_B_E() - self.Wd_B_C()) / self.L_nB())) * self.tau_rec_n_B()
        print("tau_tr.n,B=" + str(tau) + " s")
        return tau
    
    #cutoff frequency
    def f_T(self):
        ft = 1/(2 * math.pi * self.tau_tr_B())
        print("f_T=" + "{:E}".format(ft) + " Hz")
        return ft

Nd_E = 2 * math.pow(10, 20) # Emmitter
Na_B = 2 * math.pow(10, 18) # Base
Nd_C = 2 * math.pow(10, 16) # Collector
A = ((20) * math.pow(10, -4)) * ((20) * math.pow(10, -4)) # um^2 to cm^2

W_E = (2) * math.pow(10, -4) # Emmitter Width to um
W_B = (0.75) * math.pow(10, -4) # Base Width to um
W_C = (10) * math.pow(10, -4) # Collector Width to um

V_BE = 0.65
V_BC = -8

bjt = BJT(Nd_E, Na_B, Nd_C, A, W_B, W_C, W_E, V_BC, V_BE)
bjt.D_nB()

#public_method_names = [method for method in dir(bjt) if callable(getattr(bjt, method)) if not method.startswith('_')]  # 'private' methods start from _
#for method in public_method_names:
    #getattr(bjt, method)()  # call