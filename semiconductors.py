######################################################
#                                                    #
#           BASIC SEMICONDUCTOR EQUATIONS            #
#                                                    #
######################################################

import math

Q = 1.602177*math.pow(10, -19)
ni = 1.07*math.pow(10, 10) #Intrinsic carrier concentration
kB = 1.380649*math.pow(10, -23)#0.00008617333262
T = 300
epsilon_r = 11.68 #rel permitivity of Si
epsilon_o = 8.854*math.pow(10, -14) #permitivity of vacuum

class SEMI:
    def __init__(self, Nd, Na):
        self.Nd = Nd
        self.Na = Na

    def pP(self):
        pp = (0.5)*((self.Na - self.Nd) + math.sqrt(math.pow((self.Na - self.Nd), 2) + 4*ni*ni))
        print("p_P=" + str(pp) + " cm^-3")
        return pp
    
    def nN(self):
        nn = (0.5)*((self.Nd - self.Na) + math.sqrt(math.pow((self.Nd - self.Na), 2) + 4*ni*ni))
        print("n_N=" + str(nn) + " cm^-3")
        return nn

    def pN(self):
        pn = (2*(ni**2))/(self.Nd + math.sqrt((self.Nd)*(self.Nd) + 4*(ni**2)))
        #pn = (ni**2) / self.pP()
        print("p_N=" + str(pn) + " cm^-3")
        return pn
    
    def nP(self):
        np = (ni**2) / self.pP()
        print("n_P=" + str(np) + " cm^-3")
        return np

    def tau_rec_n(self):
        tau = 1/((7.8*math.pow(10,-13)) * self.Nd + (1.8*math.pow(10,-31))*(self.Nd*self.Nd))
        print("tau_rec,N=" + str(tau) + " s")
        return tau
    
    def tau_rec_p(self):
        tau = 1/((3.45*10**-12) * self.Na + (9.5*10**-32)*(self.Na*self.Na))
        print("tau_rec,P=" + str(tau) + " s")
        return tau
    
    def mu_pP(self):
        mu = 49.7 + 418.3/(1 + math.pow(((self.Na + self.Nd)/(1.6*math.pow(10, 17))), 0.7))
        print("mu_p,P=" + str(mu) + " cm^2/V.s")
        return mu

    #  Minority electron low-field bulk mobility
    def mu_nP(self):
        mu = 232 + 1180/(1 + math.pow(((self.Na + self.Nd)/(8*math.pow(10, 16))), 0.9))
        print("mu_n,P=" + str(mu) + " cm^2/V.s")
        return mu
    
    def mu_nN(self):
        mu = 92 + 1268/(1 + math.pow(((self.Na + self.Nd)/(1.3*math.pow(10, 17))), 0.91))
        print("mu_n,N=" + str(mu) + " cm^2/V.s")
        return mu

    #  Minority holes low-field bulk mobility
    def mu_pN(self):
        mu = 130 + 370/(1 + math.pow(((Na + Nd)/(8*math.pow(10, 17))), 1.25));
        print("mu_p,N=" + str(mu) + " cm^2/V.s")
        return mu
    
    def rho_N(self):
        rho = 1/(Q * (self.mu_nN() * self.nN() + self.mu_pN() * self.pN()))
        print("rho_N=" + str(rho) + " ohm.cm")
        return rho

    # Resistivity of a p-type Si (Ohm.cm)
    def rho_P(self):
        rho = 1/(Q * (self.mu_nP() * self.nP() + self.mu_pP() * self.pP()))
        print("rho_P=" + str(rho) + " ohm.cm")
        return rho
    
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
    
Na = 0#1*math.pow(10, 17)
Nd = 5*math.pow(10, 16)

sil = SEMI(Nd, Na)
sil.L_rec_pN()