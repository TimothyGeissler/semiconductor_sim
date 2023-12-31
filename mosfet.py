import math
#import numpy as np
#import matplotlib
#import matplotlib.pyplot as plt

Q = 1.602177 * math.pow(10, -19)
ni = 1.07 * math.pow(10, 10)  # Intrinsic carrier concentration
kB = 1.380649 * math.pow(10, -23)  # 0.00008617333262
T = 300
epsilon_si = 11.68  # rel permitivity of Si
epsilon_ox = 3.9
epsilon_o = 8.854 * math.pow(10, -14)  # permitivity of vacuum


class MOSFET:

    def __init__(self, Na, X_ox, W_ch, L_ch, Q_it, Q_f):
        self.Na = Na
        self.X_ox = X_ox
        self.W_ch = W_ch
        self.L_ch = L_ch

        self.Q_it = Q_it
        self.Q_f = Q_f

    # Forward Transconductance (sat)
    # Approx for lambda*vds < 1
    def gm_sat(self, vgs):
        g = mu_n_ch * self.C_ox() * (self.W_ch / self.L_ch) * (vgs - self.V_TN())
        print("g_m_sat=" + "{:E}".format(g) + " mS")
        return g

    # saturation-range drain-to-source resistance
    def r_ds_sat(self, Va, Vgs):
        r = abs(Va) / self.Id_sat(Vgs)
        print("r_ds_sat=" + str(r) + " Ohm")
        return r

    #Cutoff frequency
    def f_T(self, vgs):
        f = ((3 * mu_n_ch) / (4 * math.pi * self.L_ch**2)) * (vgs - self.V_TN())
        print("f_T=" + "{:E}".format(f) + "Hz")
        return f
    # Saturation drain current (Vbs = 0) - SAH
    def Id_sat_lvl1(self, Vgs):
        i_d = mu_n_ch * self.C_ox() * (self.W_ch / (2 * self.L_ch)) * ((Vgs - self.V_TN()) ** 2)
        print("Id_sat=" + "{:E}".format(i_d) + " A")
        return i_d
    
    def VD_sat_lvl2(self, vgs):
        vdsat = vgs - self.V_FB() - 2 * self.phi_fb() + ((Q * epsilon_o * epsilon_si * self.Na) / (self.C_ox() ** 2)) * (1 - math.sqrt(1 + ((2 * self.C_ox()**2)/ (Q * epsilon_o * epsilon_si * self.Na)) * (vgs - self.V_FB() - 0)))
        print("V_DSat=" + str(vdsat) + " V - lvl2")
        return vdsat

    # Saturation drain current (Vbs = 0) - Ihantola-Moll
    def Id_sat_lvl2(self, vgs):
        id = mu_n_ch * self.C_ox() * (self.W_ch / self.L_ch) * ((vgs - self.V_FB() - 2 * self.phi_fb() - (0.5) * self.VD_sat_lvl2(vgs)) * self.VD_sat_lvl2(vgs) - (2/3) * (math.sqrt(2 * Q * epsilon_o * epsilon_si * self.Na) / self.C_ox()) * ((2 * self.phi_fb() + self.VD_sat_lvl2(vgs) - 0) ** (3/2) - (2 * self.phi_fb() - 0) ** (3/2)))
        print("I_D,sat=" + str(id) + " A - lvl2")
        return id

    # bulk body-effect coefficient
    def gamma_bn(self):
        gamma = math.sqrt(2 * Q * epsilon_si * epsilon_o * self.Na) / self.C_ox()
        print("gamma_bn=" + str(gamma) + " V^(1/2)")
        return gamma

    # Oxide capacitance
    def C_ox(self):
        c = (epsilon_ox * epsilon_o) / self.X_ox
        print("C_ox=" + str(c) + " F/cm^2")
        return c

    # Gate-oxide capacitance
    def C_gate_ox(self):
        c = self.C_ox() * self.W_ch * self.L_ch
        print("C_ox (gate-oxide cap) =" + str(c) + "F")
        return c

    # Output Curve Sah model (lvl2)
    def output_curve_lvl2(self, vtn):
        plot = [[0 for j in range(5)] for i in range(5)]  # Fix the dimensions of the array
        x = np.linspace(0, 5, 5)
        for vgs in range(1, 6):  # Adjust the range to start from 0
            for vds in range(1, 6):
                vdsat = vgs - vtn
                if vds >= vdsat:  # saturation region
                    i_d = mu_n_ch * self.C_ox() * (self.W_ch / (2 * self.L_ch)) * ((vgs - vtn) ** 2)
                else:  # linear region
                    i_d = mu_n_ch * self.C_ox() * (self.W_ch / self.L_ch) * ((vgs - vtn) * vds - 0.5 * (vds ** 2))
                # Adjust indices to fix the IndexError
                plot[vgs-1][vds-1] = i_d

        # np.savetxt("plot.txt", plot, delimiter=',')
        fig, ax = plt.subplots(figsize=(10, 6))
        count = 0
        for vgs in plot:
            label = 'VGS=' + str(count) + "V"
            ax.plot(x, vgs, label=label)
            count+=1

        ax.legend()
        ax.set_xlabel('V_DS')
        ax.set_ylabel('I_d')
        ax.set_title('I_d(V_DS) Output Curve - Sah Model')
        plt.show()
        return 0

    # SQRT output Curve Sah model (lvl2)
    def output_curve_lvl2_sqrt(self, vtn):
        plot = [0 for j in range(10)]  # Fix the dimensions of the array
        x = np.linspace(0, 10, 10)
        for vgs in range(0, 10):
            vdsat = vgs - vtn
            if 3 >= vdsat:  # saturation region
                i_d = mu_n_ch * self.C_ox() * (self.W_ch / (2 * self.L_ch)) * ((vgs - vtn) ** 2)
            else:  # linear region
                i_d = mu_n_ch * self.C_ox() * (self.W_ch / self.L_ch) * ((vgs - vtn) * 3 - 0.5 * (3 ** 2))
            # Adjust indices to fix the IndexError
            plot[vgs - 1] = math.sqrt(i_d)

        # np.savetxt("plot.txt", plot, delimiter=',')
        fig, ax = plt.subplots(figsize=(10, 6))

        ax.plot(x, plot, label="V_DS = 3V")
        ax.legend()
        ax.set_xlabel('V_DS')
        ax.set_ylabel('sqrt(I_d)')
        ax.set_title('I_d(V_DS) Output Curve - Sah Model')
        plt.show()
        return 0

    # Output curve I-M Model (lvl1)
    def output_curve_lvl1(self, vtn):
        plot = [[0 for j in range(5)] for i in range(5)]  # Fix the dimensions of the array
        x = np.linspace(0, 5, 5)
        for vgs in range(1, 6):  # Adjust the range to start from 0
            for vds in range(1, 6):
                vdsat = vgs - self.V_FB() - 2 * self.phi_fb() + ((Q * epsilon_si * epsilon_o * self.Na) / (self.C_ox()**2))*(1 - math.sqrt(1 + ((2 * self.C_ox()**2) / (Q * epsilon_si * epsilon_o * self.Na)) * (vgs - self.V_FB())))
                if vds >= vdsat:  # saturation region
                    i_d = mu_n_ch * self.C_ox() * (self.W_ch / self.L_ch) * ((vgs - self.V_FB() - 2 * self.phi_fb() - 0.5 * vdsat) * vdsat) - (2/3)*(math.sqrt(2 * Q * epsilon_si * epsilon_o * self.Na)) / (self.C_ox()) * ((2*self.phi_fb() + vdsat)**(3/2) - (2*self.phi_fb())**(3/2))
                else:  # linear region
                    i_d = mu_n_ch * self.C_ox() * (self.W_ch / self.L_ch) * ((vgs - self.V_FB() - 2 * self.phi_fb() - 0.5 * vds) * vds) - (2/3)*((math.sqrt(2 * Q * epsilon_si * epsilon_o * self.Na)) / (self.C_ox())) * ((2*self.phi_fb() + vds)**(3/2) - (2*self.phi_fb())**(3/2))
                # Adjust indices to fix the IndexError
                plot[vgs-1][vds-1] = -i_d

        # np.savetxt("plot.txt", plot, delimiter=',')
        fig, ax = plt.subplots(figsize=(10, 6))
        count = 0
        for vgs in plot:
            label = 'VGS=' + str(count) + "V"
            ax.plot(x, vgs, label=label)
            count+=1

        ax.legend()
        ax.set_xlabel('V_DS')
        ax.set_ylabel('I_d')
        ax.set_title('I_d(V_DS) Output Curve - Ihantola-Moll Model')
        plt.show()
        return 0

    # Threshold voltage @ V_BS = 0
    def V_TN(self):
        vtn = self.V_FB() + 2 * self.phi_fb() + self.gamma_bn() * math.sqrt(2 * self.phi_fb())
        print("V_TN=" + str(vtn) + " V")
        return vtn

    # Flatband Voltage (V_BS = 0)
    def V_FB(self):
        vfb = self.V_OX() + self.phi_pm()
        print("V_FB=" + str(vfb) + " V")
        return vfb

    # Voltage along oxide (longitudinal) (V_BS = 0)
    def V_OX(self):
        v = -((self.Q_f + self.Q_it) / self.C_ox())
        print("V_OX=" + str(v) + " V")
        return v

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

    # Threshold voltage, Q_it = Q_f = 0
    # def V_tn(self):
    # v = self.phi_pm() - ((self.Q_f + self.Q_it) / self.C_ox()) + 2*self.phi_fb() + math.sqrt(4 * Q * epsilon_o * epsilon_si * self.phi_fb())/self.C_ox()
    # print("V_TN=" + str(v) + " V")
    # return v

    # Oxide Voltage @ flatband
    # def V_ox(self):
    # v = -((self.Q_f + self.Q_it) / self.C_ox())
    # print("V_ox(V_FB)=" + str(v) + " V")
    # return v

    # Flatband voltage
    # def V_fb(self):
    # v = self.phi_pm() + self.V_ox()
    # print("V_fb=" + str(v) + " V")
    # return v

    # Mobility of holes
    def mu_pP(self):
        mu = 49.7 + 418.3 / (1 + math.pow(((self.Na) / (1.6 * math.pow(10, 17))), 0.7))
        print("mu_pP=" + str(mu) + " cm^2/Vs")
        return mu

    # Mobility of electrons
    def mu_nP(self):
        mu = 232 + 1180 / (1 + math.pow(((self.Na) / (8 * math.pow(10, 16))), 0.9))
        print("mu_nP=" + str(mu) + " cm^2/Vs")
        return mu
    
    # Channel electron mobility
    def mu_nch(self, theta, engma, mu_ltf):
        mu = mu_ltf / (1 + ((self.vgs - self.V_TN()) / (theta * self.X_ox()))**engma)
        print("mu_n,ch=" + str(mu) + " cm^2/Vs")
        return mu

    # Majority carrier concentration
    def pP(self):
        p = (0.5) * (self.Na + math.sqrt((self.Na ** 2) + 4 * (ni) * (ni)))
        print("p_P=" + str(p) + " cm^-1")
        return p

    # Minority carrier concentration
    def nP(self):
        n = (2 * (ni) * (ni)) / (self.Na + math.sqrt((self.Na ** 2) + 4 * (ni) * (ni)))
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

    # Gate-to-bulk capacitance @ high Freq, Vtn
    def C_gb_HF(self):
        c = (self.C_ox() * (self.C_sc_HF() + self.Q_it)) / (self.C_ox() + self.C_sc_HF() + self.Q_it)
        print("C_gb,HF(V_Tn)=" + str(c) + " F/cm^2")
        return c

    # Gate-to-bulk capacitance @ high freq, inversion
    def C_gb_HF_inv(self):
        c = math.pow((self.X_ox / (epsilon_o * epsilon_ox)) + (self.W_d_max() / (epsilon_o * epsilon_si)), -1)
        print("C_gb,HF_inv(V_GB)=" + str(c) + " F/cm^2")
        return c
    
    # Gate to source overlap capacitance
    def C_gs_ov(self, L_ov):
        c = (epsilon_o * epsilon_ox * self.W_ch * L_ov) / self.X_ox
        print("C_gs.ov=" + str(c) + "F/cm^2")
        return c
    
    # saturation bias range gate-to-source capacitance
    def c_gs_sat(self, L_ov):
        c = (2/3)*self.C_ox() * self.W_ch * self.L_ch + self.C_ox() * self.W_ch * L_ov
        print("C_gs=" + str(c) + "F")
        return c

    # Channel-length-modulation coefficient
    def lambda_n(self, va):
        l = -1/va
        print("lambda_n=" + str(l) + " V^-1")
        return l



Na = 4 * math.pow(10, 17)
X_ox = (50) * math.pow(10, -8) #angstrom to cm
W_ch = (15) * math.pow(10, -4) #um to cm
L_ch = (5) * math.pow(10, -4)

mu_n_ch = 250
N_f = 3 * math.pow(10, 10)
Q_it = 0
Q_f = Q * N_f

mos = MOSFET(Na, X_ox, W_ch, L_ch, Q_it, Q_f)
mos.Id_sat_lvl2(1.5)