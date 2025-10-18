import numpy as np
import constantes as ct

## Question 1.a/
def PIRE(P,G): # en dB
    """ Puissance isotrope rayonnée équivalente en dB """
    EIRP = []
    for i in P:
        EIRP.append(10*np.log10(i)+G)
    return np.array(EIRP)

## Question 1.b/
def D(Re, altitude_sat,Angle_elevation): #en m
    """ Distance satellite en m en se basant sur le rayon de la Terre,
    L'altitude de la station, l'altitude du satellite et l'angle d'élévation"""
    Dist =[]
    for i in Angle_elevation:
        phi = 90-i-np.arcsin((Re/(Re+altitude_sat))*np.cos(np.deg2rad(i))) 
        Dist.append(np.sqrt(Re**2+(Re+altitude_sat)**2-2*Re*(Re+altitude_sat)*np.cos(np.deg2rad(phi)))*10**(-3)) # en km
    return np.array(Dist)

## Question 1.c/
def LFS (Re, altitude_sat, Angle_elevation ,c, f):
    """ Pertes en espace libre en dB"""
    lamb = c/f
    D = D(Re, altitude_sat, Angle_elevation)
    perte_es = []
    for i in D:
        perte_es.append(20*np.log10(lamb/(4*np.pi*i)))
    return np.array(perte_es)

## Question 1.d/
def att_pluie(kH, kV, aH, aV, E, polar, R):
    """Atténuation spécifique pluie [dB/km] pour polarisation circulaire (~45°).
       E : angle d'élévation [°]
       polar : angle de polarisation [°]
       R : taux de pluie [mm/h]
    """
    pluie = []
    for i in E:
        i_rad = np.deg2rad(i)
        pol_rad = np.deg2rad(polar)
        k_coef = 0.5 * (kH + kV + (kH - kV)*(np.cos(i_rad))**2 * np.cos(2*pol_rad))
        alpha_coef = 0.5 * (kH*aH + kV*aV + (kH*aH - kV*aV)*(np.cos(i_rad))**2 * np.cos(2*pol_rad))
        pluie.append(10*np.log10(k_coef * (R**alpha_coef)))
    return np.array(pluie)

## Question 1.e/
def Le(kH, kV, aH, aV, h_s, E, R, polar, f, Lat):
    """Calcul de la longueur effective de cellule de pluie"""
    h_R = 4500 + 360 # h0 + 360 m
    Ls=[]
    Lg=[]
    for e in E:
        erad = np.deg2rad(e)
        Ls1 = (h_R-h_s)/np.sin(erad) #en m
        Ls.append(Ls1) #en m
        Lg.append(Ls*np.cos(erad)) #en m 
    gamma_r = att_pluie(kH, kV, aH, aV, E, polar, R)
    r0_01 = []
    for i in range(len(gamma_r)):
        r0_01.append(1/(1+0.78*np.sqrt((Lg[i]*gamma_r[i])/f)-0.38*(1-np.exp(-2*Lg[i]))))
    L_E = []
    for j in range(len(r0_01)):
        zeta = np.arctan((h_R-h_s)/(Lg[j]*r0_01[j]))
        erad = np.deg2rad(E[j])
        if zeta>E:
            L_R = (Lg[j]*r0_01[j])/np.cos(erad) #en m
        else:
            L_R = (h_R-h_s)/np.sin(erad) #en m
        if Lat < 36 :
            chi = 36 - abs(Lat)
        else:
            chi = 0
        v0_01 = 1/(1+np.sqrt(np.sin(erad))*(31*(1-np.exp(-(180*E/np.pi/(1+chi)))))*np.sqrt(L_R*gamma_r[j])/f**2-0.45)
        L_E = L_R * v0_01 #en m
    return np.array(L_E)

## Question 1.f/
def A001(kH, kV, aH, aV, h_s, E, R, polar, f, Lat):
    """Affaiblissement prévu dépassé pour 0,01 % d'une année moyenne"""
    L_E = Le(kH, kV, aH, aV, h_s, E, R, polar, f, Lat)
    Attpluie = att_pluie(kH, kV, aH, aV, E, polar, R)
    A_001 = []
    for i in range(len(L_E)):
        A_001.append(Attpluie[i] * L_E[i])
    return np.array(A001)

## Question 1.g/
def A1(kH, kV, aH, aV, h_s, E, R, polar, f, Lat): 
    " Atténuation due à la pluie dépassée pendant 1%, Beta = 0 comme p=1%"
    A_001 = A001(kH, kV, aH, aV, h_s, E, R, polar, f, Lat)
    A1=[]
    for i in A_001: 
        A1.append(i*(1/0.01)**(-0.655-0.045*np.log(i)))
    return np.array(A1)

## Question 1.h/    
def perte_polar (Ar_stat, Ar_sat,polar):
    """ Perte de polarisations où Ar_stat est le ratio axial du récepteur,
    Ar_sat est le ratio axial de l'émetteur et polar, l'angle de polarisation"""
    pol = []
    for i in Ar_sat:
        pol.append(10*np.log10(0.5 + (4*Ar_stat*i + (Ar_stat**2-1)*(i**2-1)*np.cos(2*polar))/(2*(Ar_stat**2+1)*(i**2+1))))
    return np.array(pol)

## Question 1.i/
def perte_depoint (c, f, D, dep):
    lamb = c/f
    theta_3dB = 70*lamb/D
    depoint=-12*(dep/theta_3dB)**2
    return depoint

## Question 1.j/
def Gain_pertedep(c, f, D, efficacite, dep):
    """ Gain de l'antenne en prenant en compte la perte de dépointage"""
    lamb = c/f
    Dmax = ((np.pi*D)/lamb)**2
    G = 10*np.log10(efficacite*Dmax)
    depoint = perte_depoint(c, f, D, dep)
    G_depoint = G + depoint
    return G_depoint
    
## Question 1.k/
def Pre(c, f, D, efficacite, dep, P, G, kH, kV, aH, aV, E, Ar_stat, Ar_sat, polar, R, Angle_elevation, Re, h_sat):
    """
    Calcule la puissance reçue (porteuse) en dBW au niveau du récepteur.
    """
    EIRP = PIRE(P, G)
    G_depoint = Gain_pertedep(c, f, D, efficacite, dep)
    Attpluie = att_pluie(kH, kV, aH, aV, E, polar, R)
    p_pol = perte_polar(Ar_stat, Ar_sat, polar)
    Dist = D(Re, h_sat, Angle_elevation)
    lamb = c/f
    p_EL = LFS(Dist, lamb)
    P_reçue = []
    for i in range(len(EIRP)):    
        P_reçue.append(EIRP[i] + G_depoint[i] + p_EL[i] - Attpluie[i] + p_pol[i] + G_depoint[i])
    return np.array(P_reçue)

## Question 1.l/
def temp_entree_recep(TLNA, Tmx, GLNA, Tif, Gmx):
    """ Température de bruit totale à l'entrée du récepteur en K"""
    Trx = TLNA +Tmx/(10**(GLNA/10)) + Tif/(10**((GLNA+Gmx)/10))
    return Trx

## Question 1.m/
def temperature_bruit_entrée_antenne(Tsky,Tgd, Tm, P):
    """
    Calcule la température équivalente de bruit en sortie de l'antenne (en K).
    """
    # Calcul de T_a_in
    L_atm = 10 ** (P / 10)
    Ta_in = (Tsky / L_atm) + Tm * (1 - 1 / L_atm) + Tgd
    return (Ta_in)

## Question 1.n/
def temperature_bruit_sortie_antenne(c, f, D, efficacite, Tsky, Tgd, Tm, P):
    # Calcul de T_a_out Gain 
    lamb = c/f
    Dmax = ((np.pi*D)/lamb)**2
    G = 10*np.log10(efficacite*Dmax)
    GRx_ratio = 10 ** (G / 10)
    Ta_in = temperature_bruit_entrée_antenne(Tsky, Tgd, Tm, P)
    Ta_out = (1 / 720) * Ta_in * (GRx_ratio ** 2)
    return Ta_out

## Question 1.o/
def temperature_entrée_récepteur(c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable): 
    Ta_out = temperature_bruit_sortie_antenne(c, f, D, efficacite, Tsky, Tgd, Tm, P)
    Trx = Tf + Ta_out / (10**(L_cable/10))
    return Trx

## Question 1.p/
def n0(k_bolt,c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable): 
    "Densité de puissance de bruit spectrale"
    Trx = temperature_entrée_récepteur(c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable)
    n0 = k_bolt * Trx 
    return n0

## Question 1.q/
def C_N0_ratio(dep, G, kH, kV, aH, aV, E, Ar_stat, Ar_sat, polar, R, Angle_elevation, Re, h_sat, k_bolt,c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable):
    Prx = Pre(c, f, D, efficacite, dep, P, G, kH, kV, aH, aV, E, Ar_stat, Ar_sat, polar, R, Angle_elevation, Re, h_sat)
    N0 = n0(k_bolt,c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable)
    ratio =[]
    for p in Prx:
        ratio.append(p - 10*np.log10(N0))
    return np.array(ratio)

## Question 1.r/
def Eb_N0_ratio(Rb, dep, G, kH, kV, aH, aV, E, Ar_stat, Ar_sat, polar, R, Angle_elevation, Re, h_sat, k_bolt,c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable):
    CN0_ratio = C_N0_ratio(dep, G, kH, kV, aH, aV, E, Ar_stat, Ar_sat, polar, R, Angle_elevation, Re, h_sat, k_bolt,c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable)
    EbN0_ratio = []
    for i in CN0_ratio:
        EbN0_ratio.append(i - 10*np.log10(Rb))
    return np.array(EbN0_ratio)

## Question 1.s/
def bande_passante(alpha, Rb):
    "Bande passante du signal"
    B = 3*Rb * (1 + alpha)
    return B

## Question 1.t/
def eff_spectrale(c_05, Rb, alpha):
    "Efficacité spectrale"
    eta = c_05 * Rb / bande_passante(alpha, Rb)
    return eta

## Question 1.u/
def marge(Rb, alpha):
    "Marge bande passante"
    B = bande_passante(alpha, Rb)
    marge = B/2 - 3*Rb*(1+alpha)
    return marge

## Question 1.v/
def debit_max(alpha, Rb, marge_impose):
    "Débit maximum"
    B = bande_passante(alpha, Rb)
    debit_max = 3*(B-2*marge_impose)/(1+alpha)
    return debit_max

#Question 2.a
def trace_courbes_paramétriques(E,debit,Ptx,marge): 
    # Courbes pour les différents angles d'élévation
    plt.figure(figsize=(8,5))
    for i in E:
        plt.plot(marge,debit[i])
    plt.xlabel("Débit d'information pour les différentes angles d'élévation (Mbits/s")
    plt.ylabel("Marge système restante (dB)")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Courbes pour différentes puissances Ptx
    plt.figure(figsize=(8,5))
    for j in Ptx:
        plt.plot(marge,debit[j])
    
    plt.xlabel("Débit d'information pour les différentes puissances(Mbit/s)")
    plt.ylabel("Marge système restante (dB)")
    plt.legend()
    plt.grid(True)
    plt.show()
    
if __name__ == "__main__":
    P = ct.Ptx
    Dist = D(ct.Re, ct.h_sat, ct.Angle_elevation)
    lam = lamb(ct.c, ct.f)
    G = ct.Gain_bord
    EIRP = PIRE(P,G)
    Perte_EL = LFS(Dist, lam)
    Att_spe_pluie = att_pluie(ct.kH, ct.kV, ct.aH, ct.aV, ct.Angle_elevation, ct.Angle_ellipse, ct.I_precip)
    polar = perte_polar(ct.AR_stat, ct.AR_sat, ct.Angle_ellipse)
    depoint = perte_depoint(lam, ct.Diametre_antenne, ct.Depointage)
    Ls1 = Ls(ct.hr, ct.h_stat, ct.Angle_elevation)
    Lg1= Lg(Ls1, ct.Angle_elevation)
    r001= fc(Lg1, ct.I_precip, ct.f)
    zetha= f_zetha(ct.hr, ct.h_stat, Lg1, r001)
    Lr1 = Lr(Lg1, r001, ct.Angle_elevation, zetha, ct.hr, ct.h_sat)
    v001 = f_v001(ct.Angle_elevation, Lr1, Att_spe_pluie, ct.f)
    Le1 = Le(Lr1, v001)
    A_001 = A001(Att_spe_pluie, Le1)
    A_1 = A1(A_001)
    #P_reçue = Pre(EIRP, ct.Gmx, Perte_EL, A_1, polar, depoint)
    print('PIRE=', EIRP)
    print('Distance=', Dist)
    print('Pertes espace libre=', Perte_EL)
    print('Atténuation spé pluie=', Att_spe_pluie)
    print('pertes polarisation =', polar)
    print('pertes dépointage=', depoint)
    print("atténuation due à la pluie dépassée pendant 0,01%=", A_001)
    #print("Puissance de la porteuse reçue par le récepteur est:",P_reçue)
    #print ("La température de bruit en sortie de l'antenne est:",Ta_out)
    #print("La température équivalente de bruit totale à l'entrée du récepteur est",Trx)
    #print(" La densité de puissance de bruit reçue par le récepteur est", Drx)
