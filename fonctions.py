import numpy as np
import matplotlib.pyplot as plt
import constantes as ct

## Question 1.a/
def PIRE(P,G): # en dB
    """ Puissance isotrope rayonnée équivalente en dB
     P est à rentrer en W et G en dB """
    EIRP = 10*np.log10(P)+G
    return EIRP

## Question 1.b/
def D(Re, altitude_sat,Angle_elevation): #en km
    """ Distance satellite en m en se basant sur le rayon de la Terre RE,
    l'altitude du satellite et l'angle d'élévation en degré"""
    phi = 90-Angle_elevation-np.arcsin((Re/(Re+altitude_sat))*np.cos(np.deg2rad(Angle_elevation))) 
    Dist = np.sqrt(Re**2+(Re+altitude_sat)**2-2*Re*(Re+altitude_sat)*np.cos(np.deg2rad(phi)))*10**(-3) # en km
    return Dist

## Question 1.c/
def LFS (Re, altitude_sat, Angle_elevation ,c, f):
    """ Pertes en espace libre en dB
    Re le rayon de la terre
    l'altitude duc satellite 
    Angle_elevation en degré
    c la vitesse de la lumière
    f la fréquence en Hz
    """
    lamb = c/f
    Dist = D(Re, altitude_sat, Angle_elevation)
    perte_es = 20*np.log10(lamb/(4*np.pi*Dist*10**3)) # D en m
    return perte_es

## Question 1.d/
def att_pluie(kH, kV, aH, aV, E, polar, R):
    """Atténuation spécifique pluie [dB/km] pour polarisation circulaire (~45°).
       E : angle d'élévation [°]
       polar : angle de polarisation [°]
       R : taux de pluie [mm/h]
    """
    i_rad = np.deg2rad(E)
    pol_rad = np.deg2rad(polar)
    k_coef = 0.5 * (kH + kV + (kH - kV)*(np.cos(i_rad))**2 * np.cos(2*pol_rad))
    alpha_coef = 0.5 * (kH*aH + kV*aV + (kH*aH - kV*aV)*(np.cos(i_rad))**2 * np.cos(2*pol_rad))/k_coef
    pluie = k_coef * (R**alpha_coef)
    return pluie

## Question 1.e/
def Le(kH, kV, aH, aV, h_s, E, R, polar, f, Lat):
    """Calcul de la longueur effective de cellule de pluie"""
    h_R = 4500 + 360 # h0 + 360 m
    erad = np.deg2rad(E)
    Ls = (h_R-h_s)/np.sin(erad) #en m
    Lg = Ls*np.cos(erad) #en m 
    gamma_r = att_pluie(kH, kV, aH, aV, E, polar, R)
    r0_01 = 1/(1+0.78*np.sqrt((Lg*gamma_r)/f)-0.38*(1-np.exp(-2*Lg)))
    zeta = np.arctan((h_R-h_s)/(Lg*r0_01))
    if zeta>E:
        L_R = (Lg*r0_01)/np.cos(erad) #en m
    else:
        L_R = (h_R-h_s)/np.sin(erad) #en m
    if Lat < 36 :
        chi = 36 - abs(Lat)
    else:
        chi = 0
    v0_01 = 1/(1+np.sqrt(np.sin(erad))*(31*(1-np.exp(-(180*E/np.pi/(1+chi)))))*np.sqrt(L_R*gamma_r)/f**2-0.45)
    L_E = L_R * v0_01 #en m
    return L_E

## Question 1.f/
def A001(kH, kV, aH, aV, h_s, E, R, polar, f, Lat):
    """Affaiblissement prévu dépassé pour 0,01 % d'une année moyenne"""
    L_E = Le(kH, kV, aH, aV, h_s, E, R, polar, f, Lat)
    Attpluie = att_pluie(kH, kV, aH, aV, E, polar, R)
    A_001 = (Attpluie * L_E*10**(-3))
    return A_001

## Question 1.g/
def A1(kH, kV, aH, aV, h_s, E, R, polar, f, Lat): 
    " Atténuation due à la pluie dépassée pendant 1%, Beta = 0 comme p=1%"
    A_001 = A001(kH, kV, aH, aV, h_s, E, R, polar, f, Lat)
    A_1 = A_001*(1/0.01)**(-(0.655-0.045*np.log(A_001)))
    return A_1

## Question 1.h/    
def perte_polar (Ar_stat, Ar_sat,polar):
    """ Perte de polarisations où Ar_stat est le ratio axial du récepteur,
    Ar_sat est le ratio axial de l'émetteur et polar, l'angle de polarisation"""
    pol = 10*np.log10(0.5 + (4*Ar_stat*Ar_sat + (Ar_stat**2-1)*(Ar_sat**2-1)*np.cos(2*polar))/(2*(Ar_stat**2+1)*(Ar_sat**2+1)))
    return pol

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
def Pre(c, f, Diam, efficacite, dep, P, G, kH, kV, aH, aV, E, Ar_stat, Ar_sat, polar, R, Angle_elevation, Re, h_sat):
    """
    Calcule la puissance reçue (porteuse) en dBW au niveau du récepteur.
    """
    EIRP = PIRE(P, G)
    G_depoint = Gain_pertedep(c, f, Diam, efficacite, dep)
    Attpluie = att_pluie(kH, kV, aH, aV, E, polar, R)
    p_pol = perte_polar(Ar_stat, Ar_sat, polar)
    Dist = D(Re, h_sat, Angle_elevation)
    lamb = c/f
    p_EL = LFS(Dist, lamb, Angle_elevation, c, f)  
    P_reçue = EIRP + G_depoint - Attpluie + p_EL + p_pol 
    return P_reçue

## Question 1.l/
def temp_entree_recep(TLNA, Tmx, GLNA, Tif, Gmx):
    """ Température de bruit totale à l'entrée du récepteur en K"""
    Trx = TLNA +Tmx/(10**(GLNA/10)) + Tif/(10**((GLNA+Gmx)/10))
    return Trx

## Question 1.m/
def temperature_bruit_entree_antenne(Tsky,Tgd, Tm, Atm):
    """
    Calcule la température équivalente de bruit en sortie de l'antenne (en K).
    """
    # Calcul de T_a_in
    L_atm = 10 ** (Atm / 10)
    Ta_in = (Tsky / L_atm) + Tm/ L_atm**2 + Tgd
    return (Ta_in)

## Question 1.n/
def temperature_bruit_sortie_antenne(c, f, D, efficacite, Tsky, Tgd, Tm, P):
    # Calcul de T_a_out Gain 
    """Calcul de la Température de bruit en sortie d’antenne """
    theta_3dB = 70*c/(f*D)
    G_max = 36000/(theta_3dB**2)
    GRx = efficacite*G_max
    GRxdB = 10*np.log10(GRx)
    Ta_in = temperature_bruit_entree_antenne(Tsky, Tgd, Tm, P)
    Ta_out = 1/720*Ta_in*GRxdB*(theta_3dB)**2
    return Ta_out

## Question 1.o/
def temperature_eq_recep(c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable): 
    Ta_out = temperature_bruit_sortie_antenne(c, f, D, efficacite, Tsky, Tgd, Tm, P)
    Trx = Tf + Ta_out / (10**(L_cable/10))
    return Trx

## Question 1.p/
def n0(k_bolt,c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable): 
    "Densité de puissance de bruit spectrale"
    Trx = temperature_eq_recep(c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable)
    n0 = k_bolt * Trx 
    return n0

## Question 1.q/
def C_N0_ratio(dep, G, kH, kV, aH, aV, E, Ar_stat, Ar_sat, polar, R, Angle_elevation, Re, h_sat, k_bolt,c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable):
    Prx = Pre(c, f, D, efficacite, dep, P, G, kH, kV, aH, aV, E, Ar_stat, Ar_sat, polar, R, Angle_elevation, Re, h_sat)
    N0 = n0(k_bolt,c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable)
    ratio = (Prx - 10*np.log10(N0))
    return ratio

## Question 1.r/
def Eb_N0_ratio(Rb, dep, G, kH, kV, aH, aV, E, Ar_stat, Ar_sat, polar, R, Angle_elevation, Re, h_sat, k_bolt,c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable):
    CN0_ratio = C_N0_ratio(dep, G, kH, kV, aH, aV, E, Ar_stat, Ar_sat, polar, R, Angle_elevation, Re, h_sat, k_bolt,c, f, D, efficacite, Tsky, Tgd, Tm, P,Tf,L_cable)
    EbN0_ratio = CN0_ratio - 10*np.log10(Rb)
    return EbN0_ratio

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
    marge = B/2 - 3*Rb*(1+alpha)/2
    return marge

## Question 1.v/
def debit_max(alpha, Rb, marge_impose):
    "Débit maximum"
    B = bande_passante(alpha, Rb)
    debit_max = 3*(B-2*marge_impose)/(1+alpha)
    return debit_max

#Question 2.a
def trace_courbes_parametriques(E, Ptx, alpha, Rb_values):
    # Courbes pour les différents angles d'élévation
    plt.figure(figsize=(8,5))
    for i in E:
        marge_vals = [marge(Rb, alpha) for Rb in Rb_values]
        plt.plot(Rb_values, marge_vals)
    
    plt.xlabel("Débit d'information Rb (Mbit/s)")
    plt.ylabel("Marge système restante (dB)")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Courbes pour différentes puissances Ptx
    plt.figure(figsize=(8,5))
    for j in Ptx:
        marge_vals = [marge(Rb, alpha) + j/10 for Rb in Rb_values] 
        plt.plot(Rb_values, marge_vals)
    
    plt.xlabel("Débit d'information Rb (Mbit/s)")
    plt.ylabel("Marge système restante (dB)")
    plt.legend()
    plt.grid(True)
    plt.show()
if __name__ == "__main__":
    EIRP = PIRE(ct.P_test,ct.Gain_bord[9])
    Dist = D(ct.Re, ct.h_sat, ct.E_test)
    Loss_el = LFS(ct.Re, ct.h_sat, ct.E_test ,ct.c, ct.f)
    Att_spe_pluie = att_pluie(ct.kH, ct.kV, ct.aH, ct.aV, ct.E_test, ct.Angle_ellipse, ct.I_precip)
    Long_eff_pluie = Le(ct.kH, ct.kV, ct.aH, ct.aV, ct.h_stat, ct.E_test, ct.I_precip, ct.Angle_ellipse, ct.f, ct.Lat)
    A_001 = A001(ct.kH, ct.kV, ct.aH, ct.aV, ct.h_stat, ct.E_test, ct.I_precip, ct.Angle_ellipse, ct.f, ct.Lat)
    A_1 = A1(ct.kH, ct.kV, ct.aH, ct.aV, ct.h_stat, ct.E_test, ct.I_precip, ct.Angle_ellipse, ct.f, ct.Lat)
    p_pol = perte_polar (ct.AR_stat, ct.AR_sat[9], ct.Angle_ellipse)
    p_depoint = perte_depoint (ct.c, ct.f, ct.Diametre_antenne, ct.Depointage)
    G_depoint = Gain_pertedep(ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Depointage)
    P_reçue = Pre(ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Depointage, ct.P_test, ct.Gain_bord[9], ct.kH, ct.kV, ct.aH, ct.aV, ct.E_test, ct.AR_stat, ct.AR_sat[9], ct.Angle_ellipse, ct.I_precip, ct.E_test, ct.Re, ct.h_sat)
    Trxin = temp_entree_recep(ct.Tlna, ct.Tmx, ct.Glna, ct.Tif, ct.Gmx)
    Tain = temperature_bruit_entree_antenne(ct.Tsky, ct.Tgd, ct.Tm, ct.Perte_atm)
    Taout = temperature_bruit_sortie_antenne(ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test)
    Trx = temperature_eq_recep(ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx)
    N0 = n0(ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx)
    CN0_ratio = C_N0_ratio(ct.Depointage, ct.Gain_bord[9], ct.kH, ct.kV, ct.aH, ct.aV, ct.E_test, ct.AR_stat, ct.AR_sat[9], ct.Angle_ellipse, ct.I_precip, ct.E_test, ct.Re, ct.h_sat, ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx)
    EbN0_ratio = Eb_N0_ratio(ct.B_test, ct.Depointage, ct.Gain_bord[9], ct.kH, ct.kV, ct.aH, ct.aV, ct.E_test, ct.AR_stat, ct.AR_sat[9], ct.Angle_ellipse, ct.I_precip, ct.E_test, ct.Re, ct.h_sat, ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx)
    B = bande_passante(ct.Roll_off_factor, ct.B_test)
    eta = eff_spectrale(ct.Code_correcteur, ct.B_test, ct.Roll_off_factor)
    m = marge(ct.B_test, ct.Roll_off_factor)
    deb_max = debit_max(ct.Roll_off_factor, ct.B_test, ct.Marge)
    print("\n\nRésultats des questions :")
    print("Question 1.a/ \nEIRP =", EIRP, "dBW")
    print("Question 1.b/ \nDistance satellite station =", Dist, "km")
    print("Question 1.c/ \nPertes en espace libre =", Loss_el, "dB")
    print("Question 1.d/ \nAtténuation spécifique pluie =", Att_spe_pluie, "dB/km")
    print("Question 1.e/ \nLongueur effective de cellule de pluie =", Long_eff_pluie, "m")
    print("Question 1.f/ \nAffaiblissement prévu dépassé pour 0,01 % d'une année moyenne =", A_001, "dB")
    print("Question 1.g/ \nAtténuation due à la pluie dépassée pendant 1% =", A_1, "dB")
    print("Question 1.h/ \nPerte de polarisation =", p_pol, "dB")
    print("Question 1.i/ \nPerte de dépointage =", p_depoint, "dB")
    print("Question 1.j/ \nGain de l'antenne avec perte de dépointage =", G_depoint, "dB")
    print("Question 1.k/ \nPuissance reçue au niveau du récepteur =", P_reçue, "dBW")
    print("Question 1.l/ \nTempérature de bruit totale à l'entrée du récepteur =", Trxin, "K")
    print("Question 1.m/ \nTempérature équivalente de bruit en sortie de l'antenne =", Tain, "K")
    print("Question 1.n/ \nTempérature équivalente de bruit en sortie de l'antenne =", Taout, "K")
    print("Question 1.o/ \nTempérature équivalente de bruit au récepteur =", Trx, "K")
    print("Question 1.p/ \nDensité de puissance de bruit spectrale =", N0, "W/Hz")
    print("Question 1.q/ \nRapport C/N0 =", CN0_ratio, "dBHz")
    print("Question 1.r/ \nRapport Eb/N0 =", EbN0_ratio, "dB")
    print("Question 1.s/ \nBande passante du signal =", B, "Hz")
    print("Question 1.t/ \nEfficacité spectrale =", eta, "bit/s/Hz")
    print("Question 1.u/ \nMarge bande passante =", m, "Hz")
    print("Question 1.v/ \nDébit maximum =", deb_max, "bit/s")
