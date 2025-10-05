import numpy as np
import constantes as ct

def D(Re, altitude_sat,Angle_elevation): #en m
    """ Distance satellite en m en se basant sur le rayon de la Terre,
    L'altitude de la station, l'altitude du satellite et l'angle d'élévation"""
    Dist =[]
    for i in Angle_elevation:
        phi = 90-i-np.arcsin((Re/(Re+altitude_sat))*np.cos(i))
        Dist.append(np.sqrt(Re**2+(Re+altitude_sat)**2-2*Re*(Re+altitude_sat)*np.cos(phi)))
    return np.array(Dist)

def lamb(c,f):
    """ Longueur d'onde en fonction de la fréquence et de la vitesse de la lumière """
    return c/f

def Gain(D, lamb): # en dB
    """ Gain bord en dB """
    G = []
    for i in D:
        Ae =  0.56*np.pi*i**2/4
        G.append(10*np.log10(4*np.pi*Ae/lamb**2))
    return np.array(G)

def PIRE(P,G): # en dB
    """ Puissance isotrope rayonnée équivalente en dB """
    EIRP = []
    for i in P:
        for j in G:
            EIRP.append(i+j)
    return np.array(EIRP)

def LFS (D,lamb):
    """ Pertes en espace libre en dB"""
    perte_es = []
    for i in D:
        perte_es.append(20*np.log10(lamb/(4*np.pi*i)))
    return np.array(perte_es)

def att_pluie (kH, kV, aH, aV, E, polar, R):
    """ angle de polartisation de 45° donc polarisation circulaire 
    consistant avec la polarisation utilisée pour l'espace"""
    pluie = []
    for i in E:
        k_coef = (kH + kV + (kH - kV)*(np.cos(i))**2*np.cos(2*polar))/2
        alpha_coef = (kV * aV + kH * aH + (kH*aH-kV*aV)*(np.cos(i))**2*np.cos(2*polar))/2*k_coef
        pluie.append(10*np.log10(k_coef*R**(alpha_coef)))
    return np.array(pluie)

def Ls(hr,hs,E): 
    "Pour E>5° longueur du trajet oblique en dessous de la hauteur de pluie"
    Ls = []
    for i in E:
        Ls.append((hr - hs)/np.sin(i))
    return Ls

def Lg(Ls,E): 
    "Projection horizontale Lg de la longueur du trajet oblique"
    Lg = []
    for i in E: 
        Lg.append(Ls*np.cos(i))
    return Lg

def fc (Lg,pluie,f): 
    "facteur de réduction horizontale pour 0,01% du temps"
    r001 = 1 / (1+0.78*np.sqrt(Lg*pluie/f) - 0.38*(1 - np.exp(-2*Lg)))
    return r001

def f_zetha (hr,hs,Lg, r001): 
    "zetha pour facteur de réduction"
    zetha = np.arctan((hr-hs)/Lg*r001)
    return zetha

def Lr(Lg,r001,E,zetha,hr,hs): 
    if zetha>E: 
        Lr = Lg*r001/np.cos(E)
    else: 
        Lr = (hr-hs)/np.sin(E)
    return Lr

" angle de polarisation de 45° donc Xhi=0"

def f_v001 (E,Lr,pluie,f): 
    "facteur d'ajustement v001 pour 0,01% du temps"
    v001 = 1/(1+np.sqrt(np.sin(E))*(31*(1-np.exp(-E))*np.sqrt(Lr*pluie)/f**2 - 0.45))
    return v001

def Le (Lr,v001): 
    "longueur effective du trajet"
    Le = Lr * v001
    return Le

def A001 (pluie,Le): 
    "Affaiblissement prévu dépassé pour 0,01% d'une année moyenne"
    A001= pluie*Le
    return A001

def perte_polar (Ar_stat, Ar_sat,polar):
    """ Perte de polarisations où Ar_stat est le ratio axial du récepteur,
    Ar_sat est le ratio axial de l'émetteur et polar, l'angle de polarisation"""
    pol = []
    for i in Ar_sat:
        pol.append(10*np.log10(0.5 + (4*Ar_stat*i + (Ar_stat**2-1)*(i**2-1)*np.cos(2*polar))/(2*(Ar_stat**2+1)*(i**2+1))))
    return np.array(pol)

def perte_depoint (lamb, diam_ant, E):
    depoint = []
    theta_3dB = 70*lamb/diam_ant
    for i in E:
        depoint.append(-12*(90/theta_3dB)**2)
    return np.array(depoint)

if __name__ == "__main__":
    P = ct.Ptx
    Dist = D(ct.Re, ct.h_sat, ct.Angle_elevation)
    lam = lamb(ct.c, ct.f)
    G = Gain(Dist, lam)
    EIRP = PIRE(P,G)
    Perte_EL = LFS(Dist, lam)
    Att_spe_pluie = att_pluie(ct.kH, ct.kV, ct.aH, ct.aV, ct.Angle_elevation, ct.Angle_ellipse, ct.I_precip)
    polar = perte_polar(ct.AR_stat, ct.AR_sat, ct.Angle_ellipse)
    depoint = perte_depoint(lam, ct.Diametre_antenne, ct.Angle_elevation)
    print('PIRE=', EIRP)
    print('Distance=', Dist)
    print('Pertes=', Perte_EL)
    print('Atténuatiion spé pluie=', Att_spe_pluie)
    print('pertes polarisation =', polar)
    print('pertes dépointage=', depoint)
    print("atténuation due à la pluie dépassée pendant 0,01%=",A001)
