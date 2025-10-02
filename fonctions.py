import numpy as np


def D(Re, altitude_sat,Angle_elevation): #en m
    """ Distance satellite en m en se basant sur le rayon de la Terre,
    L'altitude de la station, l'altitude du satellite et l'angle d'élévation"""
    phi = 90-Angle_elevation-np.arcsin((Re/(Re+altitude_sat))*np.cos(Angle_elevation))
    return np.sqrt(Re**2+(Re+altitude_sat)**2-2*Re*(Re+altitude_sat)*np.cos(phi))

def lamb(c,f):
    """ Longueur d'onde en fonction de la fréquence et de la vitesse de la lumière """
    return c/f

def Gain(D, lamb): # en dB
    """ Gain bord en dB """
    Ae =  0.56*np.pi*D**2/4
    return 10*np.log(4*np.pi*Ae/lamb**2)

def PIRE(P,G): # en dB
    """ Puissance isotrope rayonnée équivalente en dB """
    return P+G

def LFS (D,lamb):
    """ Pertes en espace libre en dB"""
    return 20*np.log(lamb/(4*np.pi*D))

def att_pluie (kH, kV, aH, aV, E, polar, R):
    """ angle de polartisation de 45° donc polarisation circulaire 
    consistant avec la polarisation utilisée pour l'espace"""
    k_coef = (kH + kV + (kH - kV)*(np.cos(E))**2*np.cos(2*polar))/2
    alpha_coef = (kV * aV + kH * aH + (kH*aH-kV*aV)*(np.cos(E))**2*np.cos(2*polar))/2*k_coef(E, polar, kH, kV)
    return k_coef(E, polar, kH, kV)*R**(alpha_coef(kH, kV, aH, aV, E, polar))

def perte_polar (Ar_stat, Ar_sat,polar):
    """ Perte de polarisations où Ar_stat est le ratio axial du récepteur,
    Ar_sat est le ratio axial de l'émetteur et polar, l'angle de polarisation"""
    return 0.5 + (4*Ar_stat*Ar_sat + (Ar_stat**2-1)*(Ar_sat**2-1)*np.cos(2*polar))/(2*(Ar_stat**2+1)*(Ar_sat**2+1))

def perte_depoint ():
    return 