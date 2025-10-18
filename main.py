import numpy as np
import constantes as ct
import matplotlib.pyplot as plt
from fonctions import *


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
    P_reçue = Pre(ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Depointage, ct.P_test, ct.Gain_bord[9], ct.kH, ct.kV, ct.aH, ct.aV, ct.E_test, ct.AR_stat, ct.AR_sat[9], ct.Angle_ellipse, ct.I_precip, ct.h_stat, ct.Lat, Loss_el)
    Trxin = temp_entree_recep(ct.Tlna, ct.Tmx, ct.Glna, ct.Tif, ct.Gmx)
    Tain = temperature_bruit_entree_antenne(ct.Tsky, ct.Tgd, ct.Tm, ct.Perte_atm)
    Taout = temperature_bruit_sortie_antenne(ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test)
    Trx = temperature_eq_recep(ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx)
    dens_spec = n0(ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx)
    CN0_ratio = C_N0_ratio(ct.Depointage, ct.Gain_bord[9], ct.kH, ct.kV, ct.aH, ct.aV, ct.E_test, ct.AR_stat, ct.AR_sat[9], ct.Angle_ellipse, ct.I_precip, ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx, ct.h_stat, ct.Lat, Loss_el)
    EbN0_ratio = Eb_N0_ratio(ct.B_test, ct.Depointage, ct.Gain_bord[9], ct.kH, ct.kV, ct.aH, ct.aV, ct.E_test, ct.AR_stat, ct.AR_sat[9], ct.Angle_ellipse, ct.I_precip, ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx, ct.h_stat, ct.Lat, Loss_el)
    B = bande_passante(ct.Roll_off_factor, ct.B_test)
    eta = eff_spectrale(ct.Code_correcteur, ct.B_test, ct.Roll_off_factor)
    m = marge(ct.Rapport_EB_N0, ct.B_test, ct.Depointage, ct.Gain_bord[9], ct.kH, ct.kV, ct.aH, ct.aV, ct.E_test, ct.AR_stat, ct.AR_sat[9], ct.Angle_ellipse, ct.I_precip, ct.E_test, ct.Re, ct.h_sat, ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx, ct.h_stat, ct.Lat, Loss_el)
    deb_max = debit_max(P_reçue, ct.Rapport_EB_N0, ct.Marge, dens_spec)
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
    print("Question 1.p/ \nDensité de puissance de bruit spectrale =", dens_spec, "W/Hz")
    print("Question 1.q/ \nRapport C/N0 =", CN0_ratio, "dBHz")
    print("Question 1.r/ \nRapport Eb/N0 =", EbN0_ratio, "dB")
    print("Question 1.s/ \nBande passante du signal =", B, "Hz")
    print("Question 1.t/ \nEfficacité spectrale =", eta)
    print("Question 1.u/ \nMarge =", m, "dB")
    print("Question 1.v/ \nDébit maximum =", deb_max/10**6, "Mbit/s")

    trace_courbes_parametriques(ct.Rapport_EB_N0, ct.Depointage, ct.Gain_bord[9], ct.kH, ct.kV, ct.aH, ct.aV, ct.Angle_elevation, ct.E_test, ct.AR_stat, ct.AR_sat[9], ct.Angle_ellipse, ct.I_precip, ct.Re, ct.h_sat, ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.Ptx, ct.P_test ,ct.Tf,ct.Lfrx, ct.h_stat, ct.Lat, Loss_el, ct.Debit)

    m_15_60_10 = marge(ct.Rapport_EB_N0, 10e6, ct.Depointage, ct.Gain_bord[6], ct.kH, ct.kV, ct.aH, ct.aV, 60, ct.AR_stat, ct.AR_sat[6], ct.Angle_ellipse, ct.I_precip, 60, ct.Re, ct.h_sat, ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx, ct.h_stat, ct.Lat, Loss_el)
    m_15_60_2 = marge(ct.Rapport_EB_N0, 2e6, ct.Depointage, ct.Gain_bord[6], ct.kH, ct.kV, ct.aH, ct.aV, 60, ct.AR_stat, ct.AR_sat[6], ct.Angle_ellipse, ct.I_precip, 60, ct.Re, ct.h_sat, ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx, ct.h_stat, ct.Lat, Loss_el)
    
    print("Question 2.b/ \nmarge pour puissance de 1.5 W, angle d'élévation de 60° et débit (maximal) de 10Mbit/s =", m_15_60_10, "dB",
          "              \nmarge pour puissance de 1.5 W, angle d'élévation de 60° et débit (minimal) de 2Mbit/s =", m_15_60_2, "dB")

    m_025_60_10 = marge(ct.Rapport_EB_N0, 10e6, ct.Depointage, ct.Gain_bord[6], ct.kH, ct.kV, ct.aH, ct.aV, 60, ct.AR_stat, ct.AR_sat[6], ct.Angle_ellipse, ct.I_precip, 60, ct.Re, ct.h_sat, ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, 0.25, ct.Tf,ct.Lfrx, ct.h_stat, ct.Lat, Loss_el)
    m_025_60_2 = marge(ct.Rapport_EB_N0, 2e6, ct.Depointage, ct.Gain_bord[6], ct.kH, ct.kV, ct.aH, ct.aV, 60, ct.AR_stat, ct.AR_sat[6], ct.Angle_ellipse, ct.I_precip, 60, ct.Re, ct.h_sat, ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, 0.25, ct.Tf,ct.Lfrx, ct.h_stat, ct.Lat, Loss_el)
    
    print("Question 2.c/ \nmarge pour puissance de 0.25 W, angle d'élévation de 60° et débit (maximal) de 10Mbit/s =", m_025_60_10, "dB",
          "              \nmarge pour puissance de 0.25 W, angle d'élévation de 60° et débit (minimal) de 2Mbit/s =", m_025_60_2, "dB")
    
    m_QPSK = marge(9.59, 2e6, ct.Depointage, ct.Gain_bord[6], ct.kH, ct.kV, ct.aH, ct.aV, 60, ct.AR_stat, ct.AR_sat[6], ct.Angle_ellipse, ct.I_precip, 60, ct.Re, ct.h_sat, ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx, ct.h_stat, ct.Lat, Loss_el)
    m_8PSK = marge(12.97, 2e6, ct.Depointage, ct.Gain_bord[6], ct.kH, ct.kV, ct.aH, ct.aV, 60, ct.AR_stat, ct.AR_sat[6], ct.Angle_ellipse, ct.I_precip, 60, ct.Re, ct.h_sat, ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx, ct.h_stat, ct.Lat, Loss_el)
    m_16PSK = marge(17.43, 2e6, ct.Depointage, ct.Gain_bord[6], ct.kH, ct.kV, ct.aH, ct.aV, 60, ct.AR_stat, ct.AR_sat[6], ct.Angle_ellipse, ct.I_precip, 60, ct.Re, ct.h_sat, ct.k,ct.c, ct.f, ct.Diametre_antenne, ct.Efficacite_antenne, ct.Tsky, ct.Tgd, ct.Tm, ct.P_test,ct.Tf,ct.Lfrx, ct.h_stat, ct.Lat, Loss_el)

    print("Question 2.d/ \nmarge pour puissance de 1.5 W, angle d'élévation de 60° et débit 2 Mbit/s et une modulation en QPSK =", m_QPSK, "dB",
          "\nmarge pour puissance de 1.5 W, angle d'élévation de 60° et débit de 2 Mbit/s et une modulation en 8PSK =", m_8PSK, "dB",
          "\nmarge pour puissance de 1.5 W, angle d'élévation de 60° et débit de 2 Mbit/s et une modulation en 16PSK =", m_16PSK, "dB")
    
    eff_QPSK = eff_spectrale(ct.Code_correcteur, ct.B_test, ct.Roll_off_factor, 4)
    eff_8PSK = eff_spectrale(ct.Code_correcteur, ct.B_test, ct.Roll_off_factor)
    eff_16PSK = eff_spectrale(ct.Code_correcteur, ct.B_test, ct.Roll_off_factor, 16)

    print("Question 2.e/ \nefficacite pour puissance de 1.5 W, angle d'élévation de 60° et débit 2 Mbit/s et une modulation en QPSK =", eff_QPSK,
          "\nefficacite pour puissance de 1.5 W, angle d'élévation de 60° et débit de 2 Mbit/s et une modulation en 8PSK =", eff_8PSK, 
          "\nefficacite pour puissance de 1.5 W, angle d'élévation de 60° et débit de 2 Mbit/s et une modulation en 16PSK =", eff_16PSK)
