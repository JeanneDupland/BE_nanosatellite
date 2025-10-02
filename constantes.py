import numpy as np

### Paramamètres nano satellite

## Paramètre bord
f = 8.2e9 # en Hz
h_sat = 420e3 # en m
Ptx = np.linspace(0.25, 1.5, 100)
Perte_bord = 0.5 # en dB

## Caractéristique antenne
Elevation_sol = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90]
Depointage_bord = [85, 68, 62, 54, 46, 37, 28, 19, 9, 0]
Gain_bord = [-7.5, -7.2, -5.7, -3.2, -2.3, -1.1, 0.8, 3.3, 6.1, 7.4]
AR_sat = [16.5, 16.5, 18.8, 20, 14.7, 7.2, 2.5, 5.8, 5.5, 4.5]


## Modulation 8-PSK
Rapport_EB_N0 = 12.97 # en dB
Debit = np.linspace(2,10,100) # en Mbit/s
Roll_off_factor = 0.25 
Code_correcteur = 0.5 
Marge = 3  # en dB
Taux_dispo = 99 # en %

## Station sol Kourou
Lat = 5.1 # en degres nord
Altitude = 300 # en metres
Type_antenne = 'Parabolique'
Diametre_antenne = 11.125 # en metres
Angle_ellipse = 45 # en degres
AR_stat = 2
Angle_elevation = np.linspace(5,90,100)# en degres
Depointage = 0.1 # en degres
Efficacite_antenne = 0.65
Tlna = 150 # en K
Tmx = 850 # en K
Tif = 400 # en K
Tf = 290 # en K
Glna = 50 # en dB
Gmx = -10 # en dB
Lfrx = -0.5 # en dB

## Canal de propagation Kourou
I_precip = 85 # en mm/h
Perte_atm = 0.2 # en dB
Tsky = 20 # en K
Tgd = 45 # en K
Tm = 275 # en K

## Constantes physiques
c = 2.99792458e8 # en m/s
Re = 6371e3 # en metres
k = 1.38066e-23 # en J/K

## Pertes pluie
