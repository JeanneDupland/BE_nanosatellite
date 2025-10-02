## Modulation 8-PSK
Rapport_EB = 12.97 # en dB
Debit = [i for i in range(2, 11)] # en Mbit/s
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
Axial_ration = 2
Angle_elevation = [j for j in range(5, 91)] # en degres
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
Intensite_precipitations = 85 # en mm/h
Perte_atmospherique = 0.2 # en dB
Tsky = 20 # en K
Tgd = 45 # en K
Tm = 275 # en K

## Constantes physiques
c = 2.99792458e8 # en m/s
Re = 6371e3 # en metres
k = 1.38066e-23 # en J/K

## Constantes diverses
Ptx = 1.5 # en W
E = 90 # en degres
B = 2 # en Mbit/s