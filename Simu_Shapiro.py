import matplotlib.pyplot as plt
import numpy as np

x = [0.25, 0.25, 1.25, 0.5, 1, 0.25, 0.6, 0, -0.6, -0.25, -1, -0.5, -1.25, -0.25, -0.25, 0.25]
y = [0, 0.5, 0.5, 1, 1, 1.5, 1.5, 2, 1.5 , 1.5, 1, 1, 0.5, 0.5, 0, 0]

################ VARIABLES ############## (toutes les distances sont en kilomètres)
gamma = 1 #constante de RG
c = 300000 #célérité de la lumière
masseAstre = 54168468
rayonAstre = 5165498
rayonSchw = 3 #rayon de Schwarzschild
vitTerre = 30 #vitesse de la Terre sur son orbite (en km/s)
b = 1120000 #paramètre d'impact
r1 = 125400000 #distance Terre - Astre
r2 = 116165132 #distance Astre - Sonde



############## CALCUL DECALAGE FREQUENTIEL ########
freqShift = -(rayonSchw/c)*(1+gamma)*(1/b)*vitTerre

print(freqShift)

#plt.plot(x, y, '-.', color = "green", lw = 2)
#plt.title("Mon beau sapin")
#plt.axis('equal')
#plt.xlabel("C'est Noel")
#plt.ylabel("Vive le vent")
#plt.show()
#plt.close()

