import matplotlib.pyplot as plt
import numpy as np

x = [0.25, 0.25, 1.25, 0.5, 1, 0.25, 0.6, 0, -0.6, -0.25, -1, -0.5, -1.25, -0.25, -0.25, 0.25]
y = [0, 0.5, 0.5, 1, 1, 1.5, 1.5, 2, 1.5 , 1.5, 1, 1, 0.5, 0.5, 0, 0]

################ VARIABLES ############## (toutes les distances sont en kilomètres)
gamma = 1 #constante de RG
c = 2.99792458*10**8 #célérité de la lumière (en m/s)
G = 6.6738480*10**-11 #constante gravitationnelle (en m**3 kg**-1 s**-2)
masseAstre = 1.989*10**30 #masse du soleil par défault, peut être modifié par la suite (en kg)
rayonSchw = (2*G*masseAstre)/(c**2) #rayon de Schwarzschild (en mètres)
vitTerre = 30 #vitesse de la Terre sur son orbite (en km/s)
b = 1120000 #paramètre d'impact (distance la plus courte entre le centre de l'astre de le rayonnement qui l'approche) 
r1 = 125400000 #distance Terre - Astre (en km)
r2 = 116165132 #distance Astre - Sonde (en km)



############## CALCUL DECALAGE FREQUENTIEL ########
freqShift = -(rayonSchw/c)*(1+gamma)*(1/b)*vitTerre

print("Masse de l'astre : "+str(masseAstre))
print("Rayon de Schwarzschild : "+str(rayonSchw))
print("Décalage fréquentiel : "+str(freqShift))

#plt.plot(x, y, '-.', color = "green", lw = 2)
#plt.title("Mon beau sapin")
#plt.axis('equal')
#plt.xlabel("C'est Noel")
#plt.ylabel("Vive le vent")
#plt.show()
#plt.close()

