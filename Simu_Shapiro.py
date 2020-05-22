import matplotlib.pyplot as plt
import numpy as np

x = [0.25, 0.25, 1.25, 0.5, 1, 0.25, 0.6, 0, -0.6, -0.25, -1, -0.5, -1.25, -0.25, -0.25, 0.25]
y = [0, 0.5, 0.5, 1, 1, 1.5, 1.5, 2, 1.5 , 1.5, 1, 1, 0.5, 0.5, 0, 0]


############################### VARIABLES ############################## (toutes les distances sont en kilomètres)
gamma = 1 #constante de RG
c = 2.99792458*10**8 #célérité de la lumière (en m/s)
G = 6.6738480*10**-11 #constante gravitationnelle (en m**3 kg**-1 s**-2)
masseAstre = 1.989*10**30 #masse du soleil par défault, peut être modifié par la suite (en kg)
rayonSchw = (2*G*masseAstre)/(c**2) #rayon de Schwarzschild (en mètres)
vitTerre = 30 #vitesse de la Terre sur son orbite (en km/s)
b = 1120000 #paramètre d'impact (distance la plus courte entre le centre de l'astre de le rayonnement qui l'approche) 
r1 = 149597870 #distance Terre - Astre (en km)
r2 = 849000000 #distance Astre - Sonde (en km)
vitSonde = 47000 #vitesse de la sonde Cassini à partir du 30 décembre 2000 (en m/s)
dureeObs = 30 #durée de la période d'observation (en jours) 
coorSonde = [(0,0)]*720 #tableau de coordonnées de points appartenant à la trajectoire de la sonde
coorSondeInit = (50000000,850000000)#coordonnées initiales de la sonde au début de la période d'observations
freqMesure = 24 #nombre de mesures prises par jour d'observation
deltaPointMesure = ((3600*24)/freqMesure)*vitSonde/1000 #distance parcourue entre deux mesures par la sonde (en km)


################ CALCUL COORDONNEES SONDE ############################
for i in range (0,dureeObs*freqMesure):
    newCoor = [0,0]
    newCoor[0] = coorSondeInit[0]-(i*np.sqrt((deltaPointMesure**2)/2))
    newCoor[1] = coorSondeInit[1]+(i*np.sqrt((deltaPointMesure**2)/2))
    coorSonde[i] = newCoor


######################## CALCUL RETARD SHAPIRO #######################
retardSh = (1+gamma)*(rayonSchw/c)*(np.log((4*r1*r2)/(b**2))+1)


##################### CALCUL DECALAGE FREQUENTIEL ####################
freqShift = -(rayonSchw/c)*(1+gamma)*(1/b)*vitTerre


##################### AFFICHAGE DES RESULTATS ########################
print("Masse de l'astre : "+str(masseAstre))
print("Rayon de Schwarzschild : "+str(rayonSchw))
print("Retard Shapiro : "+str(retardSh))
print("Décalage fréquentiel : "+str(freqShift))
print("Distance entre mesures : "+str(deltaPointMesure))
print("Coordonnées point de mesure 1 : ("+str(coorSonde[0][0])+","+str(coorSonde[0][1])+")")
print("Coordonnées point de mesure 2 : ("+str(coorSonde[1][0])+","+str(coorSonde[1][1])+")")
print("Distance parcourue par la sonde au cours de la période d'observation : "+str(dureeObs*24*3600*vitSonde/1000)+" km")
print("Distance entre le premier et le dernier point de mesure : "+str(np.sqrt(((coorSonde[719][0]-coorSondeInit[0])**2+(coorSonde[719][1]-coorSondeInit[1])**2)))+" km")

#plt.plot(x, y, '-.', color = "green", lw = 2)
#plt.title("Mon beau sapin")
#plt.axis('equal')
#plt.xlabel("C'est Noel")
#plt.ylabel("Vive le vent")
#plt.show()
#plt.close()

