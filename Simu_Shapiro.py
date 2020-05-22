###############################################################
##                                                           ##
## ISEP 2020 - Sciences du numérique - Sujet : Effet Shapiro ##
##       Auteurs : Baptiste BONNAUDET & Martin LE GOFF       ##
##                                                           ##
###############################################################

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
vitSonde = 47000 #vitesse de la sonde Cassini à partir du 30 décembre 2000 (en m/s)
dureeObs = 30 #durée de la période d'observation (en jours) 
freqMesure = 24 #nombre de mesures prises par jour d'observation
distTerreAstre = 149597870 #distance Terre - Astre (en km)
coorSonde = [[0,0]]*dureeObs*freqMesure #tableau de coordonnées de points appartenant à la trajectoire de la sonde
coorAstre = [0,0] #coordonnées de l'astre dans un repère orthonormé dont il est l'origine 
coorTerre = [0,-distTerreAstre] #coordonnées de la Terre si l'on considère un repère orthonormé 
coorSondeInit = (50000000,850000000)#coordonnées initiales de la sonde au début de la période d'observations
deltaPointMesure = ((3600*24)/freqMesure)*vitSonde/1000 #distance parcourue entre deux mesures par la sonde (en km)



################ CALCUL COORDONNEES SONDE ############################
for i in range (0,dureeObs*freqMesure):
    newCoor = [0,0]
    newCoor[0] = coorSondeInit[0]-(i*np.sqrt((deltaPointMesure**2)/2))
    newCoor[1] = coorSondeInit[1]+(i*np.sqrt((deltaPointMesure**2)/2))
    coorSonde[i] = newCoor


######################## CALCUL RETARD SHAPIRO #######################
#d'après la formule pour un aller-retour à proximité d'un astre massif
retardSh = [0]*len(coorSonde)

for x in len(coorSonde):
    r1 = np.sqrt((coorAstre[0]-coorTerre[0])**2+(coorAstre[1]-coorTerre[0])**2) #distance Terre - Astre (en km)
    r2 = np.sqrt((coorSonde[x][0]-coorAstre[0])**2+(coorSonde[x][1]-coorAstre[1])**2) #distance Astre - Sonde (en km)

    b = 1120000 #paramètre d'impact (distance la plus courte entre le centre de l'astre de le rayonnement qui l'approche)

    retard[x] = (1+gamma)*(rayonSchw/c)*(np.log((4*r1*r2)/(b**2))+1)


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

