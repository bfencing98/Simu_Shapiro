###############################################################
##                                                           ##
## ISEP 2020 - Sciences du numérique - Sujet : Effet Shapiro ##
##       Auteurs : Baptiste BONNAUDET & Martin LE GOFF       ##
##                                                           ##
###############################################################

import matplotlib.pyplot as plt
import numpy as np


############################### VARIABLES ############################## 
gamma = 1 #constante de RG
c = 2.99792458*10**8 #célérité de la lumière (en m/s)
G = 6.6738480*10**-11 #constante gravitationnelle (en m**3 kg**-1 s**-2)
masseAstre = 1.989*10**30 #masse du soleil par défault, peut être modifié par la suite (en kg)
rayonSchw = (2*G*masseAstre)/(c**2) #rayon de Schwarzschild (en mètres)
vitTerre = 30000 #vitesse de la Terre sur son orbite (en m/s)
vitSonde = 47000 #vitesse de la sonde Cassini à partir du 30 décembre 2000 (en m/s)
dureeObs = 120 #durée de la période d'observation (en jours) 
freqMesure = 96 #nombre de mesures prises par jour d'observation
distTerreAstre = 149597870 #distance Terre - Astre (en km)
coorAstre = [0,0] #coordonnées de l'astre dans un repère orthonormé dont il est l'origine 
coorTerre = [[0,0]]*dureeObs*freqMesure #tableau de coordonnées de points appartenant à la trajectoire de la Terre
coorTerreInit = [-100000000,-111263303] #coordonnées de la Terre si l'on considère un repère orthonormé 
coorSonde = [[0,0]]*dureeObs*freqMesure #tableau de coordonnées de points appartenant à la trajectoire de la sonde
coorSondeInit = (500000000,827616305) #coordonnées initiales de la sonde au début de la période d'observations
coeffTrajSonde = 4.4 #coeffeicient de la trajectoire de la sonde
deltaPosSonde = ((3600*24)/freqMesure)*vitSonde/1000 #distance parcourue entre deux mesures par la sonde (en m)
coeffTrajTerre = 1.43 #coeffeicient de la trajectoire de la Terre
deltaPosTerre = ((3600*24)/freqMesure)*vitTerre/1000 #distance parcourue entre deux mesures par la terre (en m)
retardSh = [0]*len(coorSonde)
freqShift = [0]*len(coorSonde)
xAxe = [0]*len(coorSonde)


################ CALCUL COORDONNEES SONDE ############################
for i in range (0,len(coorSonde)):
    newCoorSonde = [0,0]
    newCoorSonde[0] = coorSondeInit[0]-(i*np.sqrt((deltaPosSonde**2)/(coeffTrajSonde**2)))
    newCoorSonde[1] = coorSondeInit[1]+(i*np.sqrt((deltaPosSonde**2)-((deltaPosSonde**2)/(coeffTrajSonde**2))))
    coorSonde[i] = newCoorSonde

    xAxe[i] = i*(dureeObs/(dureeObs*freqMesure))

    newCoorTerre = [0,0]
    newCoorTerre[0] = coorTerreInit[0]+(i*np.sqrt((deltaPosTerre**2)-((deltaPosTerre**2)/(coeffTrajTerre**2))))
    newCoorTerre[1] = coorTerreInit[1]-(i*np.sqrt((deltaPosTerre**2)/(coeffTrajTerre**2)))
    coorTerre[i] = newCoorTerre


for x in range (0,len(coorSonde)):#

    ########################## CALCUL PARAMETRES #####################
    r1 = np.sqrt((coorAstre[0]-coorTerre[x][0])**2+(coorAstre[1]-coorTerre[x][1])**2) #distance Terre - Astre (en km)
    r2 = np.sqrt((coorSonde[x][0]-coorAstre[0])**2+(coorSonde[x][1]-coorAstre[1])**2) #distance Astre - Sonde (en km)
    vectDir = (coorSonde[x][0]-coorTerre[x][0],coorSonde[x][1]-coorTerre[x][1])
    #print(vectDir)
    coeffA =  vectDir[1]
    #print(coeffA)
    coeffB = -1*vectDir[0]
    #print(coeffB)
    coeffC = coorTerre[x][0]*vectDir[1]+coorTerre[x][1]*vectDir[0]
    #print(coeffC)
    b = -1*((coeffA*coorAstre[0])+(coeffB*coorAstre[1])+(coeffC))/(np.sqrt((coeffA**2)+(coeffB**2))) # 
    #print(b)
    bbis = 1120000 #paramètre d'impact (distance la plus courte entre le centre de l'astre de le rayonnement qui l'approche)

    if x == 0 :
        b0 = b

    ######################## CALCUL RETARD SHAPIRO ###################
    #d'après la formule pour un aller-retour à proximité d'un astre massif
    retardSh[x] = (1+gamma)*(rayonSchw/c)*(np.log((4*r1*r2)/(b**2))+1)


    ##################### CALCUL DECALAGE FREQUENTIEL ####################
    freqShift[x] = -(rayonSchw/c)*(1+gamma)*(1/b)*(vitTerre/1000)

    if np.abs(freqShift[x]) > np.abs(freqShift[x-1]) :
        bMin = b


##################### AFFICHAGE DES RESULTATS ########################
print("Masse de l'astre : "+str(masseAstre))
print("Rayon de Schwarzschild : "+str(rayonSchw))
print("Retard Shapiro : "+str(retardSh[0]*1000000)+"ms avec r1 : "+str(np.sqrt((coorAstre[0]-coorTerre[0][0])**2+(coorAstre[1]-coorTerre[0][1])**2))+" r2 : "+str(np.sqrt((coorSonde[0][0]-coorAstre[0])**2+(coorSonde[0][1]-coorAstre[1])**2))+" b : "+str(b0))
#print("Décalage fréquentiel : "+str(freqShift))
print("Coordonnées point de mesure 1 : ("+str(coorSonde[0][0])+","+str(coorSonde[0][1])+")")
print("Coordonnées point de mesure 2 : ("+str(coorSonde[1][0])+","+str(coorSonde[1][1])+")")
print("Distance parcourue par la sonde au cours de la période d'observation : "+str(np.sqrt((coorSonde[0][0]-coorSonde[len(coorSonde)-1][0])**2+(coorSonde[0][1]-coorSonde[len(coorSonde)-1][1])**2))+" km")
print("Distance parcourue par la Terre : "+str(np.sqrt((coorTerre[0][0]-coorTerre[len(coorTerre)-1][0])**2+(coorTerre[0][1]-coorTerre[len(coorTerre)-1][1])**2))+" km")
print("paramètre bMin : "+str(bMin))
print("Coordonnées point de mesure 0 : ("+str(coorSonde[0][0])+","+str(coorSonde[0][1])+")")
print("Coordonnées point de mesure 1 : ("+str(coorSonde[1][0])+","+str(coorSonde[1][1])+")")
print("Distance entre mesure 0 et 1 : ("+str(np.sqrt((coorSonde[0][0]-coorSonde[1][0])**2+(coorSonde[0][1]-coorSonde[1][1])**2))+" km")
print("Coordonnées du dernier point de mesure : "+str(coorSonde[len(coorSonde)-1][0])+","+str(coorSonde[len(coorSonde)-1][1])+")")
print(deltaPosTerre)
print(deltaPosSonde)

plt.plot(xAxe, freqShift, '-', color = "blue", lw = 1)
plt.title("Décalage fréquentiel subi aux alentours d'un astre massif\n")
plt.xlabel("Jours d'observations")
plt.ylabel("Décalage fréquentiel")
plt.show()
plt.close()

