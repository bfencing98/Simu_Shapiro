###############################################################
##                                                           ##
## ISEP 2020 - Sciences du numérique - Sujet : Effet Shapiro ##
##       Auteurs : Baptiste BONNAUDET & Martin LE GOFF       ##
##                                                           ##
###############################################################

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style

import numpy as np

import tkinter as tk
from tkinter import ttk

############################### VARIABLES ############################## (toutes les distances sont en kilomètres)
gamma = 1 #constante de RG
c = 2.99792458*10**8 #célérité de la lumière (en m/s)
G = 6.6738480*10**-11 #constante gravitationnelle (en m**3 kg**-1 s**-2)
masseSoleil = 1.989*10**30 #masse du soleil par défault, peut être modifié par la suite (en kg)
masseAstre = masseSoleil #masse du soleil par défault, peut être modifié par la suite (en kg)
rayonSchw = (2*G*masseAstre)/(c**2) #rayon de Schwarzschild (en mètres)
vitTerre = 30 #vitesse de la Terre sur son orbite (en km/s)
vitSonde = 47000 #vitesse de la sonde Cassini à partir du 30 décembre 2000 (en m/s)
dureeObs = 30 #durée de la période d'observation (en jours) 
freqMesure = 1 #nombre de mesures prises par jour d'observation
distTerreAstre = 149597870 #distance Terre - Astre (en km)
coorSonde = [[0,0]]*dureeObs*freqMesure #tableau de coordonnées de points appartenant à la trajectoire de la sonde
coorAstre = [0,0] #coordonnées de l'astre dans un repère orthonormé dont il est l'origine 
coorTerre = [0,-distTerreAstre] #coordonnées de la Terre si l'on considère un repère orthonormé 
coeffTrajSonde = 3 #coefficient de la trajectoire de la sonde
deltaPosSonde = ((3600*24)/freqMesure)*vitSonde/1000 #distance parcourue entre deux mesures par la sonde (en km)
coorSondeInit = ((0+((int(len(coorSonde))/2)*np.sqrt((deltaPosSonde**2)/coeffTrajSonde))),(1114402.130+((int(len(coorSonde))/2)*np.sqrt((deltaPosSonde**2)-((deltaPosSonde**2)/coeffTrajSonde))))) #coordonnées initiales de la sonde au début de la période d'observations
retardSh = [0]*len(coorSonde)
freqShift = [0]*len(coorSonde)
xAxe = [0]*len(coorSonde)
bMin = 0
LARGE_FONT= ("Verdana", 15)
NORMAL_TEXT= ("Verdana", 10)
style.use('ggplot')

def calculFS():
    global rayonSchw
    global freqShift
    global retardSh
    global xAxe
    ################# RECALCUL DES PARAMETRES ############################
    rayonSchw = (2*G*masseAstre)/(c**2) #rayon de Schwarzschild (en mètres)
    coorSonde = [[0,0]]*dureeObs*freqMesure #tableau de coordonnées de points appartenant à la trajectoire de la sonde
    deltaPosSonde = ((3600*24)/freqMesure)*vitSonde/1000 #distance parcourue entre deux mesures par la sonde (en km)
    coorSondeInit = ((0+((int(len(coorSonde))/2)*np.sqrt((deltaPosSonde**2)/coeffTrajSonde))),(1114402.130+((int(len(coorSonde))/2)*np.sqrt((deltaPosSonde**2)-((deltaPosSonde**2)/coeffTrajSonde))))) #coordonnées initiales de la sonde au début de la période d'observations
    retardSh = [0]*len(coorSonde)
    freqShift = [0]*len(coorSonde)
    xAxe = [0]*len(coorSonde)
       
    ################ CALCUL COORDONNEES SONDE ############################
    for i in range (0,(int(len(coorSonde)))):
        newCoorSonde = [0,0]
        #if (i>= 0):
        newCoorSonde[0] = coorSondeInit[0]-(i*np.sqrt((deltaPosSonde**2)/coeffTrajSonde))
        newCoorSonde[1] = coorSondeInit[1]+(i*np.sqrt((deltaPosSonde**2)-((deltaPosSonde**2)/coeffTrajSonde)))
        #else:
        #    newCoorSonde[0] = coorSondeInit[0]-(i*np.sqrt((deltaPosSonde**2)/coeffTrajSonde))
        #    newCoorSonde[1] = coorSondeInit[1]+(i*np.sqrt((deltaPosSonde**2)-((deltaPosSonde**2)/coeffTrajSonde)))
        coorSonde[i] = newCoorSonde
        xAxe[i] = (i-(int(len(coorSonde)/2)))*(dureeObs/(dureeObs*freqMesure))
        #print("Coordonnées point de mesure "+str(i)+" : ("+str(newCoorSonde[0])+";"+str(newCoorSonde[1])+")")


    for x in range (0,len(coorSonde)):#

        ########################## CALCUL PARAMETRES #####################
        r1 = np.sqrt((coorAstre[0]-coorTerre[0])**2+(coorAstre[1]-coorTerre[1])**2) #distance Terre - Astre (en km)
        r2 = np.sqrt((coorSonde[x][0]-coorAstre[0])**2+(coorSonde[x][1]-coorAstre[1])**2) #distance Astre - Sonde (en km)
        vectDir = (coorSonde[x][0]-coorTerre[0],coorSonde[x][1]-coorTerre[1])
        coeffA =  vectDir[1]
        coeffB = -1*vectDir[0]
        coeffC = coorTerre[0]*vectDir[1]+coorTerre[1]*vectDir[0]
        b = -1*((coeffA*coorAstre[0])+(coeffB*coorAstre[1])+(coeffC))/(np.sqrt((coeffA**2)+(coeffB**2)))


        ######################## CALCUL RETARD SHAPIRO ###################
        #d'après la formule pour un aller-retour à proximité d'un astre massif
        if (b!=0):
            retardSh[x] = (1+gamma)*(rayonSchw/c)*(np.log((4*r1*r2)/(b**2))+1)
        else:
            retardSh[0] = 99999999999

        ##################### CALCUL DECALAGE FREQUENTIEL ####################
        if (b!=0):
            freqShift[x] = -(rayonSchw/c)*(1+gamma)*(1/b)*vitTerre
        else:
            retardSh[0] = 99999999999

        if freqShift[x] > freqShift[x-1] :
            bMin = b


    ##################### AFFICHAGE DES RESULTATS ########################
    print("Masse de l'astre : "+str(masseAstre))
    print("Rayon de Schwarzschild : "+str(rayonSchw))
    #print("Retard Shapiro : "+str(retardSh))
    #print("Décalage fréquentiel : "+str(freqShift))
    print("Distance entre mesures : "+str(deltaPosSonde))
    print("Coordonnées point de mesure 1 : ("+str(coorSonde[0][0])+","+str(coorSonde[0][1])+")")
    print("Coordonnées point de mesure 2 : ("+str(coorSonde[1][0])+","+str(coorSonde[1][1])+")")
    print("Distance parcourue par la sonde au cours de la période d'observation : "+str(dureeObs*24*3600*vitSonde/1000)+" km")
    print("Distance entre le premier et le dernier point de mesure : "+str(np.sqrt(((coorSonde[int(len(coorSonde))-1][0]-coorSonde[0][0])**2+(coorSonde[int(len(coorSonde))-1][1]-coorSonde[0][1])**2)))+" km")
    print("paramètre bMin : "+str(bMin))
    print(len(coorSonde))

    #plt.plot(xAxe, freqShift, '-', color = "blue", lw = 1)
    #plt.title("Décalage fréquentiel subi aux alentours d'un astre massif\n")
    #plt.xlabel("Jours d'observations")
    #plt.ylabel("Décalage fréquentiel")
    #plt.show()
    #plt.close()


#calculFS()

f = Figure(figsize=(6,5), dpi=100)
a = f.add_subplot(111)

def animate(i):
    a.clear()
    a.plot(xAxe, freqShift)
    #a.set_title("Décalage fréquentiel au cours du temps d'observation")
    a.set_xlabel("Jours autour de la conjonction solaire")
    a.set_ylabel("Décalage fréquentiel")


class SimuShapiroApp(tk.Tk):

    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)
        tk.Tk.iconbitmap(self,default='simu-shapiro.ico')
        tk.Tk.wm_title(self, "Simulation - Effet Shapiro")
        
        container = tk.Frame(self)

        container.pack(side="top", fill="both", expand = True)

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        for F in (StartPage, Courbe):

            frame = F(container, self)

            self.frames[F] = frame

            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)

    def show_frame(self, cont):

        frame = self.frames[cont]
        frame.tkraise()

        
class StartPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        label = ttk.Label(self, text="Simulation du décalage fréquentiel induit par l'effet Shapiro", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        Text1= ttk.Label(self, text = "Vous pouvez modifier différents paramètres de la simulation à l'aide des curseurs ci-dessous.\n\
Les valeurs affichées par défaut sont celles correspondant à l'observation du retard Shapiro \n\
par la sonde Cassini-Huygens en juin 2002.\n\
 Après avoir modifié (ou pas) les paramètres, cliquez sur le bouton \"Générer la courbe\" pour afficher la courbe qu'aurait produit une observation  \n\
avec les paramètres spécifiés.\n", justify="center", font=NORMAL_TEXT) #Un label pour afficher du texte
        Text1.pack() # insère les Widgets dans la fenêtre

        CurseurDureeObs= tk.Scale(self, orient='horizontal', from_=10, to=60, resolution=1, tickinterval=5, length=500, label="Durée de la période d'observation", font= NORMAL_TEXT) #On définit l'objet Entry (zone de saisie) qui porte le nom CurseurDureeObs
        CurseurDureeObs.set(dureeObs)
        CurseurDureeObs.pack()

        CurseurFreqMesure= tk.Scale(self, orient='horizontal', from_=1, to=96, resolution=1, tickinterval=10, length=500, label="Fréquence des mesures (par jour)", font= NORMAL_TEXT) #On définit l'objet Entry (zone de saisie) qui porte le nom CurseurDureeObs
        CurseurFreqMesure.set(freqMesure)
        CurseurFreqMesure.pack()

        CurseurMasseAstre= tk.Scale(self, orient='horizontal', from_=0, to=250, resolution=0.5, tickinterval=50, length=500, label="Masse de l'astre (masse du Soleil = 1)", font= NORMAL_TEXT) #On définit l'objet Entry (zone de saisie) qui porte le nom CurseurDureeObs
        CurseurMasseAstre.set(1)
        CurseurMasseAstre.pack()

        Space = ttk.Label(self, text = "\n\n")
        Space.pack()

        def updateVar():
            global dureeObs
            global freqMesure
            global masseAstre
            dureeObs = CurseurDureeObs.get()
            freqMesure = CurseurFreqMesure.get()
            masseAstre = CurseurMasseAstre.get()*masseSoleil

        button = ttk.Button(self, text="Générer la courbe", 
                            command=lambda: [updateVar(),calculFS(), controller.show_frame(Courbe),Courbe.refresh])
        button.pack()


class Courbe(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="Simulation du décalage fréquentiel induit par l'effet Shapiro", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        buttonHome = ttk.Button(self, text="Retour aux paramètres",
                            command=lambda: controller.show_frame(StartPage))
        buttonHome.pack()

        Space = ttk.Label(self, text = "\n")
        Space.pack()

        self.textInfos = tk.StringVar()
        self.textInfos.set(" La courbe a été obtenue avec les paramètres suivants :\n\n\
        Nombre de jours d\'observation : "+str(dureeObs)+" jours\n\
        Fréquence des mesures : "+str(freqMesure)+" par jour\n\
        Vitesse de la Terre : "+str(vitTerre)+" km/s\n\
        Vitesse de la sonde : "+str(vitSonde/1000)+" km/s\n\
        Distance Terre-Astre : "+str(np.abs(coorTerre[1]))+" km\n\
        Masse de l'astre approché : "+str(masseAstre)+" kg\n\
        Rayon de Schwarzschild de l'astre : "+str(rayonSchw)+" mètres\n\
        ")

        self.labelInfos = ttk.Label(self, textvariable=self.textInfos, font=NORMAL_TEXT)
        self.labelInfos.update()
        self.labelInfos.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        a.plot(xAxe, freqShift)
        #a.set_title("Décalage fréquentiel au cours du temps d'observation")
        a.set_xlabel("Jours d'observations")
        a.set_ylabel("Décalage fréquentiel")
        canvas = FigureCanvasTkAgg(f, self)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.refresh()

    def refresh(self):
        #print("refresh called")
        self.textInfos.set("La courbe a été obtenue avec les paramètres suivants :\n\n\
        Nombre de jours d\'observation : "+str(dureeObs)+" jours\n\
        Fréquence des mesures : "+str(freqMesure)+" par jour\n\
        Vitesse de la Terre : "+str(vitTerre)+" km/s\n\
        Vitesse de la sonde : "+str(vitSonde/1000)+" km/s\n\
        Distance Terre-Astre : "+str(np.abs(coorTerre[1]))+" km\n\
        Masse de l'astre approché : "+str(masseAstre)+" kg\n\
        Rayon de Schwarzschild de l'astre : "+str(rayonSchw)+" mètres\n\
        ")
        self.after(1000, self.refresh)
        #refresh_button = tk.Button(self, text="Refresh", command=self.refresh)
        #refresh_button.pack()

app = SimuShapiroApp()
ani = animation.FuncAnimation(f, animate, interval=1000)
app.mainloop()