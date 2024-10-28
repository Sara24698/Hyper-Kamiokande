# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 11:32:46 2024

@author: ICTEA
"""

import numpy as np
import matplotlib.pyplot as plt
import wire
import biotsavart
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import scipy as scipy
from scipy import optimize


#Definición de parámetros

I_rectangular = 1
I_circulares = 1


Altura = 72
pos_espira_rectangular = np.arange(-32,33,2)
radio_cilindro = 34
pos_espira_circular = np.arange(-35.5,37,2)

radio_PMT_tapas = 31.95
radio_PMT= 32.4

Angulo = np.arange(0, 6.315, 0.0218*2)
z=np.arange(-32.522, 32.523, 0.707)
limite = 100*np.ones(len(Angulo))
Angulo_tapas=[]

Btotal=[]
Media = []
points2=[]
Puntos = []
Coordenadas = []
Bperp_paredes=[]
Optimizacion=[]
Longitud=[]
Bx=[]
By=[]
Bz=[]




x = np.arange(-31.815, 32, 0.707*2)
y = np.arange(-31.815, 32, 0.707)



#Posición de los PMTs

for r in range(len(x)):
	for h in range(len(y)):
		points2.append([x[r], y[h]])


for g in range(len(points2)):
	distancia = np.sqrt(points2[g][0]**2+points2[g][1]**2)
	if distancia <=32.01:
		Puntos.append([points2[g][0], points2[g][1], 32.9])
		Puntos.append([points2[g][0], points2[g][1], -32.9])


for i in range(len(z)):
	for j in range(len(Angulo)):
		Puntos.append([radio_PMT*np.cos(Angulo[j]), radio_PMT*np.sin(Angulo[j]), z[i]])



A = np.zeros((len(Puntos)*3,len(pos_espira_rectangular)+len(pos_espira_circular)))

#Matriz campo geomagneico

b = []


#df = pd.read_csv("Matriz_A.csv",  sep=',',  comment='#')

#Matriz componentes

for k in range(0, len(Puntos)):
	print(k)
	Bperp=[]
	for i in range(len(pos_espira_rectangular)):
		lado_espira_rect = np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)
		w1r = wire.Wire(path=wire.Wire.RectangularPath(dx=lado_espira_rect*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[i],0]).Translate([0,0, 0.5])
		sol = biotsavart.BiotSavart(wire=w1r)
		B= sol.CalculateB(points=[Puntos[k]])*(10**7)
        
        if Puntos[k][2] == 32.9 or Puntos[k][2] == -32.9: 
            B_perp.append(np.sqrt(B[0][0]**2+(B[0][1]+303)**2))
            b.append(np.sqrt(303**2 + 366**2))
        
        else:
            B_perp.append(np.sqrt((B[0][0]*(Puntos[k][1]/radio_PMT)-(B[0][1]+303)*(Puntos[k][0]/radio_PMT))**2+(B[0][2]-366)**2))
            b.append(np.sqrt((303*(Puntos[k][0]/radio_PMT))**2+(366)**2))

	for j in range(len(pos_espira_circular)):
		w1c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[j]])
		sol = biotsavart.BiotSavart(wire=w1c)
		B= sol.CalculateB(points=[Puntos[k]])*(10**7)
		if Puntos[k][2] == 32.9 or Puntos[k][2] == -32.9: 
            B_perp.append(np.sqrt(B[0][0]**2+(B[0][1]+303)**2))
        
        else:
            B_perp.append(np.sqrt((B[0][0]*(Puntos[k][1]/radio_PMT)-(B[0][1]+303)*(Puntos[k][0]/radio_PMT))**2+(B[0][2]-366)**2))


	A[k,:] = B_perp


b= np.array(b)
df = pd.DataFrame(A)
df.to_csv('./Matriz_A_perp', index = False)

print(scipy.optimize.lsq_linear(A, b).x)