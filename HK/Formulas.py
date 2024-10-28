# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 22:54:21 2022

@author: sara2
"""

import numpy as np
import matplotlib.pyplot as plt
import wire
import biotsavart

P =[0, 20, 0]

a_z = 36


pos_espira = np.arange(-32,34,2)


radio_cilindro = 34

Campo_Bz = []




for i in range(len(pos_espira)):

	y = P[1]-pos_espira[i]

	a_x = np.sqrt(radio_cilindro**2 - pos_espira[i]**2)


	r1 = np.sqrt((P[0]+a_x)**2+(P[2]+a_z)**2+y**2)
	r2 = np.sqrt((P[0]-a_x)**2+(P[2]+a_z)**2+y**2)
	r3 = np.sqrt((P[0]+a_x)**2+(P[2]-a_z)**2+y**2)
	r4 = np.sqrt((P[0]-a_x)**2+(P[2]-a_z)**2+y**2)

	Bx = -y*(1/(r1*(r1-P[2]-a_z))-1/(r2*(r2-P[2]-a_z))-1/(r3*(r3-P[2]+a_z))+1/(r4*(r4-P[2]+a_z)))

	By = (P[0]+a_x)/(r1*(r1-P[2]-a_z))+(P[2]+a_z)/(r1*(r1-P[0]-a_x))-(P[0]-a_x)/(r2*(r2-P[2]-a_z))-(P[2]+a_z)/(r2*(r2-P[0]+a_x))-(P[0]+a_x)/(r3*(r3-P[2]+a_z))-(P[2]-a_z)/(r3*(r3-P[0]-a_x))+(P[0]-a_x)/(r4*(r4-P[2]+a_z))+(P[2]-a_z)/(r4*(r4-P[0]+a_x))

	Bz = -y*(1/(r1*(r1-P[0]-a_x))-1/(r2*(r2-P[0]+a_x))-1/(r3*(r3-P[0]-a_x))+1/(r4*(r4-P[0]+a_x)))
	
	Campo_Bz.append(Bz)


print(np.sum(Campo_Bz))




#Campo_total = sum(B_perp)	
	

#print("El campo magnético perpendicular en el punto P es de", Campo_total, "mG")
#print(Campo_mag)





