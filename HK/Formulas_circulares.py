# -*- coding: utf-8 -*-
"""
Created on Thu May  5 11:42:21 2022

@author: ICTEA
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as symp
from sympy.abc import x
from sympy import cos, sin, pi, sqrt
from scipy import special as sp
from scipy import pi
import wire
import biotsavart
import Corrientes
import Corrientesanalíticas






#Definición de constantes

I_rectangular = 62.5#74.5111705813268#I[0]
I_circulares = 75.87308111813685#I[1]
I_topbottom = -75.87308111813685#I[1]

Altura = 72
pos_espira_rectangular = np.arange(-32,33,2)
radio_cilindro = 34
radio_tapas = 30
pos_espira_circular = np.array([-36, -36, -36, -34, -32, -30, -28, -26, -24, -22, -20, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 36, 36, 36, 36, 36])
PMTs_vertical = np.arange(-34,34,1)
limite = 100*np.ones(len(PMTs_vertical))
radio_PMT=32.4
z=np.arange(-33.4, 32.5, 0.707)

Angulo1 = 0
Angulo2 = np.pi/6
Angulo3 = np.pi/4
Angulo4 = np.pi/3
Angulo5 = np.pi/2

Angulo = np.arange(0, 6.2667, 0.02222)



points1=[]
points2=[]
points3=[]
points4=[]
points5=[]


Br_circular=[]
Bz_circular=[]


Bx1_circular=[]
By1_circular=[]

Bx2_circular=[]
By2_circular=[]

Bx3_circular=[]
By3_circular=[]

Bx4_circular=[]
By4_circular=[]

Bx5_circular=[]
By5_circular=[]

"""

for k in range(len(PMTs_vertical)):
     points1.append([radio_PMT*np.cos(Angulo1), radio_PMT*np.sin(Angulo1) , PMTs_vertical[k]])
     points2.append([radio_PMT*np.cos(Angulo2), radio_PMT*np.sin(Angulo2) , PMTs_vertical[k]])
     points3.append([radio_PMT*np.cos(Angulo3), radio_PMT*np.sin(Angulo3) , PMTs_vertical[k]])
     points4.append([radio_PMT*np.cos(Angulo4), radio_PMT*np.sin(Angulo4) , PMTs_vertical[k]])
     points5.append([radio_PMT*np.cos(Angulo5), radio_PMT*np.sin(Angulo5) , PMTs_vertical[k]])

"""

points = []
for j in range(len(Angulo)):
 points.append([radio_PMT*np.cos(Angulo[j]), radio_PMT*np.sin(Angulo[j]), -30.571])


z1=-30.571



#Espiras circulares


	
Br=[]
Bz=[]

for j in range(len(pos_espira_circular)):
 m = 4*radio_PMT*radio_cilindro/((radio_PMT+radio_cilindro)**2+(z1-pos_espira_circular[j])**2)  
 E = sp.ellipeinc(2*pi, m)
 K = sp.ellipkinc(2*pi, m)




 Factor_1Br = 2*(z1-pos_espira_circular[j])/(radio_PMT*np.sqrt((radio_cilindro+radio_PMT)**2+(z1-pos_espira_circular[j])**2))
 Factor_2Br = (radio_cilindro**2+radio_PMT**2+(z1-pos_espira_circular[j])**2)/((radio_cilindro-radio_PMT)**2+(z1-pos_espira_circular[j])**2)
 Br.append(Factor_1Br*(Factor_2Br*E-K))




 Factor_1Bz = 2/(np.sqrt((radio_cilindro+radio_PMT)**2+(z1-pos_espira_circular[j])**2))
 Factor_2Bz = (radio_cilindro**2-radio_PMT**2-(z1-pos_espira_circular)**2)/((radio_cilindro-radio_PMT)**2+(z1-pos_espira_circular[j])**2)

Bz.append(Factor_1Bz*(Factor_2Bz*E+K))



Br_suma_espiras = np.sum(Br)
Bz_suma_espiras = np.sum(Bz)


Br_circular.append(Br_suma_espiras)
Bz_circular.append(I_circulares*Bz_suma_espiras)

for g in range(len(Angulo)):
 Bx1_circular.append(I_circulares*Br_suma_espiras*np.cos(Angulo[g]))
 By1_circular.append(I_circulares*Br_suma_espiras*np.sin(Angulo[g]))


# rectangular loops I=1 A
    

w1r = wire.Wire(path=wire.Wire.RectangularPath(dx=np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[0], 0])
sol = biotsavart.BiotSavart(wire=w1r)

for i in range(1, len(pos_espira_rectangular)):
 lado_espira_rect = np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)
 w2r = wire.Wire(path=wire.Wire.RectangularPath(dx=lado_espira_rect*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[i],0])
 sol.AddWire(w2r)
 
# calculate B field at given points
B1 = sol.CalculateB(points=points)*(10**7)




#Arrays por componentes
    
Bx1=[]
By1=[]
Bz1=[]




for l in range(len(points)):
 Bx1.append(B1[l,0])
 By1.append(B1[l,1]+303)		
 Bz1.append(B1[l,2]-366)




Bx_total1 = np.add(Bx1, Bx1_circular)
By_total1 = np.add(By1, By1_circular)
Bz_total1 = np.add(Bz1, Bz_circular)
Bperp1=[]

for h in range(len(Angulo)):
 Bperp1.append(np.sqrt((Bx_total1[h]*np.sin(Angulo[h])-By_total1[h]*np.cos(Angulo[h]))**2+Bz_total1[h]**2))

print(Bperp1)
PMTs_malos = 0
Coordenadas_paredes=[]
    
for alfa in range(len(Bperp1)):
 if np.abs(Bperp1[alfa]) > 100:
  PMTs_malos = PMTs_malos+1
  Coordenadas_paredes.append([points[alfa]])

print('El número de PMTs malos es', PMTs_malos)
print(Coordenadas_paredes)
 




"""
#Figuras paredes


fig, ax = plt.subplots()
#ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
ax.scatter(PMTs_vertical, Bperp1, label='$\Delta B_{perp}$ $0 \: rad$')
ax.scatter(PMTs_vertical, Bperp2, label='$\Delta B_{perp}$ $\pi/6 \: rad$')
ax.scatter(PMTs_vertical, Bperp3, label='$\Delta B_{perp}$ $\pi/4 \: rad$')
ax.scatter(PMTs_vertical, Bperp4, label='$\Delta B_{perp}$ $\pi/3 \: rad$')
ax.scatter(PMTs_vertical, Bperp5, label='$\Delta B_{perp}$ $\pi/2 \: rad$')
ax.plot(PMTs_vertical, limite)
ax.plot(PMTs_vertical, -limite)
plt.xlabel('Height (m)')
plt.ylabel('B (mG)')
#plt.xlim(0,35)
plt.title('$\Delta B_{perp}$ at walls as a function of height and angle')
plt.legend()
plt.savefig('Bperp_walls.png')
plt.show()

fig, ax = plt.subplots()
#ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
ax.scatter(PMTs_vertical, Bz_total1, label='$\Delta B_{z}$ 0º')
ax.scatter(PMTs_vertical, Bz_total2, label='$\Delta B_{z}$ 30º')
ax.scatter(PMTs_vertical, Bz_total3, label='$\Delta B_{z}$ 45º')
ax.scatter(PMTs_vertical, Bz_total4, label='$\Delta B_{z}$ 60º')
ax.scatter(PMTs_vertical, Bz_total5, label='$\Delta B_{z}$ 90º')
ax.plot(PMTs_vertical, limite)
ax.plot(PMTs_vertical, -limite)
plt.xlabel('Height (m)')
plt.ylabel('B (mG)')
#plt.xlim(0,35)
plt.title('$\Delta B_{z}$ at walls as a function of height and angle')
plt.legend()
plt.savefig('Bz_walls.png')
plt.show()

fig, ax = plt.subplots()
#ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
ax.scatter(PMTs_vertical, Bx_total1, label='$B_{x}$ 0º')
ax.scatter(PMTs_vertical, Bx_total2, label='$B_{x}$ 30º')
ax.scatter(PMTs_vertical, Bx_total3, label='$B_{x}$ 45º')
ax.scatter(PMTs_vertical, Bx_total4, label='$B_{x}$ 60º')
ax.scatter(PMTs_vertical, Bx_total5, label='$B_{x}$ 90º')
ax.plot(PMTs_vertical, limite)
ax.plot(PMTs_vertical, -limite)
plt.xlabel('Height (m)')
plt.ylabel('B (mG)')
#plt.xlim(0,35)
plt.title('$B_{x}$ at walls as a function of height and angle')
plt.legend()
plt.savefig('Bx_walls.png')
plt.show()

fig, ax = plt.subplots()
#ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
ax.scatter(PMTs_vertical, By_total1, label='$\Delta B_{y}$ 0º')
ax.scatter(PMTs_vertical, By_total2, label='$\Delta B_{y}$ 30º')
ax.scatter(PMTs_vertical, By_total3, label='$\Delta B_{y}$ 45º')
ax.scatter(PMTs_vertical, By_total4, label='$\Delta B_{y}$ 60º')
ax.scatter(PMTs_vertical, By_total5, label='$\Delta B_{y}$ 90º')
ax.plot(PMTs_vertical, limite)
ax.plot(PMTs_vertical, -limite)
plt.xlabel('Height (m)')
plt.ylabel('B (mG)')
#plt.xlim(0,35)
plt.title('$\Delta B_{y}$ at walls as a function of height and angle')
plt.legend()
plt.savefig('By_walls.png')
plt.show()
"""
