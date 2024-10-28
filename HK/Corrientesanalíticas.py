# -*- coding: utf-8 -*-
"""
Created on Thu May  5 11:42:21 2022

@author: ICTEA
"""

def corrientescirculares():

 import matplotlib.pyplot as plt
 import numpy as np
 import sympy as symp
 from sympy.abc import x
 from sympy import cos, sin, pi, sqrt
 from scipy import special as sp
 from scipy import pi
 import wire
 import biotsavart
 import Corrientes





 #Definición de constantes



 Altura = 72
 pos_espira_rectangular = np.arange(-32,33,2)
 radio_cilindro = 34
 radio_tapas = 30
 pos_espira_circular = np.arange(-36,37,2)
 PMTs_vertical = np.arange(-34,34,1)
 limite = 100*np.ones(len(PMTs_vertical))
 radio_PMT=32.4

 Angulo1 = 0
 Angulo2 = np.pi/6
 Angulo3 = np.pi/4
 Angulo4 = np.pi/3
 Angulo5 = np.pi/2



 points1=[]




 Br_circular=[]
 Bz_circular=[]




 for v in range(len(PMTs_vertical)):
  points1.append([0, 0 , PMTs_vertical[v]])

 z1 =[]



 for l in range(len(points1)):
  z1.append(points1[l][2])




 #Espiras circulares


 for w in range(len(z1)):

  Br=[]
  Bz=[]

  for j in range(len(pos_espira_circular)):
   m = 4*radio_PMT*radio_cilindro/((radio_PMT+radio_cilindro)**2+(z1[w]-pos_espira_circular[j])**2)  
   E = sp.ellipeinc(2*pi, m)
   K = sp.ellipkinc(2*pi, m)



   Factor_1Br = 2*(z1[w]-pos_espira_circular[j])/(radio_PMT*np.sqrt((radio_cilindro+radio_PMT)**2+(z1[w]-pos_espira_circular[j])**2))
   Factor_2Br = (radio_cilindro**2+radio_PMT**2+(z1[w]-pos_espira_circular[j])**2)/((radio_cilindro-radio_PMT)**2+(z1[w]-pos_espira_circular[j])**2)
   Br.append(Factor_1Br*(Factor_2Br*E-K))


   Factor_1Bz = 2/(np.sqrt((radio_cilindro+radio_PMT)**2+(z1[w]-pos_espira_circular[j])**2))
   Factor_2Bz = (radio_cilindro**2-radio_PMT**2-(z1[w]-pos_espira_circular[j])**2)/((radio_cilindro-radio_PMT)**2+(z1[w]-pos_espira_circular[j])**2)

   Bz.append(Factor_1Bz*(Factor_2Bz*E+K))


  Br_suma_espiras = np.sum(Br)
  Bz_suma_espiras = np.sum(Bz)



  Br_circular.append(Br_suma_espiras)
  Bz_circular.append(Bz_suma_espiras)


 corriente_circular = 366/np.max(Bz_circular)

 return corriente_circular

