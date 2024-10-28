# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

rebar=pd.read_csv('./rebar_effect_.csv')

rebar = np.array(rebar)
Coordenadas=[]
Bx = []
By = []
Bz = []

for j in range(len(rebar)):
    Prueba=rebar[j][0]
    Prueba = Prueba.replace(';', ' ')
    Prueba = Prueba.split()
    Prueba = [eval(i) for i in Prueba]
    Prueba = np.array(Prueba)
    Coordenadas.append([Prueba[0], Prueba[1], Prueba[2], Prueba[3], Prueba[4], Prueba[5], Prueba[6], 1 if Prueba[2] == 32.9 else 3 if Prueba[2] == -32.9 else 2])
    

for k in range(len(Coordenadas)):
    if Coordenadas[k][3] < 100 and Coordenadas[k][3] > -100:
        Bx.append(Coordenadas[k][3])
        
    if Coordenadas[k][4] < 1000 and Coordenadas[k][4] > -1000:
        By.append(Coordenadas[k][4]-303)
        
    if Coordenadas [k][5] > -1000 and Coordenadas [k][5] < 1000:
        Bz.append(Coordenadas[k][5]+366)
    


#df = pd.DataFrame(Coordenadas, columns=['x', 'y', 'z', 'Bx', 'By', 'Bz', 'Bp', 'faceid'])
#df.to_csv('./rebar.csv', index = False)

std_Bx = []
mean_Bx = np.sum(Bx)/len(Bx)
for i in range(len(Bx)):
	std_Bx.append((Bx[i]-mean_Bx)**2)
        

    
print(np.sum(Bx)/len(Bx))
Desviacion_Bx = np.sqrt(np.sum(std_Bx)/len(Bx))    

intervalos = np.arange(-8, 10, 1) #calculamos los extremos de los intervalos
plt.hist(By, bins=intervalos, rwidth=0.85)
plt.xlabel("Bx (mG)")
plt.ylabel("Number of PMTs")
plt.ylim(0,15000)
#plt.xticks(np.arange(0, 200, 25))
#plt.title("$B_{perp}$ distribution for all the PMTs")
textstr = '\n'.join((r'$\mu=%.4f$' % (mean_Bx, ), r'$\sigma=%.2f$' % (Desviacion_Bx, )))
    
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

plt.text(5, 12000, textstr, fontsize=10, bbox=props)
plt.show()