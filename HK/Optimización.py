import numpy as np
import matplotlib.pyplot as plt
import wire
import biotsavart
import SuperKAngulos
import Corrientes


Altura = 72
pos_espira_rectangular = np.arange(-32,33,2)
radio_cilindro = 34
radio_tapas = 22.5
PMTs_vertical = np.arange(-34,34,0.3)
limite = 100*np.ones(len(PMTs_vertical))
ceros = np.zeros(len(PMTs_vertical))

radios_PMT = np.arange(0,33.4,1)
radio_PMT= 32.4
Angulo1 = 0
Angulo2 = np.pi/6
Angulo3 = np.pi/4
Angulo4 = np.pi/3
Angulo5 = np.pi/2

Angulo = np.arange(0, np.pi, 0.1)

I_rectangular = np.arange(0,100,20)
I_circular = np.arange(0,100,20)
Distancia_circulares = np.arange(1,10, 1)
Distancia_rectangulares = np.arange(1,10, 1)
PMTs_malos=[]
#Mallado de puntos

points1=[]
points2=[]
points3=[]
points4=[]
points5=[]

points_top1=[]
points_top2=[]
points_top3=[]
points_top4=[]
points_top5=[]

points_bottom1=[]
points_bottom2=[]
points_bottom3=[]
points_bottom4=[]
points_bottom5=[]


for k in range(len(PMTs_vertical)):
 points1.append([radio_PMT*np.cos(Angulo1), radio_PMT*np.sin(Angulo1) , PMTs_vertical[k]])
 points2.append([radio_PMT*np.cos(Angulo2), radio_PMT*np.sin(Angulo2) , PMTs_vertical[k]])
 points3.append([radio_PMT*np.cos(Angulo3), radio_PMT*np.sin(Angulo3) , PMTs_vertical[k]])
 points4.append([radio_PMT*np.cos(Angulo4), radio_PMT*np.sin(Angulo4) , PMTs_vertical[k]])
 points5.append([radio_PMT*np.cos(Angulo5), radio_PMT*np.sin(Angulo5) , PMTs_vertical[k]])

for r in range(len(radios_PMT)):
 points_top1.append([radios_PMT[r]*np.cos(Angulo1), radios_PMT[r]*np.sin(Angulo1) , 33])
 points_top2.append([radios_PMT[r]*np.cos(Angulo2), radios_PMT[r]*np.sin(Angulo2) , 33])
 points_top3.append([radios_PMT[r]*np.cos(Angulo3), radios_PMT[r]*np.sin(Angulo3) , 33])
 points_top4.append([radios_PMT[r]*np.cos(Angulo4), radios_PMT[r]*np.sin(Angulo4) , 33])
 points_top5.append([radios_PMT[r]*np.cos(Angulo5), radios_PMT[r]*np.sin(Angulo5) , 33])
 
 points_bottom1.append([radios_PMT[r]*np.cos(Angulo1), radios_PMT[r]*np.sin(Angulo1) , -34])
 points_bottom2.append([radios_PMT[r]*np.cos(Angulo2), radios_PMT[r]*np.sin(Angulo2) , -34])
 points_bottom3.append([radios_PMT[r]*np.cos(Angulo3), radios_PMT[r]*np.sin(Angulo3) , -34])
 points_bottom4.append([radios_PMT[r]*np.cos(Angulo4), radios_PMT[r]*np.sin(Angulo4) , -34])
 points_bottom5.append([radios_PMT[r]*np.cos(Angulo5), radios_PMT[r]*np.sin(Angulo5) , -34])



#Definición de constantes


for Ir in range(len(I_rectangular)):
 for Ic in range(len(I_circular)):
  for e in range(len(Distancia_rectangulares)):
   pos_espira_rectangular = np.arange(-32,32.1, Distancia_rectangulares[e])
   for c in range(len(Distancia_circulares)):
    pos_espira_circular = np.arange(-36,36.1, Distancia_circulares[c])
 

# rectangular loops I=1 A


    w1r = wire.Wire(path=wire.Wire.RectangularPath(dx=np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*2, dy=Altura), discretization_length=0.1, current=I_rectangular[Ir]).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[0], 0])
    sol = biotsavart.BiotSavart(wire=w1r)

    for i in range(1, len(pos_espira_rectangular)):
     lado_espira_rect = np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)
     w2r = wire.Wire(path=wire.Wire.RectangularPath(dx=lado_espira_rect*2, dy=Altura), discretization_length=0.1, current=I_rectangular[Ir]).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[i],0])
     sol.AddWire(w2r)



# circular loops I=1 A

    w1c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas, pts=20), discretization_length=0.1, current=I_circular[Ic]).Translate([0,0, pos_espira_circular[0]])
    #sol = biotsavart.BiotSavart(wire=w1c)
    sol.AddWire(w1c)

    for j in range(len(pos_espira_circular)):
     w2c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circular[Ic]).Translate([0,0,pos_espira_circular[j]])
     sol.AddWire(w2c)

    w3c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas, pts=20), discretization_length=0.1, current=I_circular[Ic]).Translate([0,0, pos_espira_circular[len(pos_espira_circular)-1]])
    sol.AddWire(w3c)





# calculate B field at given points
    B1 = sol.CalculateB(points=points1)*(10**7)
    B2 = sol.CalculateB(points=points2)*(10**7)
    B3 = sol.CalculateB(points=points3)*(10**7)
    B4 = sol.CalculateB(points=points4)*(10**7)
    B5 = sol.CalculateB(points=points5)*(10**7)

    B_top1 = sol.CalculateB(points=points_top1)*(10**7)
    B_top2 = sol.CalculateB(points=points_top2)*(10**7)
    B_top3 = sol.CalculateB(points=points_top3)*(10**7)
    B_top4 = sol.CalculateB(points=points_top4)*(10**7)
    B_top5 = sol.CalculateB(points=points_top5)*(10**7)
    
    B_bottom1 = sol.CalculateB(points=points_bottom1)*(10**7)
    B_bottom2 = sol.CalculateB(points=points_bottom2)*(10**7)
    B_bottom3 = sol.CalculateB(points=points_bottom3)*(10**7)
    B_bottom4 = sol.CalculateB(points=points_bottom4)*(10**7)
    B_bottom5 = sol.CalculateB(points=points_bottom5)*(10**7)






#Arrays por componentes

    Bx1=[]
    By1=[]
    Bz1=[]
    B_perp1=[]


    Bx2=[]
    By2=[]
    Bz2=[]
    B_perp2=[]


    Bx3=[]
    By3=[]
    Bz3=[]
    B_perp3=[]


    Bx4=[]
    By4=[]
    Bz4=[]
    B_perp4=[]


    Bx5=[]
    By5=[]
    Bz5=[]
    B_perp5=[]

    B_topx1=[]
    B_topy1=[]
    B_topz1=[]
    B_top_perp1=[]


    B_topx2=[]
    B_topy2=[]
    B_topz2=[]
    B_top_perp2=[]


    B_topx3=[]
    B_topy3=[]
    B_topz3=[]
    B_top_perp3=[]


    B_topx4=[]
    B_topy4=[]
    B_topz4=[]
    B_top_perp4=[]


    B_topx5=[]
    B_topy5=[]
    B_topz5=[]
    B_top_perp5=[]

    B_bottomx1=[]
    B_bottomy1=[]
    B_bottomz1=[]
    B_bottom_perp1=[]


    B_bottomx2=[]
    B_bottomy2=[]
    B_bottomz2=[]
    B_bottom_perp2=[]


    B_bottomx3=[]
    B_bottomy3=[]
    B_bottomz3=[]
    B_bottom_perp3=[]


    B_bottomx4=[]
    B_bottomy4=[]
    B_bottomz4=[]
    B_bottom_perp4=[]


    B_bottomx5=[]
    B_bottomy5=[]
    B_bottomz5=[]
    B_bottom_perp5=[]
    
    Numero_PMTs = 0




    for l in range(len(points1)):
     Bx1.append(B1[l,0])
     By1.append(B1[l,1]+303)		
     Bz1.append(B1[l,2]-366)
     B_perp1.append(np.sqrt((B1[l,0]*np.sin(Angulo1)-(B1[l,1]+303)*np.cos(Angulo1))**2+(B1[l,2]-366)**2))
 
     Bx2.append(B2[l,0])
     By2.append(B2[l,1]+303)		
     Bz2.append(B2[l,2]-366)
     B_perp2.append(np.sqrt((B2[l,0]*np.sin(Angulo2)-(B2[l,1]+303)*np.cos(Angulo2))**2+(B2[l,2]-366)**2))
 
     Bx3.append(B3[l,0])
     By3.append(B3[l,1]+303)		
     Bz3.append(B3[l,2]-366)
     B_perp3.append(np.sqrt((B3[l,0]*np.sin(Angulo3)-(B3[l,1]+303)*np.cos(Angulo3))**2+(B3[l,2]-366)**2))
 
     Bx4.append(B4[l,0])
     By4.append(B4[l,1]+303)		
     Bz4.append(B4[l,2]-366)
     B_perp4.append(np.sqrt((B4[l,0]*np.sin(Angulo4)-(B4[l,1]+303)*np.cos(Angulo4))**2+(B4[l,2]-366)**2))
 
     Bx5.append(B5[l,0])  
     By5.append(B5[l,1]+303)		
     Bz5.append(B5[l,2]-366)
     B_perp5.append(np.sqrt((B5[l,0]*np.sin(Angulo5)-(B5[l,1]+303)*np.cos(Angulo5))**2+(B5[l,2]-366)**2))
     
    
    for q in range(len(B_perp1)):
     if np.abs(B_perp1[q]) > 100:
      Numero_PMTs = Numero_PMTs + 1
     
     if np.abs(B_perp2[q]) > 100:
      Numero_PMTs = Numero_PMTs + 1
      
     if np.abs(B_perp3[q]) > 100:
      Numero_PMTs = Numero_PMTs + 1
      
     if np.abs(B_perp4[q]) > 100:
      Numero_PMTs = Numero_PMTs + 1
      
     if np.abs(B_perp5[q]) > 100:
      Numero_PMTs = Numero_PMTs + 1



    for s in range(len(points_top1)):
     B_topx1.append(B_top1[s,0])
     B_topy1.append(B_top1[s,1]+303)		
     B_topz1.append(B_top1[s,2]-366)
     B_top_perp1.append(np.sqrt(B_top1[s,0]**2+(B_top1[s,1]+303)**2))


     B_topx2.append(B_top2[s,0])
     B_topy2.append(B_top2[s,1]+303)		
     B_topz2.append(B_top2[s,2]-366)
     B_top_perp2.append(np.sqrt(B_top2[s,0]**2+(B_top2[s,1]+303)**2))


     B_topx3.append(B_top3[s,0])
     B_topy3.append(B_top3[s,1]+303)		
     B_topz3.append(B_top3[s,2]-366)
     B_top_perp3.append(np.sqrt(B_top3[s,0]**2+(B_top3[s,1]+303)**2))


     B_topx4.append(B_top4[s,0])
     B_topy4.append(B_top4[s,1]+303)		
     B_topz4.append(B_top4[s,2]-366)
     B_top_perp4.append(np.sqrt(B_top4[s,0]**2+(B_top4[s,1]+303)**2))


     B_topx5.append(B_top5[s,0])
     B_topy5.append(B_top5[s,1]+303)		
     B_topz5.append(B_top5[s,2]-366)
     B_top_perp5.append(np.sqrt(B_top5[s,0]**2+(B_top5[s,1]+303)**2))
 
     B_bottomx1.append(B_bottom1[s,0])
     B_bottomy1.append(B_bottom1[s,1]+303)		
     B_bottomz1.append(B_bottom1[s,2]-366)
     B_bottom_perp1.append(np.sqrt(B_bottom1[s,0]**2+(B_bottom1[s,1]+303)**2))


     B_bottomx2.append(B_bottom2[s,0])
     B_bottomy2.append(B_bottom2[s,1]+303)		
     B_bottomz2.append(B_bottom2[s,2]-366)
     B_bottom_perp2.append(np.sqrt(B_bottom2[s,0]**2+(B_bottom2[s,1]+303)**2))


     B_bottomx3.append(B_bottom3[s,0])
     B_bottomy3.append(B_bottom3[s,1]+303)		
     B_bottomz3.append(B_bottom3[s,2]-366)
     B_bottom_perp3.append(np.sqrt(B_bottom3[s,0]**2+(B_bottom3[s,1]+303)**2))


     B_bottomx4.append(B_bottom4[s,0])
     B_bottomy4.append(B_bottom4[s,1]+303)		
     B_bottomz4.append(B_bottom4[s,2]-366)
     B_bottom_perp4.append(np.sqrt(B_bottom4[s,0]**2+(B_bottom4[s,1]+303)**2))


     B_bottomx5.append(B_bottom5[s,0])
     B_bottomy5.append(B_bottom5[s,1]+303)		
     B_bottomz5.append(B_bottom5[s,2]-366)
     B_bottom_perp5.append(np.sqrt(B_bottom5[s,0]**2+(B_bottom5[s,1]+303)**2))
     
     
    for o in range(len(B_top_perp1)):
     if np.abs(B_top_perp1[o]) > 100:
      Numero_PMTs = Numero_PMTs + 1
     
     if np.abs(B_top_perp2[o]) > 100:
      Numero_PMTs = Numero_PMTs + 1
      
     if np.abs(B_top_perp3[o]) > 100:
      Numero_PMTs = Numero_PMTs + 1
      
     if np.abs(B_top_perp4[o]) > 100:
      Numero_PMTs = Numero_PMTs + 1
      
     if np.abs(B_top_perp5[o]) > 100:
      Numero_PMTs = Numero_PMTs + 1
     
     if np.abs(B_bottom_perp1[o]) > 100:
      Numero_PMTs = Numero_PMTs + 1
     
     if np.abs(B_bottom_perp2[o]) > 100:
      Numero_PMTs = Numero_PMTs + 1
      
     if np.abs(B_bottom_perp3[o]) > 100:
      Numero_PMTs = Numero_PMTs + 1
      
     if np.abs(B_bottom_perp4[o]) > 100:
      Numero_PMTs = Numero_PMTs + 1
      
     if np.abs(B_bottom_perp5[o]) > 100:
      Numero_PMTs = Numero_PMTs + 1

    PMTs_malos.append(Numero_PMTs)


Mejor_situacion = np.amin(PMTs_malos)
Donde = np.where(PMTs_malos == np.amin(PMTs_malos))

print(Mejor_situacion, Donde)

