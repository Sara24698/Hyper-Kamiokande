#! /home/sara/anaconda3/bin/python
from classy import Class

def Angulos():
    # calculate magnetic fields (in tesla, T)
    # with the Biot-Savart law
    # the output is points (coordinates x,y,z) and B(1,B2,B3) for the 3 cases (B1,B2,B3)
    import numpy as np
    import matplotlib.pyplot as plt
    import wire
    import biotsavart




    

    
    #Definición de constantes
    
    
    I_rectangular = 62.5#74.5111705813268#I[0]
    I_circulares = -75.87308111813685#I[1]
    I_topbottom = -75.87308111813685#I[1]
    
    Altura = 72
    pos_espira_rectangular = np.arange(-32,33,2)
    radio_cilindro = 34
    radio_tapas = 27
    radio_tapas2 = 20
    pos_espira_circular = np.arange(-35.5,37,2)
    PMTs_vertical = np.arange(-33.4,32.5,0.707)
    ceros = np.zeros(len(PMTs_vertical))
    
    radio_PMT_tapas = 31.95
    radio_PMT= 32.4
    Angulo1 = 0
    Angulo2 = np.pi/6
    Angulo3 = np.pi/4
    Angulo4 = np.pi/3
    Angulo5 = np.pi/2
    
    Angulo = np.arange(0, 6.315, 0.0218)
    z=np.arange(-32.522, 32.523, 0.707)
    limite = 100*np.ones(len(PMTs_vertical))
    Angulo_tapas=[]

    PMTs_malos_total = []
    PMTs_malos_paredes=[]
    Media = []
    Media_tapas=[]
    points2=[]
    PMTs_top=[]
    PMTs_bottom=[]
    Coordenadas_top=[]
    Coordenadas_bottom=[]
    Coordenadas_paredes=[]
    radios_PMT = np.arange(0, 31.95, 0.5)




    
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

   
	
    
   #Puntos en el eje Y
   
    for k in range(len(PMTs_vertical)):
     points1.append([radio_PMT*np.cos(Angulo1), radio_PMT*np.sin(Angulo1) , PMTs_vertical[k]])
     points2.append([radio_PMT*np.cos(Angulo2), radio_PMT*np.sin(Angulo2) , PMTs_vertical[k]])
     points3.append([radio_PMT*np.cos(Angulo3), radio_PMT*np.sin(Angulo3) , PMTs_vertical[k]])
     points4.append([radio_PMT*np.cos(Angulo4), radio_PMT*np.sin(Angulo4) , PMTs_vertical[k]])
     points5.append([radio_PMT*np.cos(Angulo5), radio_PMT*np.sin(Angulo5) , PMTs_vertical[k]])
   
    for r in range(len(radios_PMT)):
     points_top1.append([radios_PMT[r]*np.cos(Angulo1), radios_PMT[r]*np.sin(Angulo1) , 32.9])
     points_top2.append([radios_PMT[r]*np.cos(Angulo2), radios_PMT[r]*np.sin(Angulo2) , 32.9])
     points_top3.append([radios_PMT[r]*np.cos(Angulo3), radios_PMT[r]*np.sin(Angulo3) , 32.9])
     points_top4.append([radios_PMT[r]*np.cos(Angulo4), radios_PMT[r]*np.sin(Angulo4) , 32.9])
     points_top5.append([radios_PMT[r]*np.cos(Angulo5), radios_PMT[r]*np.sin(Angulo5) , 32.9])
     
     points_bottom1.append([radios_PMT[r]*np.cos(Angulo1), radios_PMT[r]*np.sin(Angulo1) , -32.9])
     points_bottom2.append([radios_PMT[r]*np.cos(Angulo2), radios_PMT[r]*np.sin(Angulo2) , -32.9])
     points_bottom3.append([radios_PMT[r]*np.cos(Angulo3), radios_PMT[r]*np.sin(Angulo3) , -32.9])
     points_bottom4.append([radios_PMT[r]*np.cos(Angulo4), radios_PMT[r]*np.sin(Angulo4) , -32.9])
     points_bottom5.append([radios_PMT[r]*np.cos(Angulo5), radios_PMT[r]*np.sin(Angulo5) , -32.9])

    


# rectangular loops I=1 A


     w1r = wire.Wire(path=wire.Wire.RectangularPath(dx=np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*2, dy=Altura), discretization_length=0.1, current=1.3*I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[0], 0]).Translate([0,0, 0.5])
     sol = biotsavart.BiotSavart(wire=w1r)



     for i in range(1, len(pos_espira_rectangular)-1):
      if i ==12:
       continue
      lado_espira_rect = np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)
      w2r = wire.Wire(path=wire.Wire.RectangularPath(dx=lado_espira_rect*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[i],0]).Translate([0,0, 0.5])
      sol.AddWire(w2r)

     w3r = wire.Wire(path=wire.Wire.RectangularPath(dx=np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*2, dy=Altura), discretization_length=0.1, current=1.47*I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[-1],0]).Translate([0,0, 0.5])
     sol.AddWire(w3r)


     w4r = wire.Wire(path=wire.Wire.RectangularPath(dx=np.sqrt(radio_cilindro**2 - pos_espira_rectangular[12]**2)*2, dy=Altura), discretization_length=0.1, current=64).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[12],0]).Translate([0,0, 0.5])
     sol.AddWire(w4r)




"""

# circular loops I=1 A


w1c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas, pts=20), discretization_length=0.1, current=-50).Translate([0,0, pos_espira_circular[0]])
#sol = biotsavart.BiotSavart(wire=w1c)
sol.AddWire(w1c)
Longitud.append(2*np.pi*radio_tapas)
Longitud.append(2*np.pi*radio_tapas)




w9c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
w10c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
w11c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
w12c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
w13c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
w14c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
sol.AddWire(w9c)
sol.AddWire(w10c)
sol.AddWire(w11c)
sol.AddWire(w12c)
sol.AddWire(w13c)
sol.AddWire(w14c)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)




w15c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
w16c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
w18c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
sol.AddWire(w15c)
sol.AddWire(w16c)
sol.AddWire(w18c)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)


for j in range(3, len(pos_espira_circular)-5):
w17c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[j]])
sol.AddWire(w17c)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)



w19c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-70).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-4]])
sol.AddWire(w19c)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)


w20c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-73).Translate([0,0,pos_espira_circular[1]])
sol.AddWire(w20c)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)

w21c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-60).Translate([0,0,pos_espira_circular[2]])
sol.AddWire(w21c)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)

w22c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-70).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-5]])
sol.AddWire(w22c)
Longitud.append(2*np.pi*radio_cilindro)
Longitud.append(2*np.pi*radio_cilindro)


#Parte de abajo
w1e = wire.Wire(path=wire.Wire.EllipticalPath(rx=9.05, ry=19.05, pts=20), discretization_length=0.1, current=-35).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0,-24.8,0])
sol.AddWire(w1e)
Longitud.append(2*np.pi*np.sqrt((9.05**2+19.05**2)/2))

w2e = wire.Wire(path=wire.Wire.EllipticalPath(rx=14.765, ry=24.77, pts=20), discretization_length=0.1, current=-90).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0,-19.205,0])
sol.AddWire(w2e)
Longitud.append(2*np.pi*np.sqrt((14.765**2+24.77**2)/2))
Longitud.append(2*np.pi*np.sqrt((14.765**2+24.77**2)/2))
Longitud.append(2*np.pi*np.sqrt((14.765**2+24.77**2)/2))

w3e = wire.Wire(path=wire.Wire.EllipticalPath(rx=20.48, ry=30.49, pts=20), discretization_length=0.1, current=-148).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0,-11.91,0])
sol.AddWire(w3e)
Longitud.append(2*np.pi*np.sqrt((20.48**2+30.49**2)/2))
Longitud.append(2*np.pi*np.sqrt((20.48**2+30.49**2)/2))
Longitud.append(2*np.pi*np.sqrt((20.48**2+30.49**2)/2))
Longitud.append(2*np.pi*np.sqrt((20.48**2+30.49**2)/2))

w4e = wire.Wire(path=wire.Wire.EllipticalPath(rx=9.05, ry=19.05, pts=20), discretization_length=0.1, current=75).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0,24.8,0])
sol.AddWire(w4e)
Longitud.append(2*np.pi*np.sqrt((9.05**2+19.05**2)/2))
Longitud.append(2*np.pi*np.sqrt((9.05**2+19.05**2)/2))

#Parte de arriba

w5e = wire.Wire(path=wire.Wire.EllipticalPath(rx=20.48, ry=30.49, pts=20), discretization_length=0.1, current=-150).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0, 11.91,0])
sol.AddWire(w5e)
Longitud.append(2*np.pi*np.sqrt((20.48**2+30.49**2)/2))
Longitud.append(2*np.pi*np.sqrt((20.48**2+30.49**2)/2))
Longitud.append(2*np.pi*np.sqrt((20.48**2+30.49**2)/2))
Longitud.append(2*np.pi*np.sqrt((20.48**2+30.49**2)/2))

w6e = wire.Wire(path=wire.Wire.EllipticalPath(rx=14.765, ry=24.77, pts=20), discretization_length=0.1, current=-110).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0,19.205,0])
sol.AddWire(w6e)
Longitud.append(2*np.pi*np.sqrt((14.765**2+24.77**2)/2))
Longitud.append(2*np.pi*np.sqrt((14.765**2+24.77**2)/2))
Longitud.append(2*np.pi*np.sqrt((14.765**2+24.77**2)/2))

w7e = wire.Wire(path=wire.Wire.EllipticalPath(rx=9.05, ry=19.05, pts=20), discretization_length=0.1, current=-25).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0,24.8,0])
sol.AddWire(w7e)
Longitud.append(2*np.pi*np.sqrt((9.05**2+19.05**2)/2))


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







    for l in range(len(points1)):
     Bx1.append(B1[l,0])
     By1.append(B1[l,1]+303)		
     Bz1.append(B1[l,2]-366)
     B_perp1.append(np.sqrt((B1[l,0]*np.sin(Angulo1)-(B1[l,1]+303)*np.cos(Angulo1))**2+(B1[l,2]-366)**2))
     Media.append(np.sqrt((B1[l,0]*np.sin(Angulo1)-(B1[l,1]+303)*np.cos(Angulo1))**2+(B1[l,2]-366)**2))

    for l in range(len(points2)): 
     Bx2.append(B2[l,0])
     By2.append(B2[l,1]+303)		
     Bz2.append(B2[l,2]-366)
     B_perp2.append(np.sqrt((B2[l,0]*np.sin(Angulo2)-(B2[l,1]+303)*np.cos(Angulo2))**2+(B2[l,2]-366)**2))
     Media.append(np.sqrt((B2[l,0]*np.sin(Angulo2)-(B2[l,1]+303)*np.cos(Angulo2))**2+(B2[l,2]-366)**2))

     Bx3.append(B3[l,0])
     By3.append(B3[l,1]+303)		
     Bz3.append(B3[l,2]-366)
     B_perp3.append(np.sqrt((B3[l,0]*np.sin(Angulo3)-(B3[l,1]+303)*np.cos(Angulo3))**2+(B3[l,2]-366)**2))
     Media.append(np.sqrt((B3[l,0]*np.sin(Angulo3)-(B3[l,1]+303)*np.cos(Angulo3))**2+(B3[l,2]-366)**2))

     Bx4.append(B4[l,0])
     By4.append(B4[l,1]+303)		
     Bz4.append(B4[l,2]-366)
     B_perp4.append(np.sqrt((B4[l,0]*np.sin(Angulo4)-(B4[l,1]+303)*np.cos(Angulo4))**2+(B4[l,2]-366)**2))
     Media.append(np.sqrt((B4[l,0]*np.sin(Angulo4)-(B4[l,1]+303)*np.cos(Angulo4))**2+(B4[l,2]-366)**2))

     Bx5.append(B5[l,0])
     By5.append(B5[l,1]+303)		
     Bz5.append(B5[l,2]-366)
     B_perp5.append(np.sqrt((B5[l,0]*np.sin(Angulo5)-(B5[l,1]+303)*np.cos(Angulo5))**2+(B5[l,2]-366)**2))
     Media.append(np.sqrt((B5[l,0]*np.sin(Angulo5)-(B5[l,1]+303)*np.cos(Angulo5))**2+(B5[l,2]-366)**2))


    for s in range(len(points_top1)):
     B_topx1.append(B_top1[s,0])
     B_topy1.append(B_top1[s,1]+303)		
     B_topz1.append(B_top1[s,2]-366)
     B_top_perp1.append(np.sqrt(B_top1[s,0]**2+(B_top1[s,1]+303)**2))
     Media.append(np.sqrt(B_top1[s,0]**2+(B_top1[s,1]+303)**2))


     B_topx2.append(B_top2[s,0])
     B_topy2.append(B_top2[s,1]+303)		
     B_topz2.append(B_top2[s,2]-366)
     B_top_perp2.append(np.sqrt(B_top2[s,0]**2+(B_top2[s,1]+303)**2))
     Media.append(np.sqrt(B_top2[s,0]**2+(B_top2[s,1]+303)**2))


     B_topx3.append(B_top3[s,0])
     B_topy3.append(B_top3[s,1]+303)		
     B_topz3.append(B_top3[s,2]-366)
     B_top_perp3.append(np.sqrt(B_top3[s,0]**2+(B_top3[s,1]+303)**2))
     Media.append(np.sqrt(B_top3[s,0]**2+(B_top3[s,1]+303)**2))


     B_topx4.append(B_top4[s,0])
     B_topy4.append(B_top4[s,1]+303)		
     B_topz4.append(B_top4[s,2]-366)
     B_top_perp4.append(np.sqrt(B_top4[s,0]**2+(B_top4[s,1]+303)**2))
     Media.append(np.sqrt(B_top4[s,0]**2+(B_top4[s,1]+303)**2))


     B_topx5.append(B_top5[s,0]) 
     B_topy5.append(B_top5[s,1]+303)		
     B_topz5.append(B_top5[s,2]-366)
     B_top_perp5.append(np.sqrt(B_top5[s,0]**2+(B_top5[s,1]+303)**2))
     Media.append(np.sqrt(B_top5[s,0]**2+(B_top5[s,1]+303)**2))

     B_bottomx1.append(B_bottom1[s,0])
     B_bottomy1.append(B_bottom1[s,1]+303)		
     B_bottomz1.append(B_bottom1[s,2]-366)
     B_bottom_perp1.append(np.sqrt(B_bottom1[s,0]**2+(B_bottom1[s,1]+303)**2))
     Media.append(np.sqrt(B_bottom1[s,0]**2+(B_bottom1[s,1]+303)**2))


     B_bottomx2.append(B_bottom2[s,0])
     B_bottomy2.append(B_bottom2[s,1]+303)		
     B_bottomz2.append(B_bottom2[s,2]-366)
     B_bottom_perp2.append(np.sqrt(B_bottom2[s,0]**2+(B_bottom2[s,1]+303)**2))
     Media.append(np.sqrt(B_bottom2[s,0]**2+(B_bottom2[s,1]+303)**2))


     B_bottomx3.append(B_bottom3[s,0])
     B_bottomy3.append(B_bottom3[s,1]+303)		
     B_bottomz3.append(B_bottom3[s,2]-366)
     B_bottom_perp3.append(np.sqrt(B_bottom3[s,0]**2+(B_bottom3[s,1]+303)**2))
     Media.append(np.sqrt(B_bottom3[s,0]**2+(B_bottom3[s,1]+303)**2))


     B_bottomx4.append(B_bottom4[s,0])
     B_bottomy4.append(B_bottom4[s,1]+303)		
     B_bottomz4.append(B_bottom4[s,2]-366)
     B_bottom_perp4.append(np.sqrt(B_bottom4[s,0]**2+(B_bottom4[s,1]+303)**2))
     Media.append(np.sqrt(B_bottom4[s,0]**2+(B_bottom4[s,1]+303)**2))


     B_bottomx5.append(B_bottom5[s,0])
     B_bottomy5.append(B_bottom5[s,1]+303)		
     B_bottomz5.append(B_bottom5[s,2]-366)
     B_bottom_perp5.append(np.sqrt(B_bottom5[s,0]**2+(B_bottom5[s,1]+303)**2))
     Media.append(np.sqrt(B_bottom5[s,0]**2+(B_bottom5[s,1]+303)**2))




# make figure with loops in 3D

    Media_puntos = np.sum(Media)/680
    print('La media del campo magnético es de', Media_puntos)
  



    PMTs_malos = 0

    for alfa in range(len(B_perp1)):
     if np.abs(B_perp1[alfa]) > 100:
      PMTs_malos = PMTs_malos+1

     if np.abs(B_perp2[alfa]) > 100:
      PMTs_malos = PMTs_malos+1

     if np.abs(B_perp3[alfa]) > 100:
      PMTs_malos = PMTs_malos+1

     if np.abs(B_perp4[alfa]) > 100:
      PMTs_malos = PMTs_malos+1

     if np.abs(B_perp5[alfa]) > 100:
      PMTs_malos = PMTs_malos+1    

    for beta in range(len(B_top_perp1)):
     if np.abs(B_top_perp1[beta]) > 100:
      PMTs_malos = PMTs_malos+1

     if np.abs(B_top_perp2[beta]) > 100:
      PMTs_malos = PMTs_malos+1

     if np.abs(B_top_perp3[beta]) > 100:
      PMTs_malos = PMTs_malos+1

     if np.abs(B_top_perp4[beta]) > 100:
      PMTs_malos = PMTs_malos+1

     if np.abs(B_top_perp5[beta]) > 100:
      PMTs_malos = PMTs_malos+1

     if np.abs(B_bottom_perp1[beta]) > 100:
      PMTs_malos = PMTs_malos+1

     if np.abs(B_bottom_perp2[beta]) > 100:
      PMTs_malos = PMTs_malos+1

     if np.abs(B_bottom_perp3[beta]) > 100:
      PMTs_malos = PMTs_malos+1

     if np.abs(B_bottom_perp4[beta]) > 100:
      PMTs_malos = PMTs_malos+1

     if np.abs(B_bottom_perp5[beta]) > 100:
      PMTs_malos = PMTs_malos+1


    print('El número de PMTs malos es', PMTs_malos)
    

"""
     fig = plt.figure()
     ax = fig.gca(projection='3d')
     sol.mpl3d_PlotWires(ax)
     plt.show()
    
"""
    
    #Otras figuras
    
    
    
    #Figuras paredes

    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(PMTs_vertical, B_perp1, label='$\Delta B_{perp}$ $0 \: rad$')
    ax.scatter(PMTs_vertical, B_perp2, label='$\Delta B_{perp}$ $\pi/6 \: rad$')
    ax.scatter(PMTs_vertical, B_perp3, label='$\Delta B_{perp}$ $\pi/4 \: rad$')
    ax.scatter(PMTs_vertical, B_perp4, label='$\Delta B_{perp}$ $\pi/3 \: rad$')
    ax.scatter(PMTs_vertical, B_perp5, label='$\Delta B_{perp}$ $\pi/2 \: rad$')
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
    ax.scatter(PMTs_vertical, Bz1, label='$\Delta B_{z}$ $0 \: rad$')
    ax.scatter(PMTs_vertical, Bz2, label='$\Delta B_{z}$ $\pi/6 \: rad$')
    ax.scatter(PMTs_vertical, Bz3, label='$\Delta B_{z}$ $\pi/4 \: rad$')
    ax.scatter(PMTs_vertical, Bz4, label='$\Delta B_{z}$ $\pi/3 \: rad$')
    ax.scatter(PMTs_vertical, Bz5, label='$\Delta B_{z}$ $\pi/2 \: rad$')
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
    ax.scatter(PMTs_vertical, Bx1, label='$B_{x}$ $0 \: rad$')
    ax.scatter(PMTs_vertical, Bx2, label='$B_{x}$ $\pi/6 \: rad$')
    ax.scatter(PMTs_vertical, Bx3, label='$B_{x}$ $\pi/4 \: rad$')
    ax.scatter(PMTs_vertical, Bx4, label='$B_{x}$ $\pi/3 \: rad$')
    ax.scatter(PMTs_vertical, Bx5, label='$B_{x}$ $\pi/2 \: rad$')
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
    ax.scatter(PMTs_vertical, By1, label='$\Delta B_{y}$ $0 \: rad$')
    ax.scatter(PMTs_vertical, By2, label='$\Delta B_{y}$ $\pi/6 \: rad$')
    ax.scatter(PMTs_vertical, By3, label='$\Delta B_{y}$ $\pi/4 \: rad$')
    ax.scatter(PMTs_vertical, By4, label='$\Delta B_{y}$ $\pi/3 \: rad$')
    ax.scatter(PMTs_vertical, By5, label='$\Delta B_{y}$ $\pi/2 \: rad$')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Height (m)')
    plt.ylabel('B (mG)')
    #plt.xlim(0,35)
    plt.title('$\Delta B_{y}$ at walls as a function of height and angle')
    plt.legend()
    plt.savefig('By_walls.png')
    plt.show()
    
    
    
    #Figuras en las tapas
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(radios_PMT, B_top_perp1, label='$\Delta B_{perp}$ $0 \: rad$')
    ax.scatter(radios_PMT, B_top_perp2, label='$\Delta B_{perp}$ $\pi/6 \: rad$')
    ax.scatter(radios_PMT, B_top_perp3, label='$\Delta B_{perp}$ $\pi/4 \: rad$')
    ax.scatter(radios_PMT, B_top_perp4, label='$\Delta B_{perp}$ $\pi/3 \: rad$')
    ax.scatter(radios_PMT, B_top_perp5, label='$\Delta B_{perp}$ $\pi/2 \: rad$')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Radius (m)')
    plt.ylabel('B (mG)')
    plt.xlim(0,35)
    plt.title('$\Delta B_{perp}$ at top as a function of radius and angle')
    plt.legend()
    plt.savefig('Bperp_top.png')
    plt.show()
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(radios_PMT, B_topz1, label='$\Delta B_{z}$ $0 \: rad$')
    ax.scatter(radios_PMT, B_topz2, label='$\Delta B_{z}$ $\pi/6 \: rad$')
    ax.scatter(radios_PMT, B_topz3, label='$\Delta B_{z}$ $\pi/4 \: rad$')
    ax.scatter(radios_PMT, B_topz4, label='$\Delta B_{z}$ $\pi/3 \: rad$')
    ax.scatter(radios_PMT, B_topz5, label='$\Delta B_{z}$ $\pi/2 \: rad$')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Radius (m)')
    plt.ylabel('B (mG)')
    plt.xlim(0,35)
    plt.title('$\Delta B_{z}$ at top as a function of radius and angle')
    plt.legend()
    plt.savefig('Bz_top.png')
    plt.show()
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(radios_PMT, B_topy1, label='$\Delta B_{y}$ $0 \: rad$')
    ax.scatter(radios_PMT, B_topy2, label='$\Delta B_{y}$ $\pi/6 \: rad$')
    ax.scatter(radios_PMT, B_topy3, label='$\Delta B_{y}$ $\pi/4 \: rad$')
    ax.scatter(radios_PMT, B_topy4, label='$\Delta B_{y}$ $\pi/3 \: rad$')
    ax.scatter(radios_PMT, B_topy5, label='$\Delta B_{y}$ $\pi/2 \: rad$')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Radius (m)')
    plt.ylabel('B (mG)')
    plt.xlim(0,35)
    plt.title('$\Delta B_{y}$ at top as a function of radius and angle')
    plt.legend()
    plt.savefig('By_top.png')
    plt.show()
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(radios_PMT, B_topx1, label='$B_{x}$ $0 \: rad$')
    ax.scatter(radios_PMT, B_topx2, label='$B_{x}$ $\pi/6 \: rad$')
    ax.scatter(radios_PMT, B_topx3, label='$B_{x}$ $\pi/4 \: rad$')
    ax.scatter(radios_PMT, B_topx4, label='$B_{x}$ $\pi/3 \: rad$')
    ax.scatter(radios_PMT, B_topx5, label='$B_{x}$ $\pi/2 \: rad$')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Radius (m)')
    plt.ylabel('B (mG)')
    plt.xlim(0,35)
    plt.title('$B_{x}$ at top as a function of radius and angle')
    plt.legend()
    plt.savefig('Bx_top.png')
    plt.show()
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(radios_PMT, B_bottom_perp1, label='$\Delta B_{perp}$ $0 \: rad$')
    ax.scatter(radios_PMT, B_bottom_perp2, label='$\Delta B_{perp}$ $\pi/6 \: rad$')
    ax.scatter(radios_PMT, B_bottom_perp3, label='$\Delta B_{perp}$ $\pi/4 \: rad$')
    ax.scatter(radios_PMT, B_bottom_perp4, label='$\Delta B_{perp}$ $\pi/3 \: rad$')
    ax.scatter(radios_PMT, B_bottom_perp5, label='$\Delta B_{perp}$ $\pi/2 \: rad$')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Radius (m)')
    plt.ylabel('B (mG)')
    plt.xlim(0,35)
    plt.title('$\Delta B_{perp}$ at bottom as a function of radius and angle')
    plt.legend()
    plt.savefig('Bperp_bottom.png')
    plt.show()
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(radios_PMT, B_bottomz1, label='$\Delta B_{z}$ $0 \: rad$')
    ax.scatter(radios_PMT, B_bottomz2, label='$\Delta B_{z}$ $\pi/6 \: rad$')
    ax.scatter(radios_PMT, B_bottomz3, label='$\Delta B_{z}$ $\pi/4 \: rad$')
    ax.scatter(radios_PMT, B_bottomz4, label='$\Delta B_{z}$ $\pi/3 \: rad$')
    ax.scatter(radios_PMT, B_bottomz5, label='$\Delta B_{z}$ $\pi/2 \: rad$')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Radius (m)')
    plt.ylabel('B (mG)')
    plt.xlim(0,35)
    plt.title('$\Delta B_{z}$ at bottom as a function of radius and angle')
    plt.legend()
    plt.savefig('Bz_bottom.png')
    plt.show()
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(radios_PMT, B_bottomy1, label='$\Delta B_{y}$ $0 \: rad$')
    ax.scatter(radios_PMT, B_bottomy2, label='$\Delta B_{y}$ $\pi/6 \: rad$')
    ax.scatter(radios_PMT, B_bottomy3, label='$\Delta B_{y}$ $\pi/4 \: rad$')
    ax.scatter(radios_PMT, B_bottomy4, label='$\Delta B_{y}$ $\pi/3 \: rad$')
    ax.scatter(radios_PMT, B_bottomy5, label='$\Delta B_{y}$ $\pi/2 \: rad$')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Radius (m)')
    plt.ylabel('B (mG)')
    plt.xlim(0,35)
    plt.title('$\Delta B_{y}$ at bottom as a function of radius and angle')
    plt.legend()
    plt.savefig('By_bottom.png')
    plt.show()
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(radios_PMT, B_bottomx1, label='$B_{x}$ $0 \: rad$')
    ax.scatter(radios_PMT, B_bottomx2, label='$B_{x}$ $\pi/6 \: rad$')
    ax.scatter(radios_PMT, B_bottomx3, label='$B_{x}$ $\pi/4 \: rad$')
    ax.scatter(radios_PMT, B_bottomx4, label='$B_{x}$ $\pi/3 \: rad$')
    ax.scatter(radios_PMT, B_bottomx5, label='$B_{x}$ $\pi/2 \: rad$')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Radius (m)')
    plt.ylabel('B (mG)')
    plt.xlim(0,35)
    plt.title('$B_{x}$ at bottom as a function of radius and angle')
    plt.legend()
    plt.savefig('Bx_bottom.png')
    plt.show()

	
    
    return PMTs_vertical, B_perp1

"""

Angulos()
