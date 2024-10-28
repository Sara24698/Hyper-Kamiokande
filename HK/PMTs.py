def Angulos():
    # calculate magnetic fields (in tesla, T)
    # with the Biot-Savart law
    # the output is points (coordinates x,y,z) and B(1,B2,B3) for the 3 cases (B1,B2,B3)
    import numpy as np
    import matplotlib.pyplot as plt
    import wire
    import biotsavart
    import pandas as pd
    from mpl_toolkits.mplot3d import Axes3D



    

    
    #Definicion de constantes
    
    
    #I = Corrientes.Corriente()
    I_rectangular = 62.5
    I_circulares = -35
    I_topbottom = -35

    Altura = 72
    pos_espira_rectangular = np.arange(-32,33,2)
    radio_cilindro = 34
    radio_tapas = 27
    radio_tapas2 = 20
    pos_espira_circular = np.arange(-35.5,37,1)
    PMTs_vertical = np.arange(-33.4,32.5,0.707)
    ceros = np.zeros(len(PMTs_vertical))


    radio_PMT_tapas = 31.95
    radio_PMT= 32.4
    Angulo = np.arange(0, 6.315, 0.0218)
    z=np.arange(-32.522, 32.523, 0.707)
    limite = 100*np.ones(len(Angulo))


    PMTs_malos_total = []
    PMTs_malos_paredes=[]
    Media = []
    Media_top=[]
    Media_bottom=[]
    points2=[]
    PMTs_top=[]
    PMTs_bottom=[]
    Coordenadas_top=[]
    Coordenadas_bottom=[]
    Coordenadas_paredes=[]
    Longitud=[]
    Longitud1 = []
    Longitud2 = []
    Longitud3 = []
    
    x = np.arange(-31.815, 32, 0.707)
    y = np.arange(-31.815, 32, 0.707)

    for r in range(len(x)):
     for h in range(len(y)):
      points2.append([x[r], y[h]])
    
    
    for g in range(len(points2)):
     distancia = np.sqrt(points2[g][0]**2+points2[g][1]**2)
     if distancia <=32.01:
      PMTs_top.append([points2[g][0], points2[g][1], 32.9])
      PMTs_bottom.append([points2[g][0], points2[g][1], -32.9])
      
 
 
# rectangular loops I=1 A


    w1r = wire.Wire(path=wire.Wire.RectangularPath(dx=np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*2, dy=Altura), discretization_length=0.1, current=1.3*I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[0], 0]).Translate([0,0, 0.5])
    sol = biotsavart.BiotSavart(wire=w1r)


    for i in range(1, len(pos_espira_rectangular)-1):
     lado_espira_rect = np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)
     w2r = wire.Wire(path=wire.Wire.RectangularPath(dx=lado_espira_rect*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[i],0]).Translate([0,0, 0.5])
     sol.AddWire(w2r)


    w3r = wire.Wire(path=wire.Wire.RectangularPath(dx=np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*2, dy=Altura), discretization_length=0.1, current=1.47*I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[-1],0]).Translate([0,0, 0.5])#1.47I
    sol.AddWire(w3r)


# circular loops I=1 A

    w1c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas, pts=20), discretization_length=0.1, current=-53).Translate([0,0, pos_espira_circular[0]])
    #sol = biotsavart.BiotSavart(wire=w1c)
    sol.AddWire(w1c)




    w9c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
    w10c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
    w11c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
    w12c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
    w13c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
    w14c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
    sol.AddWire(w9c)
    #sol = biotsavart.BiotSavart(wire=w9c)

    sol.AddWire(w10c)
    sol.AddWire(w11c)
    sol.AddWire(w12c)
    sol.AddWire(w13c)
    sol.AddWire(w14c)




    w15c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
    w16c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
    w18c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-80).Translate([0,0,pos_espira_circular[0]])#I_circulares
    #sol.AddWire(w15c)

    sol.AddWire(w16c)
    sol.AddWire(w18c)



    for j in range(3, len(pos_espira_circular)-5):
     if j==15:
      continue
     if j == 54:
      continue
     if j == 67:
      continue
 
     w17c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[j]])
     sol.AddWire(w17c)

    w19c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-57).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-4]])#I=70
    sol.AddWire(w19c)



    w20c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-73).Translate([0,0,pos_espira_circular[1]])
    sol.AddWire(w20c)


    w21c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-60).Translate([0,0,pos_espira_circular[2]])
    sol.AddWire(w21c)


    w22c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-85).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-5]])
    sol.AddWire(w22c)





#Parte de abajo

    w1e = wire.Wire(path=wire.Wire.EllipticalPath(rx=9.05, ry=19.05, pts=20), discretization_length=0.1, current=-40).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0,-24.7,0])
    sol.AddWire(w1e)


    w2e = wire.Wire(path=wire.Wire.EllipticalPath(rx=14.765, ry=24.77, pts=20), discretization_length=0.1, current=-90).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0,-18.7,0])#I = 80
    sol.AddWire(w2e)


    w3e = wire.Wire(path=wire.Wire.EllipticalPath(rx=23.5, ry=30.4, pts=20), discretization_length=0.1, current=-110).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0, -9.6,0])#I = 130
    sol.AddWire(w3e)


    w4e = wire.Wire(path=wire.Wire.EllipticalPath(rx=9.05, ry=19.05, pts=20), discretization_length=0.1, current=80).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0,24.7,0])#I = 80
    sol.AddWire(w4e)
 





#Parte de arriba

    w5e = wire.Wire(path=wire.Wire.EllipticalPath(rx=23.7, ry=31.4, pts=20), discretization_length=0.1, current=-160).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0, 8.5,0])
    sol.AddWire(w5e)


    w6e = wire.Wire(path=wire.Wire.EllipticalPath(rx=14.765, ry=24.77, pts=20), discretization_length=0.1, current=-120).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0,18.7,0])#I = 120
    sol.AddWire(w6e)


    w7e = wire.Wire(path=wire.Wire.EllipticalPath(rx=9.05, ry=19.05, pts=20), discretization_length=0.1, current=-20).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0,24.7,0])#I =50
    sol.AddWire(w7e)
    
   
    # calculate B field at given points
    B2 = sol.CalculateB(points=PMTs_top)*(10**7)
    B3 = sol.CalculateB(points=PMTs_bottom)*(10**7)

  

    #Arrays por componentes
    
     
    Bx2=[]
    By2=[]
    Bz2=[]
    B_perp2=[]
     
    Bx3=[]
    By3=[]
    Bz3=[]
    B_perp3=[]
    

 
    
  
    


    for q in range(len(PMTs_top)):
     Bx2.append(B2[q,0])
     By2.append(B2[q,1]+303)		
     Bz2.append(B2[q,2]-366)
     B_perp2.append(np.sqrt(B2[q,0]**2+(B2[q,1]+303)**2))
     Media_top.append(np.sqrt(B2[q,0]**2+(B2[q,1]+303)**2))
     
      
     Bx3.append(B3[q,0])
     By3.append(B3[q,1]+303)		
     Bz3.append(B3[q,2]-366)
     B_perp3.append(np.sqrt(B3[q,0]**2+(B3[q,1]+303)**2))
     Media_bottom.append(np.sqrt(B3[q,0]**2+(B3[q,1]+303)**2))
     

    PMTs_malos_top = 0
    PMTs_malos_bottom = 0
    
    
    for beta in range(len(B_perp3)):
      if np.abs(B_perp2[beta]) > 100:
       PMTs_malos_top = PMTs_malos_top+1
       #Coordenadas_top.append([PMTs_top[beta], Bx2[beta], By2[beta], Bz2[beta], B_perp2[beta]])
       

      
      if np.abs(B_perp3[beta]) > 100:
       PMTs_malos_bottom = PMTs_malos_bottom+1
       #Coordenadas_bottom.append([PMTs_bottom[beta], Bx3[beta], By3[beta], Bz3[beta], B_perp3[beta]])

     
    print('El numero de PMTs malos en top es', PMTs_malos_top)
    print('El numero de PMTs malos en bottom es', PMTs_malos_bottom)
    
    for mu in range(len(PMTs_top)):
     Coordenadas_top.append([PMTs_top[mu], Bx2[mu], By2[mu], Bz2[mu], B_perp2[mu]])
     Coordenadas_bottom.append([PMTs_bottom[mu], Bx3[mu], By3[mu], Bz3[mu], B_perp3[mu]])


   #Mallado de puntos

    for i in range(len(z)):
     points = []
     for j in range(len(Angulo)):
      points.append([radio_PMT*np.cos(Angulo[j]), radio_PMT*np.sin(Angulo[j]), z[i]])
   

# rectangular loops I=1 A


     w1r = wire.Wire(path=wire.Wire.RectangularPath(dx=np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*2, dy=Altura), discretization_length=0.1, current=1.3*I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[0], 0]).Translate([0,0, 0.5])
     sol = biotsavart.BiotSavart(wire=w1r)


     for i in range(1, len(pos_espira_rectangular)-1):
      lado_espira_rect = np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)
      w2r = wire.Wire(path=wire.Wire.RectangularPath(dx=lado_espira_rect*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[i],0]).Translate([0,0, 0.5])
      sol.AddWire(w2r)


     w3r = wire.Wire(path=wire.Wire.RectangularPath(dx=np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*2, dy=Altura), discretization_length=0.1, current=1.47*I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[-1],0]).Translate([0,0, 0.5])#1.47I
     sol.AddWire(w3r)

# circular loops I=1 A

     w1c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas, pts=20), discretization_length=0.1, current=-53).Translate([0,0, pos_espira_circular[0]])
     #sol = biotsavart.BiotSavart(wire=w1c)
     sol.AddWire(w1c)




     w9c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
     w10c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
     w11c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
     w12c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
     w13c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
     w14c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
     sol.AddWire(w9c)
     #sol = biotsavart.BiotSavart(wire=w9c)

     sol.AddWire(w10c)
     sol.AddWire(w11c)
     sol.AddWire(w12c)
     sol.AddWire(w13c)
     sol.AddWire(w14c)




     w15c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
     w16c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
     w18c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-80).Translate([0,0,pos_espira_circular[0]])#I_circulares
     #sol.AddWire(w15c)

     sol.AddWire(w16c)
     sol.AddWire(w18c)



     for j in range(3, len(pos_espira_circular)-5):
      if j==15:
       continue
      if j == 54:
       continue
      if j == 67:
       continue

      w17c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[j]])
      sol.AddWire(w17c)

     w19c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-57).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-4]])#I=70
     sol.AddWire(w19c)



     w20c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-73).Translate([0,0,pos_espira_circular[1]])
     sol.AddWire(w20c)


     w21c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-60).Translate([0,0,pos_espira_circular[2]])
     sol.AddWire(w21c)


     w22c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-85).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-5]])
     sol.AddWire(w22c)





#Parte de abajo

     w1e = wire.Wire(path=wire.Wire.EllipticalPath(rx=9.05, ry=19.05, pts=20), discretization_length=0.1, current=-40).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0,-24.7,0])
     sol.AddWire(w1e)


     w2e = wire.Wire(path=wire.Wire.EllipticalPath(rx=14.765, ry=24.77, pts=20), discretization_length=0.1, current=-90).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0,-18.7,0])#I = 80
     sol.AddWire(w2e)


     w3e = wire.Wire(path=wire.Wire.EllipticalPath(rx=23.5, ry=30.4, pts=20), discretization_length=0.1, current=-110).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0, -9.6,0])#I = 130
     sol.AddWire(w3e)


     w4e = wire.Wire(path=wire.Wire.EllipticalPath(rx=9.05, ry=19.05, pts=20), discretization_length=0.1, current=80).Rotate(axis=(0,0,1), deg=90).Translate([0,0,-35.5]).Translate([0,24.7,0])#I = 80
     sol.AddWire(w4e)
 





#Parte de arriba

     w5e = wire.Wire(path=wire.Wire.EllipticalPath(rx=23.7, ry=31.4, pts=20), discretization_length=0.1, current=-160).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0, 8.5,0])
     sol.AddWire(w5e)


     w6e = wire.Wire(path=wire.Wire.EllipticalPath(rx=14.765, ry=24.77, pts=20), discretization_length=0.1, current=-120).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0,18.7,0])#I = 120
     sol.AddWire(w6e)


     w7e = wire.Wire(path=wire.Wire.EllipticalPath(rx=9.05, ry=19.05, pts=20), discretization_length=0.1, current=-20).Rotate(axis=(0,0,1), deg=90).Translate([0,0,36.5]).Translate([0,24.7,0])#I =50
     sol.AddWire(w7e)


    
    # calculate B field at given points
     B1 = sol.CalculateB(points=points)*(10**7)


    

    #Arrays por componentes
    
     Bx1=[]
     By1=[]
     Bz1=[]
     B_perp1=[]
     

    

     for l in range(len(points)):
      Bx1.append(B1[l,0])
      By1.append(B1[l,1]+303)		
      Bz1.append(B1[l,2]-366)
      B_perp1.append(np.sqrt((B1[l,0]*np.sin(Angulo[l])-(B1[l,1]+303)*np.cos(Angulo[l]))**2+(B1[l,2]-366)**2))
      Media.append(np.sqrt((B1[l,0]*np.sin(Angulo[l])-(B1[l,1]+303)*np.cos(Angulo[l]))**2+(B1[l,2]-366)**2))

      


    # make figure with loops in 3D
    
     PMTs_malos = 0
     
     for alfa in range(len(B_perp1)):
      if np.abs(B_perp1[alfa]) > 100:
       PMTs_malos = PMTs_malos+1
       #Coordenadas_paredes.append([points[alfa], Bx1[alfa], By1[alfa], Bz1[alfa], B_perp1[alfa]])

     
    
     print('El numero de PMTs malos es', PMTs_malos)
     
     PMTs_malos_paredes.append(PMTs_malos)
    
     for gamma in range(len(points)):
      Coordenadas_paredes.append([points[gamma], Bx1[gamma], By1[gamma], Bz1[gamma], B_perp1[gamma]])
     
  
    print(np.sum(Longitud1))
    print(np.sum(Longitud2))
    print(np.sum(Longitud3))
    print(np.sum(Longitud1)+np.sum(Longitud2)+np.sum(Longitud3))
    Longitudes = np.sum(Longitud)
    Media_puntos = (np.sum(Media)+ np.sum(Media_top)+np.sum(Media_bottom))/(26602+2*len(PMTs_top))
    print(np.sum(Media_top)/len(PMTs_top), np.sum(Media_bottom)/len(PMTs_top), np.sum(Media)/26602)
    print('La media de Bperp es de', Media_puntos)
    print('La media del campo magnetico en las paredes es de', np.sum(Media)/26602)
    print('La media en las tapas es', (np.sum(Media_top)+np.sum(Media_bottom))/(2*len(PMTs_top)))
    print('El numero de PMTs malos en top es', PMTs_malos_top)
    print('El numero de PMTs malos en bottom es', PMTs_malos_bottom)
    
    PMTs_final = np.sum(PMTs_malos_paredes) + PMTs_malos_top + PMTs_malos_bottom
    print('El numero de PMTs malos total es', PMTs_final)
    print('El numero de PMTs en la pared es', len(z)*len(Angulo), 'en cada una de las tapas', len(PMTs_top), 'y en total en el detector hay', len(z)*len(Angulo)+2*len(PMTs_top))
    print('El porcentaje de PMTs malos es', PMTs_final*100/(len(z)*len(Angulo)+2*len(PMTs_top)))
    print('La longitud de cable que necesitamos es', Longitudes, 'm')
    
    Array_media = np.concatenate((Media, Media_top, Media_bottom))
    RMS = np.sqrt(np.mean(np.square(Array_media)))
    print('La media cuadrática es', RMS)
    
    
    df = pd.DataFrame(Coordenadas_paredes)
    df.to_excel(excel_writer = "/home/sara/Sara/Coordenadasparedes.xlsx")
    
    df = pd.DataFrame(Coordenadas_top)
    df.to_excel(excel_writer = "/home/sara/Sara/Coordenadastop.xlsx")
    
    df = pd.DataFrame(Coordenadas_bottom)
    df.to_excel(excel_writer = "/home/sara/Sara/Coordenadasbottom.xlsx")


    data = np.concatenate((B_perp1, B_perp2, B_perp3))
    

    intervalos = np.arange(0, 210, 10) #calculamos los extremos de los intervalos
    plt.hist(x=data, bins=intervalos, color='#F2AB6D', rwidth=0.85)
    plt.xlabel("Remaining magnetic field perpendicular to PMT (mG)")
    plt.ylabel("Number of PMTs")
    plt.ylim(0,14000)
    plt.xticks(np.arange(0, 200, 25))
    plt.title("$B_{perp}$ distribution for all the PMTs")
    plt.show()
    
    Coord_x=[]
    Coord_y=[]
    Coord_z=[]
    
    x_top=[]
    y_top=[]
    z_top=[]
    
    x_bottom=[]
    y_bottom=[]
    z_bottom=[]
    

    
    
    for indice in range(len(Coordenadas_paredes)):
     Coord_x.append(Coordenadas_paredes[indice][0])	 
     Coord_y.append(Coordenadas_paredes[indice][1])	 
     Coord_z.append(Coordenadas_paredes[indice][2])
     
     
    for indice_top in range(len(Coordenadas_top)):
     x_top.append(Coordenadas_top[indice_top][0])	 
     y_top.append(Coordenadas_top[indice_top][1])	 
     z_top.append(Coordenadas_top[indice_top][2])
     
    
    for indice_bottom in range(len(Coordenadas_bottom)):
     x_bottom.append(Coordenadas_bottom[indice_bottom][0])	 
     y_bottom.append(Coordenadas_bottom[indice_bottom][1])	 
     z_bottom.append(Coordenadas_bottom[indice_bottom][2])	
    

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.scatter(Coord_x, Coord_y, Coord_z,label='PMTs walls')
    ax.scatter(x_bottom, y_bottom, z_bottom, label ='PMTs bottom')
    #plt.legend()
    
    
    
    
    # Cylinder
    
    x_cilindro = np.linspace(-34, 34, 1000)
    z_cilindro = np.linspace(-35.5, 36.5, 1000)
    Xc, Zc=np.meshgrid(x_cilindro, z_cilindro)
    Yc = np.sqrt(34**2-Xc**2)

    rstride = 20
    cstride = 10
    ax.plot_surface(Xc, Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride)
    ax.plot_surface(Xc, -Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.show()

    
    #fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    #ax.scatter(Angulo, B_perp1, label='$\Delta B_{perp}$ h=-33.4')
    #ax.plot(Angulo, limite)
    #ax.plot(Angulo, -limite)
    #plt.xlabel('Angle (rad)')
    #plt.ylabel('B (mG)')
    #plt.xlim(0,35)
    #plt.title('$\Delta B_{perp}$ at walls as a function angle')
    #plt.legend()
    #plt.savefig('Bperp_walls.png')
    #plt.show()
    

    
    return PMTs_vertical, B_perp1



Angulos()
