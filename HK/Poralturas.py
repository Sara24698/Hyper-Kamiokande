
def Angulos():
    # calculate magnetic fields (in tesla, T)
    # with the Biot-Savart law
    # the output is points (coordinates x,y,z) and B(1,B2,B3) for the 3 cases (B1,B2,B3)
    import numpy as np
    import matplotlib.pyplot as plt
    import wire
    import biotsavart
    import Corrientes



    

    
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
    limite = 100*np.ones(len(Angulo))
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
    
    x = np.arange(-31.95, 32, 0.71)
    y = np.arange(-31.95, 32, 0.71)
    
    #for r in range(len(x)):
     #for h in range(len(y)):
      #points2.append([x[r], y[h]])
    
    
    #for g in range(len(points2)):
     #distancia = np.sqrt(points2[g][0]**2+points2[g][1]**2)
     #if distancia <=31.95:
      #PMTs_top.append([points2[g][0], points2[g][1], 32.351])
      #PMTs_bottom.append([points2[g][0], points2[g][1], -33.5])
    

    for r in range(len(x)):
     for h in range(len(y)):
      points2.append([x[r], y[h]])
    
    
    for g in range(len(points2)):
     distancia = np.sqrt(points2[g][0]**2+points2[g][1]**2)
     if distancia <=32.01:
      PMTs_top.append([points2[g][0], points2[g][1], 32.9])
    
    #for i in range(len(z)):
    #points = []
    #for j in range(len(Angulo)):
     #points.append([radio_PMT*np.cos(Angulo[j]), radio_PMT*np.sin(Angulo[j]), 0])
		
    

    
    

    
    # rectangular loops I=1 A


    w1r = wire.Wire(path=wire.Wire.RectangularPath(dx=np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*2, dy=Altura), discretization_length=0.1, current=1.3*I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[0], 0])
    sol = biotsavart.BiotSavart(wire=w1r)
   
    
    for i in range(1, len(pos_espira_rectangular)-1):
     lado_espira_rect = np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)
     w2r = wire.Wire(path=wire.Wire.RectangularPath(dx=lado_espira_rect*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[i],0])
     sol.AddWire(w2r)
     
    w3r = wire.Wire(path=wire.Wire.RectangularPath(dx=np.sqrt(radio_cilindro**2 - pos_espira_rectangular[0]**2)*2, dy=Altura), discretization_length=0.1, current=1.4*I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[-1], 0])
    sol.AddWire(w3r)
    



    # circular loops I=1 A
    
    
    w1c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0, pos_espira_circular[0]])
    w2c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0, pos_espira_circular[0]])
    w25c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0, pos_espira_circular[0]])
    w3c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas2, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0, pos_espira_circular[0]])
    w4c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas2, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0, pos_espira_circular[0]])
    #sol = biotsavart.BiotSavart(wire=w1c)
    sol.AddWire(w1c)
    sol.AddWire(w2c)
    sol.AddWire(w25c)
    sol.AddWire(w3c)
    sol.AddWire(w4c)
     
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
     
    w15c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
    w16c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
    w18c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
    sol.AddWire(w15c)
    sol.AddWire(w16c)
    sol.AddWire(w18c)
     
    
    for j in range(4, len(pos_espira_circular)-3):
     w17c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[j]])
     sol.AddWire(w17c)

    w7c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas2, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0, pos_espira_circular[len(pos_espira_circular)-1]])
    w8c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas2, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0, pos_espira_circular[len(pos_espira_circular)-1]])
    sol.AddWire(w7c)
    sol.AddWire(w8c)

    
    #w1e = wire.Wire(path=wire.Wire.EllipticalPath(rx=99.04/2, ry=17, pts=20), discretization_length=0.1, current=I_circulares).Rotate(axis=(0,1,0), deg=225).Rotate(axis=(0,0,1), deg=270)
    #sol.AddWire(w1e)
    
    w20c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[1]])
    sol.AddWire(w20c)
    w21c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-67).Translate([0,0, pos_espira_circular[len(pos_espira_circular)-3]])
    sol.AddWire(w21c)
    w22c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=-80).Translate([0,0,pos_espira_circular[3]])
    sol.AddWire(w22c)


    
    
    
    # calculate B field at given points
    B1 = sol.CalculateB(points=PMTs_top)*(10**7)

    

    #Arrays por componentes
    
    Bx1=[]
    By1=[]
    Bz1=[]
    B_perp1=[]
    

 
    
  
    

    for l in range(len(PMTs_top)):
     Bx1.append(B1[l,0])
     By1.append(B1[l,1]+303)		
     Bz1.append(B1[l,2]-366)
     B_perp1.append(np.sqrt(B1[l,0]**2+(B1[l,1]+303)**2))
     #Media.append(np.sqrt((B1[l,0]*np.sin(Angulo[l])-(B1[l,1]+303)*np.cos(Angulo[l]))**2+(B1[l,2]-366)**2))
     





    # make figure with loops in 3D
    
    PMTs_malos = 0
    Coordenadas_paredes=[]
    
    for alfa in range(len(B_perp1)):
     if np.abs(B_perp1[alfa]) > 100:
      PMTs_malos = PMTs_malos+1
      Coordenadas_paredes.append([points[alfa]])
    
    print('El número de PMTs malos es', PMTs_malos)
    print(Coordenadas_paredes)
     
    PMTs_malos_total.append(PMTs_malos)
    
    
    Media_puntos = np.sum(Media)/len(points)
    print('La media del campo magnético es de', Media_puntos)
     
    PMTs_final = np.sum(PMTs_malos_total)
    print('El número de PMTs malos total es', PMTs_final)
    
    
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    #ax.scatter(Angulo, Bx1, label='$B_{x}$')
    #ax.scatter(Angulo, By1, label='$\Delta B_{y}$')
    #ax.scatter(Angulo, Bz1, label='$\Delta B_{z}$')
    ax.scatter(Angulo, B_perp1, label='$\Delta B_{perp}$')
    ax.plot(Angulo, limite)
    #ax.plot(Angulo, -limite)
    plt.xlabel('Angle (rad)')
    plt.ylabel('B (mG)')
    #plt.xlim(0,35)
    plt.title('$B_{perp}$ at walls on bottom as a function of angle')
    #plt.legend()
    #plt.savefig('Bperp_walls.png')
    plt.show()
    

    
    return PMTs_vertical, B_perp1




Angulos()
