def SuperK():
    # calculate magnetic fields (in tesla, T)
    # with the Biot-Savart law
    # the output is points (coordinates x,y,z) and B(1,B2,B3) for the 3 cases (B1,B2,B3)
    import numpy as np
    import matplotlib.pyplot as plt
    import wire
    import biotsavart



    

    
    #Definición de constantes
    
    I_rectangular = 150.61
    I_circular = -122.22
    n=7
    I_topbottom = -122.22
    
    Altura = 42
    pos_espira_rectangular =np.array([-20, -14.9, -10.7, -6.4, -2.13, 2.13, 6.4, 10.7, 14.9, 20])
    radio_cilindro = 20
    radio_tapas = 26/2
    pos_espira_circular = np.array([-21, -18.5, -14.5, -10.5, -6.5, -2.5, 1.5, 5.5, 9.5, 13.5, 17.5, 20.4])
    PMTs_vertical = np.arange(-18.4, 18.6, 0.2)
    PMTs_horizontal = np.arange(-17.05, 17.15, 0.1)
    limite = 100*np.ones(len(PMTs_horizontal))
    ceros = np.zeros(len(PMTs_vertical))
    
    #radio_PMT = 
    Angulo1 = 0
    Angulo2 = np.pi/6
    Angulo3 = np.pi/4
    Angulo4 = np.pi/3
    Angulo5 = np.pi/2
    
    Angulo = np.arange(0, np.pi, 0.1)

    
    


    
   #Mallado de puntos
   
    points1=[]
    points2=[]
    



   
 
   
	
    
   #Puntos en el eje Y
   
    for v in range(len(PMTs_vertical)):
     points2.append([0, 0 , PMTs_vertical[v]])
   
    for k in range(len(PMTs_horizontal)):
     points1.append([0, PMTs_horizontal[k] , 0])




    
    

    
    # rectangular loops I=1 A
    
    w1r = wire.Wire(path=wire.Wire.RectangularPath(dx=8.738563955*2, dy=35.35), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,-19.65,0])
    sol = biotsavart.BiotSavart(wire=w1r)
    w11r = wire.Wire(path=wire.Wire.RectangularPath(dx=8.738563955*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,-19.65,0])
    sol.AddWire(w11r)
    
    for i in range(1, len(pos_espira_rectangular)-1):
     lado_espira_rect = np.sqrt(radio_cilindro**2 - pos_espira_rectangular[i]**2)
     w2r = wire.Wire(path=wire.Wire.RectangularPath(dx=lado_espira_rect*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos_espira_rectangular[i],0])
     sol.AddWire(w2r)

    
    w3r = wire.Wire(path=wire.Wire.RectangularPath(dx=8.738563955*2, dy=35.35), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,19.65,0])
    sol.AddWire(w3r)
    w33r = wire.Wire(path=wire.Wire.RectangularPath(dx=8.738563955*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,19.65,0])
    sol.AddWire(w33r)
    
    
    
    # circular loops I=1 A
    
    w1c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas, pts=20), discretization_length=0.1, current=I_topbottom).Translate([0,0, -21])
    #sol = biotsavart.BiotSavart(wire=w1c)
    sol.AddWire(w1c)
    
    for j in range(len(pos_espira_circular)):
     w2c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circular).Translate([0,0,pos_espira_circular[j]])
     sol.AddWire(w2c)

    w3c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_tapas, pts=20), discretization_length=0.1, current=I_topbottom).Translate([0,0, 21])
    sol.AddWire(w3c)
    
    # calculate B field at given points
    B1= sol.CalculateB(points=points1)*(10**7)
    B2= sol.CalculateB(points=points2)*(10**7)
    
    


    


    #Arrays por componentes
    
    Bx1=[]
    By1=[]
    Bz1=[]
    B_perp1=[]
    
    Bx2=[]
    By2=[]
    Bz2=[]
    B_perp2=[]
    
    
    for t in range(len(points2)):
     Bx2.append(B2[t,0])
     By2.append(B2[t,1]+318)		
     Bz2.append(B2[t,2])#-318)
     B_perp2.append(np.sqrt(B2[t,0]**2 + (B2[t,1]+318)**2))

    
    


    
    
    for l in range(len(points1)):
     Bx1.append(B1[l,0])
     By1.append(B1[l,1])#+318)		
     Bz1.append(B1[l,2]-318)
     B_perp1.append(np.sqrt(B1[l,0]**2 + (B1[t,1]+318)**2))
     
     

    



    # make figure with loops in 3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    sol.mpl3d_PlotWires(ax)
    plt.show()
    
    
    #Otras figuras
    
    print(By1)
    
    


    
    fig, ax = plt.subplots()
    #ax.quiver(y, z, By, Bz) # to plot B vectors, not interesting 
    ax.scatter(PMTs_horizontal, By1, label='$\Delta B_{y}$ in y axis')
    #ax.scatter(PMTs_vertical, Bz2, label='$\Delta B_{z}$ in z axis')
    #ax.plot(PMTs_horizontal, limite)
    #ax.plot(PMTs_horizontal, -limite)
    plt.xlabel('Distance (m)')
    plt.ylabel('B (mG)')
    #plt.xlim(0,35)
    plt.title('$\Delta B_{z}$ in its respective axis')
    plt.legend()
    plt.show()
    
    
	
    
    return B2




SuperK()
