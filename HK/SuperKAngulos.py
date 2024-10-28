def SuperKAngulos():
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
    limite = 100*np.ones(len(PMTs_vertical))
    ceros = np.zeros(len(PMTs_vertical))
    
    
    radios_PMT = np.arange(0, 17.15, 0.1)
    radio_PMT = 17.05
    Angulo1 = 0
    Angulo2 = np.pi/6
    Angulo3 = np.pi/4
    Angulo4 = np.pi/3
    Angulo5 = np.pi/2
    
    Angulo = np.arange(0, np.pi, 0.1)
    


    
   #Mallado de puntos
   
    points1=[]
    points2=[]
    points3=[]
    points4=[]
    points5=[]
    
    points_tapas1=[]
    points_tapas2=[]
    points_tapas3=[]
    points_tapas4=[]
    points_tapas5=[]

    

 
   
	
    
   #Puntos en el eje Y
   
    for k in range(len(PMTs_vertical)):
     points1.append([radio_PMT*np.cos(Angulo1), radio_PMT*np.sin(Angulo1) , PMTs_vertical[k]])
     points2.append([radio_PMT*np.cos(Angulo2), radio_PMT*np.sin(Angulo2) , PMTs_vertical[k]])
     points3.append([radio_PMT*np.cos(Angulo3), radio_PMT*np.sin(Angulo3) , PMTs_vertical[k]])
     points4.append([radio_PMT*np.cos(Angulo4), radio_PMT*np.sin(Angulo4) , PMTs_vertical[k]])
     points5.append([radio_PMT*np.cos(Angulo5), radio_PMT*np.sin(Angulo5) , PMTs_vertical[k]])
   
    for r in range(len(radios_PMT)):
     points_tapas1.append([radios_PMT[r]*np.cos(Angulo1), radios_PMT[r]*np.sin(Angulo1) , 18.4])
     points_tapas2.append([radios_PMT[r]*np.cos(Angulo2), radios_PMT[r]*np.sin(Angulo2) , 18.4])
     points_tapas3.append([radios_PMT[r]*np.cos(Angulo3), radios_PMT[r]*np.sin(Angulo3) , 18.4])
     points_tapas4.append([radios_PMT[r]*np.cos(Angulo4), radios_PMT[r]*np.sin(Angulo4) , 18.4])
     points_tapas5.append([radios_PMT[r]*np.cos(Angulo5), radios_PMT[r]*np.sin(Angulo5) , 18.4])

    
    

    
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
    B1 = sol.CalculateB(points=points1)*(10**7)
    B2 = sol.CalculateB(points=points2)*(10**7)
    B3 = sol.CalculateB(points=points3)*(10**7)
    B4 = sol.CalculateB(points=points4)*(10**7)
    B5 = sol.CalculateB(points=points5)*(10**7)
    
    B_tapas1 = sol.CalculateB(points=points_tapas1)*(10**7)
    B_tapas2 = sol.CalculateB(points=points_tapas2)*(10**7)
    B_tapas3 = sol.CalculateB(points=points_tapas3)*(10**7)
    B_tapas4 = sol.CalculateB(points=points_tapas4)*(10**7)
    B_tapas5 = sol.CalculateB(points=points_tapas5)*(10**7)
    
   

    


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
    
    B_tapasx1=[]
    B_tapasy1=[]
    B_tapasz1=[]
    B_tapas_perp1=[]
    
    
    B_tapasx2=[]
    B_tapasy2=[]
    B_tapasz2=[]
    B_tapas_perp2=[]
    
    
    B_tapasx3=[]
    B_tapasy3=[]
    B_tapasz3=[]
    B_tapas_perp3=[]
    
    
    B_tapasx4=[]
    B_tapasy4=[]
    B_tapasz4=[]
    B_tapas_perp4=[]
   
    
    B_tapasx5=[]
    B_tapasy5=[]
    B_tapasz5=[]
    B_tapas_perp5=[]
    


    
    
    for l in range(len(points1)):
     Bx1.append(B1[l,0])
     By1.append(B1[l,1]+318)		
     Bz1.append(B1[l,2]-318)
     B_perp1.append(np.sqrt((B1[l,0]*np.sin(Angulo1)-(B1[l,1]+318)*np.cos(Angulo1))**2+(B1[l,2]-318)**2))
     
     Bx2.append(B2[l,0])
     By2.append(B2[l,1]+318)		
     Bz2.append(B2[l,2]-366)
     B_perp2.append(np.sqrt((B2[l,0]*np.sin(Angulo2)-(B2[l,1]+318)*np.cos(Angulo2))**2+(B2[l,2]-318)**2))
     
     Bx3.append(B3[l,0])
     By3.append(B3[l,1]+318)		
     Bz3.append(B3[l,2]-318)
     B_perp3.append(np.sqrt((B3[l,0]*np.sin(Angulo3)-(B3[l,1]+318)*np.cos(Angulo3))**2+(B3[l,2]-318)**2))
     
     Bx4.append(B4[l,0])
     By4.append(B4[l,1]+318)		
     Bz4.append(B4[l,2]-318)
     B_perp4.append(np.sqrt((B4[l,0]*np.sin(Angulo4)-(B4[l,1]+318)*np.cos(Angulo4))**2+(B4[l,2]-318)**2))
     
     Bx5.append(B5[l,0])
     By5.append(B5[l,1]+318)		
     Bz5.append(B5[l,2]-318)
     B_perp5.append(np.sqrt((B5[l,0]*np.sin(Angulo5)-(B5[l,1]+318)*np.cos(Angulo5))**2+(B5[l,2]-318)**2))
     

    for s in range(len(points_tapas1)):
     B_tapasx1.append(B_tapas1[s,0])
     B_tapasy1.append(B_tapas1[s,1]+318)		
     B_tapasz1.append(B_tapas1[s,2]-318)
     B_tapas_perp1.append(np.sqrt(B_tapas1[s,0]**2+(B_tapas1[s,1]+318)**2))


     B_tapasx2.append(B_tapas2[s,0])
     B_tapasy2.append(B_tapas2[s,1]+318)		
     B_tapasz2.append(B_tapas2[s,2]-318)
     B_tapas_perp2.append(np.sqrt(B_tapas2[s,0]**2+(B_tapas2[s,1]+318)**2))


     B_tapasx3.append(B_tapas3[s,0])
     B_tapasy3.append(B_tapas3[s,1]+318)		
     B_tapasz3.append(B_tapas3[s,2]-318)
     B_tapas_perp3.append(np.sqrt(B_tapas3[s,0]**2+(B_tapas3[s,1]+318)**2))


     B_tapasx4.append(B_tapas4[s,0])
     B_tapasy4.append(B_tapas4[s,1]+318)		
     B_tapasz4.append(B_tapas4[s,2]-318)
     B_tapas_perp4.append(np.sqrt(B_tapas4[s,0]**2+(B_tapas4[s,1]+318)**2))


     B_tapasx5.append(B_tapas5[s,0])
     B_tapasy5.append(B_tapas5[s,1]+318)		
     B_tapasz5.append(B_tapas5[s,2]-318)
     B_tapas_perp5.append(np.sqrt(B_tapas5[s,0]**2+(B_tapas5[s,1]+318)**2))
     





    # make figure with loops in 3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    sol.mpl3d_PlotWires(ax)
    plt.show()
    
    
    
    
    #Otras figuras
    

    
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
    plt.show()
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(PMTs_vertical, Bz1, label='$\Delta B_{z}$ 0º')
    ax.scatter(PMTs_vertical, Bz2, label='$\Delta B_{z}$ 30º')
    ax.scatter(PMTs_vertical, Bz3, label='$\Delta B_{z}$ 45º')
    ax.scatter(PMTs_vertical, Bz4, label='$\Delta B_{z}$ 60º')
    ax.scatter(PMTs_vertical, Bz5, label='$\Delta B_{z}$ 90º')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Height (m)')
    plt.ylabel('B (mG)')
    #plt.xlim(0,35)
    plt.title('$\Delta B_{z}$ at walls as a function of height and angle')
    plt.legend()
    plt.show()
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(PMTs_vertical, Bx1, label='$B_{x}$ 0º')
    ax.scatter(PMTs_vertical, Bx2, label='$B_{x}$ 30º')
    ax.scatter(PMTs_vertical, Bx3, label='$B_{x}$ 45º')
    ax.scatter(PMTs_vertical, Bx4, label='$B_{x}$ 60º')
    ax.scatter(PMTs_vertical, Bx5, label='$B_{x}$ 90º')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Height (m)')
    plt.ylabel('B (mG)')
    #plt.xlim(0,35)
    plt.title('$B_{x}$ at walls as a function of height and angle')
    plt.legend()
    plt.show()
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(PMTs_vertical, By1, label='$\Delta B_{y}$ 0º')
    ax.scatter(PMTs_vertical, By2, label='$\Delta B_{y}$ 30º')
    ax.scatter(PMTs_vertical, By3, label='$\Delta B_{y}$ 45º')
    ax.scatter(PMTs_vertical, By4, label='$\Delta B_{y}$ 60º')
    ax.scatter(PMTs_vertical, By5, label='$\Delta B_{y}$ 90º')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Height (m)')
    plt.ylabel('B (mG)')
    #plt.xlim(0,35)
    plt.title('$\Delta B_{y}$ at walls as a function of height and angle')
    plt.legend()
    plt.show()
    
    
    
    #Figuras en las tapas
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(radios_PMT, B_tapas_perp1, label='$\Delta B_{perp}$ $0 \: rad$')
    ax.scatter(radios_PMT, B_tapas_perp2, label='$\Delta B_{perp}$ $\pi/6 \: rad$')
    ax.scatter(radios_PMT, B_tapas_perp3, label='$\Delta B_{perp}$ $\pi/4 \: rad$')
    ax.scatter(radios_PMT, B_tapas_perp4, label='$\Delta B_{perp}$ $\pi/3 \: rad$')
    ax.scatter(radios_PMT, B_tapas_perp5, label='$\Delta B_{perp}$ $\pi/2 \: rad$')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Radius (m)')
    plt.ylabel('B (mG)')
    plt.xlim(0,18)
    plt.title('$\Delta B_{perp}$ at top as a function of radius and angle')
    plt.legend()
    plt.show()
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(radios_PMT, B_tapasz1, label='$\Delta B_{z}$ 0º')
    ax.scatter(radios_PMT, B_tapasz2, label='$\Delta B_{z}$ 30º')
    ax.scatter(radios_PMT, B_tapasz3, label='$\Delta B_{z}$ 45º')
    ax.scatter(radios_PMT, B_tapasz4, label='$\Delta B_{z}$ 60º')
    ax.scatter(radios_PMT, B_tapasz5, label='$\Delta B_{z}$ 90º')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Radius (m)')
    plt.ylabel('B (mG)')
    plt.xlim(0,18)
    plt.title('$\Delta B_{z}$ at bottom as a function of radius and angle')
    plt.legend()
    plt.show()
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(radios_PMT, B_tapasy1, label='$\Delta B_{y}$ 0º')
    ax.scatter(radios_PMT, B_tapasy2, label='$\Delta B_{y}$ 30º')
    ax.scatter(radios_PMT, B_tapasy3, label='$\Delta B_{y}$ 45º')
    ax.scatter(radios_PMT, B_tapasy4, label='$\Delta B_{y}$ 60º')
    ax.scatter(radios_PMT, B_tapasy5, label='$\Delta B_{y}$ 90º')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Radius (m)')
    plt.ylabel('B (mG)')
    plt.xlim(0,18)
    plt.title('$\Delta B_{y}$ at bottom as a function of radius and angle')
    plt.legend()
    plt.show()
    
    fig, ax = plt.subplots()
    #ax.quiver(x, y, Bx, By) # to plot B vectors, not interesting 
    ax.scatter(radios_PMT, B_tapasx1, label='$B_{x}$ 0º')
    ax.scatter(radios_PMT, B_tapasx2, label='$B_{x}$ 30º')
    ax.scatter(radios_PMT, B_tapasx3, label='$B_{x}$ 45º')
    ax.scatter(radios_PMT, B_tapasx4, label='$B_{x}$ 60º')
    ax.scatter(radios_PMT, B_tapasx5, label='$B_{x}$ 90º')
    ax.plot(PMTs_vertical, limite)
    ax.plot(PMTs_vertical, -limite)
    plt.xlabel('Radius (m)')
    plt.ylabel('B (mG)')
    plt.xlim(0,18)
    plt.title('$B_{x}$ at bottom as a function of radius and angle')
    plt.legend()
    plt.show()

	
    
    return B1

SuperKAngulos()
