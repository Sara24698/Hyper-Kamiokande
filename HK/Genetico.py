import pygad
import numpy as np

solution=[-60, -60, -60] #[Ic, Iv, Itapas, radio tapas1, radio tapas2, distancia espiras]

def Fitness_HK(solution, solution_idx):
    # calculate magnetic fields (in tesla, T)
    # with the Biot-Savart law
    # the output is points (coordinates x,y,z) and B(1,B2,B3) for the 3 cases (B1,B2,B3)
    import numpy as np
    import matplotlib.pyplot as plt
    import wire
    import biotsavart




    

    
    #Definición de constantes
    
    
    #I = Corrientes.Corriente()
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
     
    
    for j in range(5, len(pos_espira_circular)-4):
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
    w22c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=solution[0]).Translate([0,0,pos_espira_circular[3]])
    sol.AddWire(w22c)
    w23c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=solution[1]).Translate([0,0,pos_espira_circular[4]])
    sol.AddWire(w23c)
    w24c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=solution[2]).Translate([0,0, pos_espira_circular[len(pos_espira_circular)-4]])
    sol.AddWire(w24c)





	
    
    
    
    
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
     Media_tapas.append(np.sqrt(B2[q,0]**2+(B2[q,1]+303)**2))
      
     Bx3.append(B3[q,0])
     By3.append(B3[q,1]+303)		
     Bz3.append(B3[q,2]-366)
     B_perp3.append(np.sqrt(B3[q,0]**2+(B3[q,1]+303)**2))
     Media_tapas.append(np.sqrt(B3[q,0]**2+(B3[q,1]+303)**2))
     

    PMTs_malos_top = 0
    PMTs_malos_bottom = 0
    
    
    for beta in range(len(B_perp3)):
      if np.abs(B_perp2[beta]) > 100:
       PMTs_malos_top = PMTs_malos_top+1
       Coordenadas_top.append(PMTs_top[beta])
       

      
      if np.abs(B_perp3[beta]) > 100:
       PMTs_malos_bottom = PMTs_malos_bottom+1
       Coordenadas_bottom.append(PMTs_bottom[beta])

     
    print('El número de PMTs malos en top es', PMTs_malos_top)
    print('El número de PMTs malos en bottom es', PMTs_malos_bottom)
    


   #Mallado de puntos

    for i in range(len(z)):
     points = []
     for j in range(len(Angulo)):
      points.append([radio_PMT*np.cos(Angulo[j]), radio_PMT*np.sin(Angulo[j]), z[i]])
      
     
     	
    

    
    

    
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
     
     w9c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
     w10c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
     w11c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
     w12c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
     w13c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
     w14c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[0]])
     sol.AddWire(w9c)
     sol.AddWire(w10c)
     sol.AddWire(w11c)
     sol.AddWire(w12c)
     sol.AddWire(w13c)
     sol.AddWire(w14c)
     
     w15c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
     w16c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
     w18c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0,pos_espira_circular[len(pos_espira_circular)-1]])
     sol.AddWire(w15c)
     sol.AddWire(w16c)
     sol.AddWire(w18c)
     
    
     for j in range(5, len(pos_espira_circular)-4):
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
     w22c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=solution[0]).Translate([0,0,pos_espira_circular[3]])
     sol.AddWire(w22c)
     w23c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=solution[1]).Translate([0,0,pos_espira_circular[4]])
     sol.AddWire(w23c)
     w24c = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=solution[2]).Translate([0,0, pos_espira_circular[len(pos_espira_circular)-4]])
     sol.AddWire(w24c)
     





     

    
    
    
    
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
       Coordenadas_paredes.append(points[alfa])

     
    
     print('El número de PMTs malos es', PMTs_malos)
     
     PMTs_malos_paredes.append(PMTs_malos)
    
    

    
    Media_puntos = (np.sum(Media)+ np.sum(Media_tapas))/(26602+2*len(PMTs_top))
    print('La media del campo magnético es de', Media_puntos)
    print('La media del campo magnético en las paredes es de', np.sum(Media)/26602)
    print('La media en las tapas es', np.sum(Media_tapas)/(2*len(PMTs_top)))
    print('El número de PMTs malos en top es', PMTs_malos_top)
    print('El número de PMTs malos en bottom es', PMTs_malos_bottom)
    
    PMTs_final = np.sum(PMTs_malos_paredes) + PMTs_malos_top + PMTs_malos_bottom
    print('El número de PMTs malos total es', PMTs_final)
    print('El número de PMTs en la pared es', len(z)*len(Angulo), 'en cada una de las tapas', len(PMTs_top), 'y en total en el detector hay', len(z)*len(Angulo)+2*len(PMTs_top))
    print('El porcentaje de PMTs malos es', PMTs_final*100/(len(z)*len(Angulo)+2*len(PMTs_top)), 'y los parámetros valen', solution[0], solution[1], solution[2])
    Porcentaje = PMTs_final*100/(len(z)*len(Angulo)+2*len(PMTs_top))
    return Porcentaje



fitness_function = Fitness_HK

num_generations = 50
num_parents_mating = 10

sol_per_pop = 20
num_genes = len(solution)

init_range_low = -90
init_range_high = -40

parent_selection_type = "tournament"
keep_parents = 1

crossover_type = "single_point"

mutation_type = "random"
mutation_percent_genes = 10

ga_instance = pygad.GA(num_generations=num_generations,
                       num_parents_mating=num_parents_mating,
                       fitness_func=fitness_function,
                       sol_per_pop=sol_per_pop,
                       num_genes=num_genes,
                       init_range_low=init_range_low,
                       init_range_high=init_range_high,
                       parent_selection_type=parent_selection_type,
                       keep_parents=keep_parents,
                       crossover_type=crossover_type,
                       mutation_type=mutation_type,
                       mutation_percent_genes=mutation_percent_genes,
                       random_mutation_min_val = 0,
                       random_mutation_max_val = 1.0)
                       
                       
ga_instance.run()


ga_instance.plot_result()



solution, solution_fitness, solution_idx = ga_instance.best_solution()
print("Parameters of the best solution : {solution}".format(solution=solution))
print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))

prediction = numpy.sum(numpy.array(function_inputs)*solution)
print("Predicted output based on the best solution : {prediction}".format(prediction=prediction))
