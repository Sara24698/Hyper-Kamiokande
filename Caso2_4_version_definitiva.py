# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:14:05 2024

@author: ICTEA
"""

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


I_rectangular = 69
I_circulares = -82


Altura = 72
radio_cilindro = 34
pos_espira_rectangular = np.arange(-32,33,2)
pos_espira_circular = np.arange(-33.1,34.2,2.4)
PMTs_vertical = np.arange(-33.4,32.5,0.707)
radio_PMT= 32.4



#Programa principal

def create_espiras():
    
    sol = biotsavart.BiotSavart()

	# rectangular loops I=1 A
    
    for pos in pos_espira_rectangular:
        lado_espira_rect = np.sqrt(radio_cilindro**2 - pos**2)
        wr = wire.Wire(path=wire.Wire.RectangularPath(dx=lado_espira_rect*2, dy=Altura), discretization_length=0.1, current=I_rectangular).Rotate(axis=(1, 0, 0), deg=-90).Translate([0,pos,0]).Translate([0,0, 0.5])
        sol.AddWire(wr)

	

	
	
	# circular loops I=1 A



    wc_grande_top = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=5*I_circulares).Translate([0,0,36.5])
    sol.AddWire(wc_grande_top)
	
	
    wc_pequeña_bottom = wire.Wire(path=wire.Wire.CircularPath(radius=25, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0, -35.5])
    sol.AddWire(wc_pequeña_bottom)

    wc_pequeña_top = wire.Wire(path=wire.Wire.CircularPath(radius=26, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0, 36.5])
    sol.AddWire(wc_pequeña_top)


    wc_grande_bottom = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=5*I_circulares).Translate([0,0,-35.5])
    sol.AddWire(wc_grande_bottom)




    for pos in pos_espira_circular:     
        wc = wire.Wire(path=wire.Wire.CircularPath(radius=radio_cilindro, pts=20), discretization_length=0.1, current=I_circulares).Translate([0,0, pos])
        sol.AddWire(wc)
        
    return sol

def generate_points():
    
    x_grid = np.arange(-31.815, 32, 0.707)
    y_grid = np.arange(-31.815, 32, 0.707)
    
    PMTs_top = [[x, y, 32.9] for x in x_grid for y in y_grid if np.sqrt(x**2 + y**2) <= 32.01]
    PMTs_bottom = [[coords[0], coords[1], -32.9] for coords in PMTs_top]

            
            
    Angulo = np.arange(0, 6.315, 0.0218)
    z=np.arange(-32.522, 32.523, 0.707)
    
    PMTs_paredes = []
    Angulos_z = []
    for i in range(len(z)):
        for j in range(len(Angulo)):
            PMTs_paredes.append([radio_PMT*np.cos(Angulo[j]), radio_PMT*np.sin(Angulo[j]), z[i]])
            Angulos_z.append(Angulo[j])
	
	
    return PMTs_top, PMTs_bottom, PMTs_paredes, Angulos_z


def calculate_magnetic_field(sol, PMTs_top, PMTs_bottom, PMTs_paredes, Angulos_z):

    
    B_top = sol.CalculateB(points=PMTs_top)*10**7
    B_bottom = sol.CalculateB(points=PMTs_bottom)*10**7
    B_paredes = sol.CalculateB(points=PMTs_paredes)*10**7


    Bx_top = B_top[:, 0]
    By_top = B_top[:, 1]+303
    Bz_top = B_top[:, 2]-366
    B_perp_top = np.sqrt(Bx_top**2 + By_top**2)
    
    Bx_bottom = B_bottom[:, 0]
    By_bottom = B_bottom[:, 1]+303
    Bz_bottom = B_bottom[:, 2]-366
    B_perp_bottom = np.sqrt(Bx_bottom**2 + By_bottom**2)
    
    Bx_paredes = B_paredes[:, 0]
    By_paredes = B_paredes[:, 1]+303
    Bz_paredes = B_paredes[:, 2]-366
    
    B_perp_paredes=[]
    for i in range(len(PMTs_paredes)):
        B_perp_paredes.append(np.sqrt((Bx_paredes[i]*np.sin(Angulos_z[i])-By_paredes[i]*np.cos(Angulos_z[i]))**2+Bz_paredes[i]**2))

    return B_perp_top, B_perp_bottom, B_perp_paredes, B_top, B_bottom, B_paredes

		
def export_data (PMTs_top, PMTs_bottom, PMTs_paredes, B_perp_top, B_perp_bottom, B_perp_paredes, B_top, B_bottom, B_paredes):	
    
    # make figure with loops in 3D
    Coordenadas=[]
    B_all = np.concatenate([B_top, B_bottom, B_paredes])
    B_perp_all = np.concatenate([B_perp_top, B_perp_bottom, B_perp_paredes])
    PMTs_all = np.concatenate([PMTs_top, PMTs_bottom, PMTs_paredes])
    

    for alfa in range(len(PMTs_all)):
        Coordenadas.append([PMTs_all[alfa][0], PMTs_all[alfa][1], PMTs_all[alfa][2], B_all[alfa][0], B_all[alfa][1], B_all[alfa][2], B_perp_all[alfa], 1 if PMTs_all[alfa][2] == 32.9 else 3 if PMTs_all[alfa][2] == -32.9 else 2])
			

    df = pd.DataFrame(Coordenadas, columns=['x', 'y', 'z', 'Bx', 'By', 'Bz', 'Bp', 'faceid'])

    # Export to CSV
    df.to_csv('pmt_data_above_100mg.csv', index=False)

    print("Data exported successfully to 'pmt_data_above_100mg.csv'.")
	
	


def results(B_perp_all):
    mean_B_perp = np.mean(B_perp_all)
    std_B_perp = np.std(B_perp_all)

    # Calculate proportion of PMTs with B_perp >= 100 mG
    threshold = 100  # mG
    proportion_above_threshold = np.sum(B_perp_all >= threshold) / len(B_perp_all)

    # Display results
    #total_pmt = len(PMTs_top) + len(PMTs_bottom) + len(PMTs_paredes)

# Print the total number of PMTs
    #print(f"Total number of PMTs being simulated: {total_pmt}")
    print(f"Mean B_perp: {mean_B_perp:.2f} mG")
    print(f"Standard Deviation of B_perp: {std_B_perp:.2f} mG")
    print(f"Proportion of PMTs with B_perp >= {threshold} mG: {proportion_above_threshold:.2%}")

    # Histogram of B_perp
    plt.hist(B_perp_all, bins=20, color='orange', alpha=0.7)
    plt.title("Distribution of $B_{\\perp}$")
    plt.xlabel("Magnetic Field Perpendicular to PMT (mG)")
    plt.ylabel("Number of PMTs")
    legend_text = (f"Mean: {mean_B_perp:.2f} mG\n"
               f"Std Dev: {std_B_perp:.2f} mG\n"
               f"PMTs > 100 mG: {proportion_above_threshold*100:.2f}%")
    plt.legend([legend_text], loc='upper right', fontsize=10, frameon=True)


    plt.tight_layout()
    plt.show()
    
    
def cylinder(B_perp_all, PMTs_all):


    threshold = 100  # mG

    # Filter the data to keep only points with B_perp > 100 mG
    mask = B_perp_all > threshold

    # Apply the mask to filter the coordinates and B_perp values
    x = PMTs_all[mask, 0]
    y = PMTs_all[mask, 1]
    z = PMTs_all[mask, 2]
    B_perp_filtered = B_perp_all[mask]
    
    # Creating the 3D scatter plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Scatter plot with color mapping for B_perp
    sc = ax.scatter(x, y, z, c=B_perp_filtered, cmap='viridis')  # You can change 'viridis' to other color maps
    
    # Adding color bar
    plt.colorbar(sc, label='B_perp')
    
    # Adding labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('PMTs with Bperp > 100 mG')
    
    # Show the plot
    plt.show()
    

if __name__ == "__main__":
    # Generate PMT coordinates for top, bottom, and walls
    PMTs_top, PMTs_bottom, PMTs_paredes, Angulos_z = generate_points()

    # Create Biot-Savart solution
    sol = create_espiras()

    # Calculate magnetic fields
    B_perp_top, B_perp_bottom, B_perp_paredes, B_top, B_bottom, B_paredes = calculate_magnetic_field(sol, PMTs_top, PMTs_bottom, PMTs_paredes, Angulos_z)

    # Combine B_perp values for all PMTs
    B_perp_all = np.concatenate([B_perp_top, B_perp_bottom, B_perp_paredes])
    
    PMTs_all = np.concatenate([PMTs_top, PMTs_bottom, PMTs_paredes])

    #Calculate the results
    results(B_perp_all)
    
    cylinder(B_perp_all, PMTs_all)
    
