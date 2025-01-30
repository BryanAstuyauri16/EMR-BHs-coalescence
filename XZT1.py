import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import os
import glob
import pandas as pd

def Plot_XZT(m, varphi):
    colors = ['Black', 'DarkBlue', 'DarkGreen', 'DarkRed', 'Purple', 'Teal', 'Gold', 'DarkOrange', 'SlateGray', 'Crimson']
    labels = []
    for i in range(11):
        labels.append(f'q = {i}')

    # directory = f"RN_Data/BH{varphi}/phi0/"   
    # directory = f"RN_Data_2.5/BH{varphi}/phi0/"   
    directory = f"Sch5/BH{varphi}/phi0/"   

    file_pattern = os.path.join(directory, "data*.txt")
    files = glob.glob(file_pattern)
    for item, file in enumerate(files): 
        P_alpha = float(file[len(directory)+4:][:-4]) 
        try:
            file_data = np.loadtxt(file)
            x = file_data[:, 0]
            y = file_data[:, 1]
            z = file_data[:, 2] 
            t = file_data[:, 3]

            # mask = (z < 0)
            # x = x[mask]
            # y = y[mask]
            # z = z[mask]
            # t = t[mask]

            q_caustic = 2.2286 * (2*m)

            if (abs(P_alpha) == q_caustic):
                ax.plot(x, z, t, color= 'DarkOrange', markersize = 1, linewidth = 1, linestyle = '--')
            else:
                ax.plot(x, z, t, color= colors[varphi], markersize = 1, linewidth = 1, alpha = 1)
            
        except Exception as e:
            print(f"No se pudo leer {file}: {e}")

    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))  # Eliminar el fondo de la pared en el eje X
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))  # Eliminar el fondo de la pared en el eje Y
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))  # Eliminar el fondo de la pared en el eje Z

    ax.set_xlabel('X')
    ax.set_ylabel('Z')
    ax.set_zlabel('T')
    ax.set_box_aspect(None, zoom=.8)
    ax.set_axis_off()

    sh = m * 20
    ax.set_zlim(0, sh)
    ax.set_ylim(-sh, sh)
    ax.set_xlim(-sh, sh)
    ax.plot([-sh*5, sh*5], [0, 0], [0, 0], alpha = .3, color='black', linewidth=.5)  
    ax.plot([0, 0], [-sh*5, sh*5], [0, 0], alpha = .3, color='black', linewidth=.5)  
    ax.plot([0, 0], [0, 0], [-sh*5, sh*5], alpha = .3, color='black', linewidth=.5)  

    foldername = 'Figures'
    path = f"{foldername}/PLOT 3D XZT"
    directory = os.path.dirname(path)

    if not os.path.exists(directory):
        os.makedirs(directory)

    ax.view_init(elev=27, azim=-26)  # elev es la altura, azim es el ángulo horizontal
    # plt.savefig(f'{path} POV 1.pdf', format='pdf', bbox_inches='tight')


m = 50

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
Plot_XZT(m, 0)

plt.show()

