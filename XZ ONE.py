import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from matplotlib.widgets import Slider, TextBox, Button
import pandas as pd
from scipy.optimize import curve_fit


def Black_water_drop(m, varphi, S):
    data = []
    try:
        directory = f"Sch5/BH{varphi}/phi0/"   
        file_pattern = os.path.join(directory, "data*.txt")
        files = glob.glob(file_pattern)

        for item, file in enumerate(files):    
            P_alpha = float(file[len(directory)+4:][:-4]) 
            try:
                file_data = np.loadtxt(file)
                x = file_data[:, 0]
                z = file_data[:, 2]
                t = file_data[:, 3]

                dt = m / (200)
                mask = (t > S - dt) & (t < S + dt)
                x_filtrado = x[mask]
                z_filtrado = z[mask]
                t_filtrado = t[mask]

                if abs(P_alpha) <= 2.67848 * (2 * m) :
                # if x.size > 0 :
                    q_df = pd.DataFrame({
                        'X': x_filtrado,
                        'Z': z_filtrado,
                        'T': t_filtrado,
                    })
                    q_df_mean = q_df[['X', 'Z', 'T']].mean().to_frame().T  
                    q_df_mean['q'] = varphi
                    q_df_mean['P_alpha'] = P_alpha

                    data.append(q_df_mean)

            except Exception as e:
                print(f"No se pudo leer {file}: {e}")     
                # pass
    except: 
        pass

    # data = [df for df in data if not df.empty]
    All_Data = pd.concat(data, ignore_index=True)
    z_min, z_max = All_Data['Z'].min(), All_Data['Z'].max() 
    return All_Data

def Select_Model(which):
    def Model1(theta, a, b, c, d):
        return a * np.abs(np.sin(b * theta))**c + d
    def Model2(theta, a, b, c, d):
        return a * (1- b*theta**2)**c + d
    def Model3(theta, a, b, c):
        return  a * np.e**(-b * theta**2) + c
    def Model4(theta, a, b, c, d):
        return a * np.cosh(b * theta)**(-c) + d
    def Model5(theta, a, b, c, d):
        return a / (1 + b * theta**2)**c + d
    models = {1: Model1, 2: Model2, 3: Model3, 4: Model4, 5: Model5}

    IC = {1: [1, 2, 2, 0], 2: [-1.125, 10, 1, 10], 3: [2, 25, 100], 4: [34, -15, .1, 82], 5: [10, 100, 100]}
    
    def Derivative4(theta, a, b, c, d):
        def drdtheta(theta, a, b, c):
            return -a*b*c*np.cosh(b*theta)*np.sinh(b*theta)
        dz = -drdtheta(theta, a, b, c)*np.cos(theta) + Model4(theta, a, b, c, d)*np.sin(theta)
        dx = drdtheta(theta, a, b, c)*np.sin(theta) + Model4(theta, a, b, c, d)*np.cos(theta)
        return dz / dx 
    Derivatives = {1: Model1, 2: Model2, 3: Model3, 4: Derivative4, 5: Model5}
    
    def DDerivative4(theta, a, b, c, d):
        def drdtheta(theta, a, b, c):
            return -a*b*c*(np.cosh(b*theta)**(-c-1))*np.sinh(b*theta)
        def ddrdtheta(theta, a, b, c):
            return a*(b**2)*c*(c+1)*(np.cosh(b*theta)**(-c-2))*np.sinh(b*theta)**2 - a*(b**2)*c*(np.cosh(b*theta)**(-c))
        dz = -drdtheta(theta, a, b, c)*np.cos(theta) + Model4(theta, a, b, c, d)*np.sin(theta)
        dx = drdtheta(theta, a, b, c)*np.sin(theta) + Model4(theta, a, b, c, d)*np.cos(theta)

        ddz = ((-ddrdtheta(theta, a, b, c)*np.cos(theta) + 2*drdtheta(theta, a, b, c)*(np.sin(theta)) + Model4(theta, a, b, c, d)*np.cos(theta))*dx - (ddrdtheta(theta, a, b, c)*np.sin(theta) + 2*drdtheta(theta, a, b, c)*np.cos(theta) - Model4(theta, a, b, c, d)*np.sin(theta))*dz) / (dx**2)
        return ddz / dx 
    DDerivatives = {1: Model1, 2: Model2, 3: Model3, 4: DDerivative4, 5: Model5}

    return models.get(which, None), IC.get(which, None), Derivatives.get(which, None), DDerivatives.get(which, None)

def Fit_Model(m, varphi, Data, z_max):
    s = 4
    r_data = np.sqrt(Data['X']**2 + Data['Z']**2)
    if s == 1:
        theta_data = np.arctan2(Data['Z'], Data['X'])
    elif s == 2:
        theta_data = np.arctan2(Data['Z'], Data['X']) - np.pi
    elif s in [3,4]:
        theta_data = np.arctan2(Data['Z'], Data['X']) + np.pi/2

    Model, p0, Der, DDer = Select_Model(s)
    params, covarianza = curve_fit(Model, theta_data, r_data, p0=p0, maxfev=5000)
    theta_fit = np.linspace(-np.pi, np.pi, 10000)
    r_fit = Model(theta_fit, *params)
    # r_fit = Model(theta_fit, p0[0], p0[1], p0[2],p0[3])
    # X_fit = r_data * np.cos(theta_data)
    # Z_fit = r_data * np.sin(theta_data)
    # X_fit = r_fit * np.cos(theta_fit)
    # Z_fit = r_fit * np.sin(theta_fit)
    if s == 1:
        X_fit = r_fit * np.cos(theta_fit)
        Z_fit = r_fit * np.sin(theta_fit)
    if s == 2:
        Z_fit = r_fit * np.cos(theta_fit)
        X_fit = -r_fit * np.sin(theta_fit)
    elif s in [3,4]:
        Z_fit = -r_fit * np.cos(theta_fit)
        X_fit = r_fit * np.sin(theta_fit)

    mask = Z_fit < z_max  
    X_fit = X_fit[mask]   
    Z_fit = Z_fit[mask]  

    Dz_dx0 = Der(0, *params)
    DDz_dx0 = DDer(0, *params)
    return X_fit, Z_fit, Dz_dx0, DDz_dx0

def calcular_area_dataframe(Data):
    z_min, z_max = Data['Z'].min(), Data['Z'].max() 
    x_min, x_max = Data['X'].min(), Data['X'].max() 

    print(z_min, z_max, x_min, x_max)
    x = Data['X'].values
    y = Data['Z'].values
    # Sshoelace formula
    area = 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
    return area

def Plot_XZ_var(m, varphi, z_max):
    colors = ['Black', 'DarkBlue', 'DarkGreen', 'DarkRed', 'Purple', 'Teal', 'Gold', 'DarkOrange', 'SlateGray', 'Crimson']

    T_init = -5* (2*m)

    Data = Black_water_drop(m, varphi, T_init)
    Data = Data[(Data['Z'] < z_max)]
    print(Data)

    X_fit, Z_fit, Dz_dx, DDz_dx = Fit_Model(m, varphi, Data, z_max)

    # line1 = ax.scatter(Data['X'], Data['Z'], color=colors[varphi], s=1, label=f"Black Water (φ={varphi})")
    # line2 = ax.scatter(x, y, color="red", s=1, label="Adjusted curve")
    line1, = ax.plot(Data['X'], Data['Z'], color= colors[varphi], linewidth = 1, markersize = 1, label=f"Black Water (φ={varphi})")
    line2, = plt.plot(X_fit, Z_fit, label="Ajuste con curve_fit", color="red", linewidth = .4)

    ax.legend(loc='upper right', title="Leyenda", fontsize=6, title_fontsize=8)
    
    ax_T = plt.axes([0.18, 0.025, 0.5, 0.06])
    slider_T = Slider(ax_T, 'Time', -10*(2*m), -4.5*(2*m), valinit=T_init, valstep=1)
    slider_T.valtext.set_visible(False)

    ax_text = plt.axes([0.71, 0.025, 0.075, 0.06])  # Espacio para el cuadro de texto
    text_box = TextBox(ax_text, 'T', initial=str(T_init))
    text_box.label.set_position((-0.2, 0.5)) 

    ax_button = plt.axes([0.775, 0.025, 0.1, 0.06])  # Espacio para el botón
    button = Button(ax_button, 'Update')

    loading_text = ax.text(0, - 6.75/2 * m, f'Derivative = {Dz_dx}, Second Derivative = {DDz_dx}', ha='center', fontsize=12, color='darkred')

    def mostrar_carga(mensaje):
        loading_text.set_text(mensaje)
        fig.canvas.draw_idle()

    def actualizar1(val):
        mostrar_carga("Cargando...")  
        S = slider_T.val
        Data = Black_water_drop(m, varphi, S)
        Data = Data[(Data['Z'] < z_max)]
        # line1.set_offsets(np.c_[Data['X'], Data['Z']])
        line1.set_xdata(Data['X'])
        line1.set_ydata(Data['Z'])

        print(calcular_area_dataframe(Data))
        X_fit, Z_fit, Dz_dx, DDz_dx = Fit_Model(m, varphi, Data, z_max)
        line2.set_xdata(X_fit)
        line2.set_ydata(Z_fit)
        mostrar_carga(f'Derivative = {Dz_dx}, Second Derivative = {DDz_dx}')  # Mostrar "Listo" después de terminar
        text_box.set_val(str(S)) 
        fig.canvas.draw_idle()

    def actualizar2(event):
        try:
            S = float(text_box.text)  
            print(f"Actualizando el valor del tiempo a {S} desde el cuadro de texto")
            slider_T.set_val(S)  
        except ValueError:
            print("Por favor, ingrese un número válido en el cuadro de texto.")

    slider_T.on_changed(actualizar1)
    button.on_clicked(actualizar2)

    sh = m * 3
    ax.set_xlim(-sh , sh)
    ax.set_ylim(-sh, z_max)
    plt.show()

# Crear figura y eje
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)

m = 50
varphi = 1
Plot_XZ_var(m, varphi, -1*m)
