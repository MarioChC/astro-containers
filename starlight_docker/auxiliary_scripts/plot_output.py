#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 15:23:12 2024

@author: mario
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
from matplotlib import rcParams
import pandas as pd
from io import StringIO

wd = "/Users/mario/kk/astro-containers/starlight_docker/starlight/shared_directory/output/"
output_file = "output_spectrum_xpos_20_ypos_26_fiber_number_0887_NGC7025_LR-V.txt"

file_path = wd + output_file

# Palabras clave para buscar en las secciones
keywords = ["Synthesis Results - Best model", "l_obs f_obs f_syn wei"]

# Abre el fichero en modo de lectura
with open(file_path, 'r') as file:
    # Lee el contenido del fichero
    content = file.read()
    
    
# Divide el contenido en secciones usando '##' como separador
sections = content.split('##')

# Elimina espacios en blanco al principio y al final de cada sección
sections = [section.strip() for section in sections]

# Filtra las secciones que contienen alguna de las palabras clave y guarda sus índices
selected_sections_indices = [(i, section) for i, section in enumerate(sections) if any(keyword.lower() in section.lower() for keyword in keywords)]

# Muestra los índices y las secciones seleccionadas
# for i, (index, section) in enumerate(selected_sections_indices):
#     print(f"Sección {index + 1} (índice {index}):\n{section}\n")
    

#######################
### Bestfit results ###
#######################

resultados_bestfit = sections[selected_sections_indices[0][0] + 1]

chi2_Nleff = float(resultados_bestfit.split('\n')[0].split()[0])
adev = float(resultados_bestfit.split('\n')[1].split()[0])

velocity = float(resultados_bestfit.split('\n')[8].split()[0]) # km/s
velocity_dispersion = float(resultados_bestfit.split('\n')[9].split()[0]) # km/s
A_V = float(resultados_bestfit.split('\n')[10].split()[0]) # mag

print(" chi2_Nleff =", chi2_Nleff, '\n', "adev =", adev, '\n', "velocity =", velocity, 'km/s\n', 
      "velocity_dispersion =", velocity_dispersion, 'km/s\n', "A_V =", A_V, "mag")


# Leer los datos de M_ini, M_corr, edad, metalicidad, (L/M)_j...


### Eliminamos el primer carácter '#' de la primera línea para obtener los nombres de las columnas
lines = '\n'.join(resultados_bestfit.split('\n')[13:]).strip().split('\n')
header = lines[0][1:].strip().split()  # Eliminamos '#' y dividimos por espacios

### Creamos un nuevo texto sin las líneas que comienzan con '#'
data_text = '\n'.join([line for line in lines[1:]])

### Creamos un DataFrame de pandas
df_results = pd.read_csv(StringIO(data_text), delim_whitespace=True, header=None, names=header)


light_weighted_age = df_results['x_j(%)']*df_results['age_j(yr)']*1e-2 # Divido por 100 porque los pesos están en %
average_light_weighted_age = light_weighted_age.sum()*1e-9 # Gyr
light_weighted_z = df_results['x_j(%)']*df_results['Z_j']*1e-2
average_light_weighted_z = light_weighted_z.sum()


#############################################
### Plotear espectro junto con el bestfit ###
#############################################

datos_ajuste = selected_sections_indices[1][1].split("[Nl_obs]\n")[1]


# Dividir la cadena en líneas
lines = datos_ajuste.strip().split('\n')

# Inicializar listas vacías para cada columna
columna1 = []
columna2 = []
columna3 = []
columna4 = []

# Recorrer cada línea y dividir por espacios
for line in lines:
    parts = line.split()
    if len(parts) == 4:  # Asegurarse de que haya 4 partes por línea
        columna1.append(float(parts[0]))  # Convertir a float el primer elemento
        columna2.append(float(parts[1]))  # Convertir a float el segundo elemento
        columna3.append(float(parts[2]))  # Convertir a float el tercer elemento
        columna4.append(float(parts[3]))  # Convertir a float el cuarto elemento
        
# Convertir listas a arrays de NumPy
array1 = np.array(columna1)
array2 = np.array(columna2)
array3 = np.array(columna3)
array4 = np.array(columna4)

# Índices de puntos enmascarados
indices_enmascarados = np.where(array4==0)[0]

# Crear una máscara booleana para los índices
mask_below = np.zeros_like(array1, dtype=bool)
mask_below[indices_enmascarados] = True




rcParams['font.family'] = 'serif'
rcParams['axes.linewidth'] = 1

plt.rc('xtick', direction = 'in')
fig = plt.figure(figsize=(8,3.7))
plt.clf()
frame1=fig.add_axes((.1,.3,.8,.63))
plt.title("STARLIGHT")  
plt.plot(array1,array2, 'k', label='Input spectrum', linewidth=0.5)
plt.plot(array1,array3, 'r', label='STARLIGHT fit', linewidth=0.5)
# plt.axvspan(5400,5500,color="yellow")
ax = plt.gca()
# Colorear la región bajo la curva según la máscara
ax.fill_between(array1, np.max(array2), y2=np.min(array2), where=mask_below, color='lightskyblue', edgecolor='none', alpha=0.7, interpolate=True)
ax.xaxis.set_tick_params(length = 5, width=1,labelsize=0)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_tick_params(length = 5, width=1,labelsize=12)
ax.yaxis.set_minor_locator(AutoMinorLocator())          
plt.ylabel("Flux", fontsize = 11)
plt.setp(ax.spines.values(), linewidth=1,zorder=100)
plt.legend(frameon = False)
frame2=fig.add_axes((.1,.15,.8,.14))   
plt.plot(array1,array2-array3, 'g', label='Residuals', linewidth=0.5)
plt.xlabel("Observed wavelength [$\mathregular{\AA}$]", fontsize = 12)  
plt.rc('axes',linewidth=1.5)
ax = plt.gca() 
# Colorear la región bajo la curva según la máscara
ax.fill_between(array1, np.max(array2-array3), y2=np.min(array2-array3), where=mask_below, color='lightskyblue', edgecolor='none', alpha=0.7, interpolate=True)
ax.xaxis.set_tick_params(length = 5, width=1,labelsize=12)
ax.yaxis.set_tick_params(length = 5, width=1,labelsize=12)
ax.xaxis.set_ticks_position('both')
plt.rc('xtick', direction = 'in')                     
plt.setp(ax.spines.values(), linewidth=1,zorder=100)
plt.subplots_adjust(left =0.1, bottom =0.2, right =0.9, top =0.99)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.savefig(wd+"starlight_spectral_fitting.pdf", dpi=600)
#############################################
################## FIN DE ###################
### Plotear espectro junto con el bestfit ###
#############################################
