#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import math as m
import random as rnd
import matplotlib.tri as tri
from scipy.spatial import Voronoi, voronoi_plot_2d

import tile as t  #CREDITS: Inspired by https://github.com/alsignoriello/periodic_voronoi/tile.py 



##############################  Medidas de la red  ##############################

Ny=10
Nx=8
Ncells=Nx*2*Ny

L=1.0
h=(L/2.0)/m.tan(30.0/180.0*np.pi)

#################################################################################



##########################  START | Genero los puntos  ##########################

# Red triangular regular
x_reg=np.zeros(Ncells)
y_reg=np.zeros(Ncells)

for fila in range(2*Ny):
    if (fila%2==0):
        for col in range (Nx):
            x_reg[fila*Nx +col] = 3*col*L +L
            y_reg[fila*Nx +col] = fila*h +h
    else:
        for col in range (Nx):
            x_reg[fila*Nx +col] = (3*col+2)*L +L/2.0
            y_reg[fila*Nx +col] = fila*h +h

points = np.column_stack((x_reg, y_reg))        
"""
plt.figure("Puntos regulares")
plt.scatter(points[:,0], points[:, 1])
"""
# Opción A: varío la posición

gamma = 1.5
T=2*h # distancia al centro del vecino más cercano de la red regular, ¿varía?

mod = np.zeros((Ncells, 2))

for i in range (Ncells):
    rx=rnd.random()
    ry=rnd.random()
    mod[i, 0] = gamma/T*(rx-0.5)
    mod[i, 1] = gamma/T*(ry-0.5)


# fin opción A


# Opción B: importo la variación en la posición de los puntos
#mod = np.loadtxt('mod_TEJIDO_PRUEBA.txt')
# fin opción B


# Condiciones periódicas: Multiplico mis puntos para tener 9 cajas iguales:
all_points, x_len, y_len = t.tile_points(points, Ncells, Nx, Ny, mod)

main_points = points + mod
"""
plt.figure("Puntos modificados")
plt.scatter(main_points[:,0], main_points[:, 1])
"""
###########################  END | Genero los puntos  ###########################



######################  START | Creo el diagrama Voronoi  #######################

# Guardo la información dada por la función Voronoi
vor = Voronoi(all_points)
all_vertices = vor.vertices
all_regions = vor.regions


# Guardo la información de la caja central y las relaciones con la numeración general
index_map_vertices = {}
index_map_regions = {}
main_vertices = set()
main_regions = []

for p in range(Ncells):
    r=vor.point_region[p] # all = index_points_regions[main]
    region_dentro=all_regions[r] # indices en all de esta región
    main_regions.append(region_dentro) # indices en all
    index_map_regions[r]=p # index_map_regions[all] = main
    
    for i in range(len(region_dentro)):
        v = region_dentro[i]
        if v not in index_map_vertices:
            index_map_vertices[v] = len(main_vertices)
            vertice = all_vertices[v]
            main_vertices.add((vertice[0], vertice[1]))

main_vertices = np.array(list(main_vertices))

#######################  END | Creo el diagrama Voronoi  ########################



##############################  START | FUNCIONES  ##############################

def busco_analoga_dentro(centroid): 
    """ Dado un punto en plano, primero encuentra el centro de la región asociada 
    más cercana. Como esta se encuentra fuera de la caja central, luego, encuentra 
    la análoga de esa región dentro de la caja central. Para ello, prueba a sumar
    o restar la longitud y altura de la caja central. Devuelve el índice de la 
    región análoga a la inicial en main_regions. """
    
    distancias = np.linalg.norm(all_points - centroid, axis=1)
    ind = np.argmin(distancias)
    fuera = all_points[ind]
    
    x_fuera = fuera[0]
    y_fuera = fuera[1]
    # Valores de n_x y n_y para probar
    valores_n = [-1, 0, 1]
    
    # Iterar sobre los valores de n_x y n_y
    for n_x in valores_n:
        for n_y in valores_n:
            # Calcular las coordenadas (x0, y0)
            x = x_fuera + n_x * x_len
            y = y_fuera + n_y * y_len
            
            # Verificar si las coordenadas (x_dentro, y_dentro) corresponden a algún punto en main_points
            for i in range(len(main_points)):
                punto=main_points[i]
                x_dentro=punto[0]
                y_dentro=punto[1]
                if (abs(x-x_dentro) < 0.0000001) and (abs(y-y_dentro) < 0.0000001):
                    return i # Devolver el índice del punto en m_points
    
    
    # Si no se encuentra ningún punto, devolver None
    return None



def busco_clave(valor, diccionario): 
    for clave, val in diccionario.items():
        if val == valor:
            return clave
    return None



def dist_entre_vertices(vertices): 
    """ Dado un array con los índices de 2 vértices, devuelve la distancia
    entre ellos en el plano. """
    
    if len(vertices) != 2:
        print("La lista contiene ",len(vertices), " índices.")
        return None
    
    v1_index, v2_index = vertices
    
    if v1_index < 0 or v1_index >= len(all_vertices):
        print("El índice ", v1_index, " está fuera de los límites.")
        return None
    
    if v2_index < 0 or v2_index >= len(all_vertices):
        print("El índice ", v2_index, " está fuera de los límites.")
        return None
    
    # Obtener las coordenadas de los dos vértices
    v1 = all_vertices[vertices[0]]
    v2 = all_vertices[vertices[1]]
    #print(v1, v2, "\n")
    
    dist = np.linalg.norm(np.array(v1) - np.array(v2))
    #print(dist)
    
    return dist



def encontrar_indice_arista(p1, p2): 
    """ Dados los centros de dos regiones, encuentra el índice en ridge_points
    de la arista que separa estas dos regiones. """
    
    for index_a, pair in enumerate(vor.ridge_points):
        if list(pair) == [p1, p2] or list(pair) == [p2, p1]:
            return index_a
    return None


    
def distancia(p1,p2):
    """ Dados los centros de dos regiones, calcula la distancia de la arista
    que las separa. """
    
    indice_arista = encontrar_indice_arista(p1, p2)
     
    if (indice_arista == None):
        return None
    
    vertices = vor.ridge_vertices[indice_arista]
    dist = dist_entre_vertices(vertices)
    
    return dist
    


def invierto_point_region(): 
    region_point = {}
    for p in range(len(vor.point_region)):
        reg = vor.point_region[p]
        region_point[reg]= p
    
    return region_point
    

    
def contribucion_vecinos():
    # Crear una lista para almacenar las regiones vecinas de cada región en la caja central
    regiones_vecinas = []
    distancias = []
    
    region_point = invierto_point_region()
    
    # Iterar sobre las regiones de la caja central
    for i in range(len(main_regions)):
        region_actual = main_regions[i]
        
        # Crear una lista para almacenar las regiones vecinas de la región actual
        vecinas_region_actual = set()
        distancias_region_actual = []
        
        count=0
        
        #Compruebo que p1 = i: #####################
        i_all = busco_clave(i,index_map_regions)
        prueba = region_point[i_all]
        if prueba == i:
            p1=i
        else:
            print("Hay un error en la línea 278")
        ############################################
        
        # Iterar sobre los vértices de la región actual i
        for vertice in region_actual: # Cojo un vértice de mi región
            # Iterar sobre todas las regiones para buscar las vecinas
            for j in range(len(all_regions)):
                if vertice in all_regions[j]: # Entro si la región j tiene el vértice
                    
                    # Si la región j está dentro de la caja central:
                    if j in index_map_regions: # Si está dentro (está en el diccionario)
                        index = index_map_regions[j]
                        if (i != index):
                            count +=1
                            if index not in vecinas_region_actual:
                                vecinas_region_actual.add(index)
                                p2 = region_point[j]
                                
                                # Guardo la distancia de la arista
                                dist = distancia(p1,p2)
                                distancias_region_actual.append([i, index, dist])
                    
                    # Si NO está dentro, busco la análoga que sí lo está:
                    else:
                        region_vertices = all_regions[j]
                        centroid = np.mean(all_vertices[region_vertices], axis=0)
                        index = busco_analoga_dentro(centroid)
                        if (i != index):
                            count +=1
                            if index not in vecinas_region_actual:
                                vecinas_region_actual.add(index)
                                p2 = region_point[j]
                                
                                # Guardo la distancia de la arista
                                dist = distancia(p1,p2)
                                distancias_region_actual.append([i, index, dist])
                    
                
        # Agregar la lista de vecinos de la región actual a la lista de regiones vecinas
        regiones_vecinas.append(list(vecinas_region_actual))
        distancias.append(distancias_region_actual)
    
    return regiones_vecinas, distancias

###############################  END | FUNCIONES  ###############################



###################  START | Guardo los vecinos y compruebo  ####################

regiones_vecinas, distancias = contribucion_vecinos()

regiones_coincidentes =0
for i in range(Ncells):
    distancias_i = distancias[i]
    n_distancias = len(distancias_i)
    n_vertices = len(main_regions[i])
    if n_vertices != n_distancias:
        print("No se han encontrado todas las distancias para la región ", i)

    dist_coincidentes = 0
    for j in range(n_distancias):
        vecina = distancias_i[j]  
        ind_vecina = vecina[1]    
        dist_arista = vecina[2]   
        distancias_vecina = distancias[ind_vecina]
        for k in range(len(distancias_vecina)):
            if distancias_vecina[k][1]==i:
                if (abs(distancias_vecina[k][2] - dist_arista)>0.000000000001):
                    print("La arista entre las regiones", i, "y", ind_vecina, "no coincide.")
                else:
                    dist_coincidentes+=1
    if dist_coincidentes == n_distancias:
        regiones_coincidentes+=1
    else:
        print("Revisar errores. \n")

if regiones_coincidentes == Ncells:
    print("Todas las distancias coinciden.")

####################  END | Guardo los vecinos y compruebo  #####################



####################### START | Escribo todo en ficheros ########################

# write points #########################
"""
f = open("20240706.1.5.12x12.all_points.txt", "w+")
f.write("# index     x               y\n")
for i in range(len(all_points)):
    f.write("%d \t %f \t %f\n" % (i, all_points[i,0],all_points[i,1]))
f.close()

f = open("20240706.1.5.12x12.main_points.txt", "w+")
f.write("# index     x               y\n")
for i in range(len(main_points)):
    f.write("%d \t %f \t %f\n" % (i, main_points[i,0],main_points[i,1]))
f.close()


f = open("20240706.1.5.12x12.mod.txt", "w+")
f.write("# index     x               y\n")
for i in range(len(mod)):
    f.write("%d \t %f \t %f\n" % (i, mod[i,0],mod[i,1]))
f.close()

"""
"""
# write vertices #########################
f = open("20240525.regular.5-8.2-all_vertices.txt", "w+")
f.write("# index     x               y\n")
for i in range(len(all_vertices)):
    f.write("%d \t %f \t %f\n" % (i, all_vertices[i,0],all_vertices[i,1]))
f.close()

f = open("20240525.regular.5-8.3-main_vertices.txt", "w+")
f.write("# index     x               y\n")
for i in range(len(main_vertices)):
    f.write("%d \t %f \t %f\n" % (i, main_vertices[i,0], main_vertices[i,1]))
f.close()

f = open("20240525.regular.5-8.4-index_map_vertices.txt", "w+")
for index, (key, value) in enumerate(index_map_vertices.items()):
    f.write(f"[{index}] {key}: {value},\n")
f.close()


# write regions #########################
f = open("20240525.regular.5-8.5-all_regions.txt", "w+")
for i in range(len(all_regions)):
    poly = all_regions[i]
    if len(poly) != 0:
        f.write(f"[{i}]\t")
        for index in poly:
            f.write("%d\t" % index)
        f.write("\n")
f.close()

f = open("20240525.regular.5-8.6-main_regions.txt", "w+")
for i in range(len(all_regions)):
    poly = all_regions[i]
    if len(poly) != 0:
        f.write(f"[{i}]\t")
        for index in poly:
            f.write("%d\t" % index)
        f.write("\n")
f.close()

f = open("20240525.regular.5-8.7-index_map_regions.txt", "w+")
for index, (key, value) in enumerate(index_map_regions.items()):
    f.write(f"[{index}] {key}: {value},\n")
f.close()

  

# write neighbours #########################
f = open("20240525.regular.5-8.8-neighbours.txt", "w+")
for ind in range(Ncells):
	f.write(f"Neighbours of cell number {ind}:\t {regiones_vecinas[ind]}\n")
f.close()

"""
"""
# write distances #########################
f = open("20240706.1.5.12x12.distances.txt", "w+")
for ind in range(Ncells):
    for j in range(len(distancias[ind])):
        f.write(f"{distancias[ind][j]}\n")
f.close()
"""


######################### END | Escribo todo en ficheros ########################








########################## START | Histograma vecinos ###########################
n_vecinos = np.zeros(Ncells)
for i in range(Ncells):
    n_vecinos[i] = np.size(regiones_vecinas[i])

"""
plt.figure()
bins = np.arange(0.5, 10.5, 1)
hist, bin_edges = np.histogram(n_vecinos, bins=bins)
bcen = (bin_edges[1:] + bin_edges[:-1]) * 0.5
plt.bar(bcen, hist, width=0.9, edgecolor='black', align='center')
plt.xlabel('Número de vecinos', fontsize=17)
plt.ylabel('Frecuencia',  fontsize=17)
#plt.title('Histograma de número de vecinos')

plt.xticks(np.arange(3, 11, 1), fontsize=14)
plt.yticks(fontsize=14)  # Para que las etiquetas del eje x sean enteros
plt.xlim(3,10)
plt.tight_layout()
plt.show()
"""
#ruta="/Users/luciaabolea/Documents/UNIVERSIDAD/CUARTO/TFG/Latex/figuras"
#plt.savefig(f"{ruta}/gamma0_hist_vecinos.png", dpi=500)


########################### END | Histograma vecinos ############################






######################### START | Plot Diagrama Voronoi #########################

# Plottear el diagrama de Voronoi de la red triangular REGULAR  
"""
#plt.figure()
#plt.scatter(x,y)

# Crear la triangulación de Delaunay
triang = tri.Triangulation(x_reg, y_reg)

plt.figure("Red hexagonal regular")
plt.triplot(triang, 'o-')

"""

# Plottear el diagrama de Voronoi COMPLETO
voronoi_plot_2d(vor, show_vertices=False)

"""
# Labels
for region_id, point_id in index_map_regions.items():
    if region_id in vor.point_region:  # Verifica si la región está presente en el diagrama Voronoi
        region_index = vor.point_region[region_id]  # Obtiene el índice de la región en vor.regions
        region = vor.regions[region_index]
        #if -1 not in region:  # Ignora las regiones que tienen un vértice "infinite"
        # Calcula el centroide de la región
        centroid = vor.points[point_id] + 0.01
        # Agrega la etiqueta en el centroide de la región
        plt.text(centroid[0], centroid[1], str(index_map_regions[region_id]), fontsize=12, ha='center', va='center',  bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.3'))
"""



# Plotear el resto de los puntos de color 'lightpink'
#plt.scatter(all_points[Ncells:, 0], all_points[Ncells:, 1], c='lightpink', s=10)

#colors = ['purple'] * Ncells + ['lightpink'] * (len(all_points) - Ncells)
#voronoi_plot_2d(vor, show_vertices=False, point_colors=colors)



# Añadir puntos de la red triangular Delaunay (IRREGULAR)
triang = tri.Triangulation(all_points[:, 0], all_points[:, 1]) 
plt.triplot(triang, 'o-',color='lightpink',linewidth=0.4, alpha=0.7)

plt.scatter(all_points[:Ncells, 0], all_points[:Ncells, 1], c='xkcd:darkish pink', s=80)


"""triang_centro = tri.Triangulation(main_points[:, 0], main_points[:, 1]) 
plt.triplot(triang_centro, 'o-',color='purple',linewidth=0.55, alpha=0.8)
"""

""" #este quitar
# Plot recuadro naranja que delimita la main box:
main_box_x0 = np.zeros(20)
main_box_y0 = np.zeros(20)
main_box_x = np.linspace(0,x_len, num=20)
main_box_y = np.linspace(0,y_len, num=20)
main_box_x_len = np.ones(20)*x_len
main_box_y_len = np.ones(20)*y_len
plt.plot(main_box_x0, main_box_y,color='orange')
plt.plot(main_box_x, main_box_y0, color='orange')
plt.plot(main_box_x_len, main_box_y,color='orange')
plt.plot(main_box_x, main_box_y_len,color='orange')
"""


# Agregar líneas de cuadrícula cada 0.25 puntos en cada eje
#plt.grid(True, which='both', linestyle='-', linewidth=0.5)
#plt.gca().set_xticks(np.arange(0, 3 * Nx * L + 0.5, 0.5))
#plt.gca().set_yticks(np.arange(0, 2 * Ny * h + 0.5, 0.5))

#plt.gca().set_xticks(np.arange(-3, 3 * Nx * L + 2.5, 0.5))
#plt.gca().set_yticks(np.arange(-2.5, 2 * Ny * h + 2.5, 0.5))



# Plot de los vértices de la región 0:
"""
reg = main_regions[0]
verts_reg = set()

print(len(reg))

for i in range(len(reg)):
    v = reg[i]
    vertice = all_vertices[v]
    verts_reg.add((vertice[0], vertice[1]))
    
verts_reg = np.array(list(verts_reg))

x_vert = verts_reg[:,0]
y_vert = verts_reg[:,1]
plt.scatter(x_vert,y_vert, c="limegreen", s=40)
"""


# Plot de todos los vértices de main_vertices 
#x_main_vert = main_vertices[:,0]
#y_main_vert = main_vertices[:,1]
#plt.scatter(x_main_vert,y_main_vert, c="limegreen", s=40)


########################## END | Plot Diagrama Voronoi ##########################





"""
# Plottear el diagrama de Voronoi CAJA CENTRAL
vor_main = Voronoi(main_points)
voronoi_plot_2d(vor_main, show_vertices=False)

plt.scatter(main_points[:, 0], main_points[:, 1], c='purple', s=80)
"""
# No es el mismo que la caja central del completo (problema con las células del borde)


















