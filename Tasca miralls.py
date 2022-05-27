# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 01:49:08 2021

@author: Sergi
"""
import numpy as np
import matplotlib.pyplot as plt

'''
This program simulates a mirror composed by many layers of dielectric materials (TiO2, BK7, MgF2). Then the transmitance and reflectance are computed
and plotted as a function of the incident light wavelength. As a complementary exercise, the properties of a mirror with the layers arranged in
reverse order are also studied.
'''

# Funcions n

def nTiO2(landa):
    return np.sqrt(5.913+0.2441/(landa**2-0.0803))

def nBK7(landa):
    return np.sqrt(1 + 1.03961212/(1 - 0.00600069867/landa**2) +
                   0.231792344/(1 - 0.0200179144/landa**2) +
                   1.01046945/(1 - 103.560653/landa**2))

def nMgF2(landa):
    return np.sqrt(1 + 0.48755108/(1 - (0.04338408/landa)**2) +
                   0.39875031/(1 - (0.09461442/landa)**2) +
                   2.3120353/(1 - (23.793604/landa)**2))

# Paràmetres

N_punts = 600
landa_i = 400 # nanòmetres
landa_f = 1000 # nanòmetres
landa_ref = 632 #nanòmetres
N_capes = 10 # ha de ser parell

'''
Tot i que originalment l'exercici proposa 10 capes, afegir més capes tendeix a
una espècie de funció esglaó rectangular.
'''

# Inicialització 

x = np.linspace(landa_i, landa_f, N_punts)
R = np.array([])
R_invers = np.array([])

# Gràfica índex de refracció

plt.figure()
plt.grid(True)
plt.xlabel('$\lambda$ (nm)')
plt.ylabel('n')
plt.plot(x, nTiO2(x/1000), label='TiO$_2$') # conversió a micròmetres
plt.plot(x, nBK7(x/1000), label='BK7')
plt.plot(x, nMgF2(x/1000), label='MgF$_2$')
plt.legend()
plt.tight_layout()

# Matriu de transferència

def funcio_capa(n, long_ona, d):
    fase = 2 * np.pi * n * d/long_ona
    m_11 = np.cos(fase)
    m_12 = (np.sin(fase)*1j) / n
    m_21 = n*np.sin(fase)*1j
    m_22 = np.cos(fase)
    matriu = np.array([[m_11, m_12], [m_21, m_22]])
    return matriu

def Reflexio(B, C):
    return np.abs((B-C)/(B+C)) # n_0=1 omès
    
#Gruix capes

d_TiO2 = landa_ref / (4*nTiO2(landa_ref/1000))
d_BK7 = landa_ref / (4*nBK7(landa_ref/1000))
d_MgF2 = landa_ref / (4*nMgF2(landa_ref/1000))


for landa in x:
    n1 = nTiO2(landa/1000)
    n2 = nMgF2(landa/1000)
    n3 = nBK7(landa/1000)
    M1 = funcio_capa(n1, landa, d_TiO2)
    M2 = funcio_capa(n2, landa, d_MgF2)
    M_tot = np.linalg.matrix_power(np.dot(M1, M2), int(N_capes/2))
    ordre_invers = np.linalg.matrix_power(np.dot(M2, M1), int(N_capes/2))
    B_inv = ordre_invers[0][0]+n3*ordre_invers[0][1]
    C_inv = ordre_invers[1][0]+n3*ordre_invers[1][1]
    B = M_tot[0][0]+n3*M_tot[0][1]
    C = M_tot[1][0]+n3*M_tot[1][1]
    R = np.append(R, Reflexio(B, C))
    R_invers = np.append(R_invers, Reflexio(B_inv, C_inv))

# Graficació R, T
plt.figure()
plt.grid(True)
plt.title('Mirall 10 capes')
plt.xlabel('$\lambda$ (nm)')
plt.plot(x, R, label='Reflectància')
plt.plot(x, 1-R, label='Transmitància')
plt.legend()
plt.tight_layout()

'''
He decidit no posar etiqueta a l'eix y ja que és adimensional i s'utilitza per
representar alhora R i T. No em semblava que hi hagués un títol prou evident
i la llegenda és prou clara per interpretar la gràfica.
'''

# Experiment per veure l'efecte de girar l'ordre de les capes
plt.figure()
plt.grid(True)
plt.title('Ordre de les capes girat')
plt.xlabel('$\lambda$ (nm)')
plt.plot(x, R_invers, label='Reflectància')
plt.plot(x, 1-R_invers, label='Transmitància')
plt.legend()
plt.tight_layout()

'''
Girar l'ordre de les capes no és tan catastròfic. El comportament és semblant
al de l'ordre òptim. R s'apropa a 1 més lentament i és menys constant en la
finestra de treball. Si no fos perquè l'altre ordre és millor, podria ser un
mirall acceptable.
'''
