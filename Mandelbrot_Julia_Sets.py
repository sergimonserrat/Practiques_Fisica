# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 02:03:14 2021

@author: Sergi
"""
import numpy as np
import matplotlib.pyplot as plt
import imageio

'''
This program computes the Mandelbrot set and generates a gif showing the convergence in the complex plane. Then a Julia set is also computed and represented
graphically.
'''

# Mandelbrot

# Dimensions
llargada = 500
amplada = 500

# Inicialització dels arrays
z_n = np.zeros([amplada, llargada])
CM = np.zeros([amplada, llargada])
VC = np.zeros([amplada, llargada])
N = 24

# Paràmetre c
c_inferior = - 2 - 1.5j
c_superior = 1 + 1.5j


# Inicialització del pla
Re, Im = np.meshgrid(np.linspace(np.real(c_inferior), np.real(c_superior),
                                 amplada), np.linspace(np.imag(c_inferior),
                                                       np.imag(c_superior),
                                                       llargada))
pla_constants = Re + Im*1j

# Generació de les figures
fig_1 = plt.figure(1)
fig_2 = plt.figure(2)
imatges = []

for n in range(1, N+1):
    mask = np.abs(z_n) > 100
    z_n = z_n*(1 - mask) + mask*100
    z_n = z_n**2 + pla_constants
    CM = np.abs(z_n) < 2
    VC = VC + CM
    
    plot1 = fig_1.add_subplot(4, 6, n)
    plot1.imshow(CM, cmap = 'gray')
    plot1.axis('off')
    
    plot2 = fig_2.add_subplot(4, 6, n)
    plot2.imshow(VC, cmap = 'jet')
    plot2.axis('off')
    
    if n < 10:
        nom_fitxer = 'CM_0' + str(n) + '.jpg'
    else:
        nom_fitxer = 'CM' + str(n) + '.jpg'

    plt.imsave(nom_fitxer, VC)
    imatges.append(imageio.imread(nom_fitxer))
    
imageio.mimsave('gif_CM.gif', imatges, duration = 1)

#%% Julia

# Inicialització dels arrays
z_julia = np.zeros([amplada, llargada])
CJ = np.zeros([amplada, llargada])
VC_Julia = np.zeros([amplada, llargada])

# Paràmetre c
cinf_julia = - 1.5 - 1.5j
csup_julia = 1.5 + 1.5j


# Inicialització del pla
Re, Im = np.meshgrid(np.linspace(np.real(cinf_julia), np.real(csup_julia),
                                 amplada), np.linspace(np.imag(cinf_julia),
                                                       np.imag(csup_julia),
                                                       llargada))
z_julia = Re + Im*1j # z_0 = c

# Generació de les figures
fig_3 = plt.figure(3)
fig_4 = plt.figure(4)
imatjulia = []

for n in range(1, N+1):
    maskjulia = np.abs(z_julia) > 10000
    z_julia = z_julia*(1 - maskjulia) + maskjulia*10000
    z_julia = z_julia**2 - 0.7269 + 0.1889*1j
    CJ = np.abs(z_julia) < 2
    VC_Julia = VC_Julia + CJ
    
    plot3 = fig_3.add_subplot(4, 6, n)
    plot3.imshow(CJ, cmap = 'gray')
    plot3.axis('off')
    
    plot4 = fig_4.add_subplot(4, 6, n)
    plot4.imshow(VC_Julia, cmap = 'jet')
    plot4.axis('off')
    
    if n < 10:
        nom_fitxer = 'CJ_0' + str(n) + '.jpg'
    else:
        nom_fitxer = 'CJ' + str(n) + '.jpg'

    plt.imsave(nom_fitxer, VC_Julia)
    imatjulia.append(imageio.imread(nom_fitxer))
    
imageio.mimsave('gif_CJ.gif', imatjulia, duration = 1)
