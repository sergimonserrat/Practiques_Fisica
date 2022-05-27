# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 13:01:11 2020

@author: Sergi
"""
import numpy as np
import matplotlib.pyplot as plt
import random

'''
This program attempts to use montecarlo techniques to calculate pi. Unsure it works.
'''

pi = 0
k = 0
error = np.array([])
tolerancia = 1e-8
while np.abs(np.pi - pi) > tolerancia:
    terme = ((1/16)**k) * ((4/(8*k+1)) - (2/(8*k+4)) - (1/(8*k+5)) - (1/8*k+6))
    pi = pi + terme
    error = np.append(error, np.abs(np.pi - pi))
    print(pi)
    k = k + 1


#%% serveix per crear equivalents a les celes del jupyter
import random
import numpy as np
L = 1
Impactes_quadrat = 0
Impactes_cercle = 0
fita = 1e-7
monte_error = 1
errorray = np.array([])
while monte_error > fita:
    x = random.uniform(0, L)
    y = random.uniform(0, L)
    Impactes_quadrat = Impactes_quadrat + 1
    r = (x**2 + y**2) ** (1/2)
    if r <= L:
        Impactes_cercle = Impactes_cercle + 1
    monte_error = abs(np.pi - 4*(Impactes_cercle/Impactes_quadrat))
    errorray = np.append(errorray, monte_error)

montepi = 4 * (Impactes_cercle/Impactes_quadrat)
print(montepi)

#%% integració simpson millor per a dades/taula de dades experimentals
"""
quadratura millor per funcions donades. quad crida el nom de la funció, no
el nombre
"""

