# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 09:37:34 2020

@author: Sergi
"""
import numpy as np
import matplotlib.pyplot as plt

""" 
python lliga variables que contenen el mateix array i modificar-ne un implica
modificar l'altre. Cal emprar np.copy per deslligar
"""


def test1(x):
    return a*x + b


def test2(x, a, b):
    return a*x + b


def test3(x, a, b):
    return a*x + b


a = 3
b = 2

print(test1(3))
print(test2(3, a, b))
print(test2(3, 3, 2))
print(test3(3, a, b))

"""
en aquest exemple es veu el funcionament de les variables. Si empréssim
classes, seria convenient mirar el funcionament de la comanda `global`
"""
#%% funcionament numpy.generic

array = np.array([3, 5, 7, 2, 0, 9])

"""
se li pot demanar `array.max()`, `array.min()`, `array.argmax()`, 
`array.argmin()`. On arg fa referència a la posició del max/min dins l'array.
Si posam `array[3:].argmax()` cercarà a partir de l'argument 3. Si hi ha més
d'un màxim, retorna el primer que trobi.
Més comandes `array.sum()`, `array.mean()`, `array.std()`, `array.size()`,
`array.shape()`, `array.dtype()`
Cercar a numpy.generic per més info.
""" 
#%% exercici permitivitat


"""
def perm(freq):
    copiar parametres del fitxer
    epsilon_0 = bla
    epsilon_inf = bla
    tau = bla
    epsilon = bla
    return epsilon
"""

"""

"""