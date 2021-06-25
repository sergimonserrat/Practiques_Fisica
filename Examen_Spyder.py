# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 10:39:14 2020

@author: Sergi
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special

# Graficació

def integrand(tau, x_prima):
    return (1/np.pi)*np.cos(x_prima*np.sin(tau))

x = np.linspace(0,20,1000)

Resolucio_special = np.array([])
Resolucio_def = np.array([])

for x_prima in x:
    Bessel_def = integrate.quad(integrand, 0, np.pi, args=x_prima)
    Resolucio_def = np.append(Resolucio_def, Bessel_def[0])
    Bessel_special = special.j0(x_prima)
    Resolucio_special = np.append(Resolucio_special, Bessel_special)

plt.figure()
# Subplot 1
plt.subplot(1, 2, 1)
plt.title('$J_0$(x) partir de scipy.special.j0')
plt.grid(True)
plt.plot(x, Resolucio_special)
plt.xlabel('x')
plt.ylabel('$J_0$(x)')
plt.tight_layout()

# Subplot 2
plt.subplot(1, 2, 2)
plt.title('$J_0$(x) partir de la definició')
plt.plot(x, Resolucio_def)
plt.grid(True)
plt.xlabel('x')
plt.ylabel('$J_0$(x)')
plt.tight_layout()

# Zeros de la funció
Zeros = []
for element in Resolucio_def:
    index = np.where(Resolucio_def == element)
    seguent = Resolucio_def[index[0]+1]
    test = element*seguent # Compara elements per trobar canvis de signe
    if test < 0:
        Zeros = np.append(Zeros, (x[index[0]]+x[index[0]+1])/2)
        # He fet la mitjana dels dos punts on hi ha el canvi de signe
        
    elif len(Zeros) == 5: # Per evitar problemes amb l'últim element de l'array
        break # Només funciona coneixent el número de zeros amb antelació
print('Zeros de la funció:', Zeros)

# Màxims i mínims
    