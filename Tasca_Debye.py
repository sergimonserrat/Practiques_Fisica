# -*- coding: utf-8 -*-
"""
@author: Sergi Monserrat Mascaró

Realitzat amb Python 3.8.3, entorn Spyder 4.1.4
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import scipy.stats as stats

# Paràmetres

a1 = 87.9
b1 = 0.404
c1 = 9.59e-4
d1 = 1.33e-6
a2 = 80.7
b2 = 4.42e-3
c2 = 1.37e-13
d2 = 651
T0 = 133

# Funció permitivitat

def epsilon(nu):
    '''
    Aquesta funció calcula la permitivitat en funció de la freqüència nu
    i la temperatura Temp, ambdues són entrades tipus float.
    Retorna un complex amb el valor de la funció calculat.
    '''
    epsilon_0 = a1 - b1*Temp + c1*Temp**2 - d1*Temp**3
    epsilon_inf = epsilon_0 - a2*np.exp(-b2*Temp)
    tau = c2*np.exp(d2/(Temp + T0))
    epsilon = epsilon_inf + (epsilon_0 - epsilon_inf)/(1 - 1j*2*np.pi*nu*tau)
    return epsilon

# Funció resta. Útil per determinar talls amb optimize solve

def resta(nu):
    '''
    Aquesta funció realitza la diferència entre la part real i la imaginària
    del complex obtingut de la funció epsilon. Quan aquesta funció és 0, tenim
    un punt de tall.
    L'entrada nu és tipus float i correspon a la freqüència.
    Retorna un float amb el valor de la diferència.
    '''
    resta = np.real(epsilon(nu)) - np.imag(epsilon(nu))
    return resta

# Arrays de temperatures i freqüències que volem com a input de la funció

Temperatures = np.linspace(0, 100, 6)
Frequencies = np.logspace(8, 12, 10000) 
# Tot i utilitzar 10000 punts per les freqüències, a nivell gràfic la figura
# generada serà una corba suau fins i tot amb només 100 punts, però a nivell
# numèric s'empitjora la precisió dels resultats.

# Array per emmagatzemar punts de tall

'''
L'array Talls_max contindrà els punts de tall propers al màxim. L'array
Talls_cua conté els punts de tall situats a la cua de la part imaginària
'''
Talls_cua = np.array([])
Talls_max = np.array([])

# Graficació i manipulació de la funció epsilon

plt.figure()
plt.grid(True)
plt.xlabel('Freqüència (Hz)')
plt.ylabel('$\epsilon_r$ i $\epsilon_i$')
plt.tight_layout()
for Temp in Temperatures:
    # Separació part real i imaginària
    epsilon_i = np.imag(epsilon(Frequencies))
    epsilon_r = np.real(epsilon(Frequencies))
    # Graficació
    plt.semilogx(Frequencies, epsilon_r, label = 'T =' + str(int(Temp)) + ' Re')
    plt.semilogx(Frequencies, epsilon_i, label = 'T =' + str(int(Temp)) + ' Im')
    # Determinació del màxim imaginari
    epsilon_imax = epsilon_i.argmax()
    nu_imax = Frequencies[epsilon_imax]
    # Determinació de talls propers al màxim
    Talls_max = np.append(Talls_max, optimize.fsolve(resta, nu_imax))
    print('T =', Temp, 'Nu màx =',  format(nu_imax, 'e'))
    print('Max =', optimize.fsolve(resta, nu_imax))
    # Determinació de talls a la cua
    Talls_cua = np.append(Talls_cua, optimize.fsolve(resta, Frequencies[-1]))
    print('Cua =', optimize.fsolve(resta, Frequencies[-1]))
    # S'aprofita l'índex [-1] que permet indexar l'array començant pel final
plt.legend() #Si no es posa la llegenda aquí, no detecta les labels

# Regressió

Regressio_max = stats.linregress(Talls_max, Temperatures)
Regressio_cua = stats.linregress(Talls_cua, Temperatures)

print('Regressió del màxim:')
print('Pendent =', Regressio_max[0])
print('Terme independent =', Regressio_max[1])
print('Correlació =', Regressio_max[2])

print('Regressió de la cua:')
print('Pendent =', Regressio_cua[0])
print('Terme independent =', Regressio_cua[1])
print('Correlació =', Regressio_cua[2])

# Graficació dels punts de tall i la regressió
plt.figure()
plt.subplot(1, 2, 1)
# Subplot dels punts de tall prop del màxim
plt.grid(True)
plt.xlabel('Freqüència de tall (Hz)')
plt.ylabel('Temperatura (ºC)')
plt.title('Punts de tall prop del màxim')
plt.plot(Talls_max, Temperatures, 'r+')
plt.plot(Talls_max, Regressio_max[0]*Talls_max+Regressio_max[1], 'b')
plt.tight_layout()
plt.subplot(1, 2, 2)
# Subplot dels punts de tall a la cua
plt.grid(True)
plt.xlabel('Freqüència de tall (Hz)')
plt.ylabel('Temperatura (ºC)')
plt.title('Punts de tall a la cua')
plt.plot(Talls_cua, Temperatures, 'm+')
plt.plot(Talls_cua, Regressio_cua[0]*Talls_cua+Regressio_cua[1], 'g')
plt.tight_layout()

