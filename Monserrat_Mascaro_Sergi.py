# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 11:13:19 2021

@author: Sergi
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.integrate as integrate

'''
This code was presented as the answer for an exam during my bachelor's degree. It consisted of a .txt file containing spectral data from
the cosmic microwave background (CMB). This program attempts to calculate the CMB temperature from the spectrum using different methods.
First, a guessing method is used that sweeps an interval of possible answers and checks which one minimises the mean square error.
Second, a curve fitting using the black body radiation equation.
Finally, the Planck function is integrated and checked against the expected result (the Stefan-Boltzmann law)
'''
# Obrir el fitxer

xf = []
yf = []
with open('cmb-data.txt', 'r') as inputFile:
    for line in inputFile.readlines():
        lineStripped = line.strip()
        lineSplitted = lineStripped.split(';')
        xf.append(float(lineSplitted[0]))
        yf.append(float(lineSplitted[1]))
        
# Conversió en arrays

nu = np.array(xf)
Intensitat = np.array(yf)
long_ona = 1/nu
        
# Dades

c = 2.99792458e10 # velocitat de la llum
h = 6.62607015e-27 # constant de Planck
kB = 1.380649e-16 # constant de Boltzmann
sigma = 5.67e-5 # constant d’Stefan-Boltzmann

# Funcions 

def cos_negre(nu, T):
    return (2*h*(nu**3)*(c**2))/(np.exp(h*c*nu/(kB*T))-1)

def integrand(long_ona, T):
    return ((2*h*c**2)/(long_ona**5)) * (1/(np.exp(h*c/(kB*T*long_ona)) - 1))

# Tempteig MSE

ini = 1
T = np.linspace(2.7, 2.8, 100)
for temperature in T:
    MSE = (1/len(nu))*np.sum(Intensitat-cos_negre(nu, temperature))**2
    if MSE < ini:
        print('T(K):', temperature, 'MSE:', MSE)
        # el darrer valor printat és el millor
    ini = MSE
    
# Ajust

popt, pcov = opt.curve_fit(cos_negre, nu, Intensitat)
print('T òptima ajust (K):', popt[0])
    
# Graficació

x = np.linspace(2, 22, 100)    

plt.figure()
plt.grid(True)
plt.xlabel('$\\nu$')
# vaig descobrir que \\ em funciona millor per representar LaTeX al meu PC
plt.ylabel('I') # unitats? l'enunciat no les especifica
plt.title('Fons radiació microones')
plt.plot(nu, Intensitat, 'r+', label='Mesures COBE')
plt.plot(x, cos_negre(x, popt[0]), 'b-', label='Ajust cos negre') 
# per defecte ja sortia una línia blava
plt.tight_layout()

StefanBoltzman = integrate.quad(integrand, long_ona[33], long_ona[0],
                                args=popt[0])
print('Resultat integral:', StefanBoltzman[0])
print('Llei T^4:', (sigma*popt[0]**4)/np.pi)
