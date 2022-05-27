# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 09:32:27 2020

@author: Sergi
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import scipy.special as special

'''
This program uses scipy to solve differential equations for pendulums, the Schrödinger equation, coupled harmonic oscillators and planetary orbits. Then plots
the solutions using matplotlib.
Detailed commentary is in catalan.
'''

# Pèndols

# Paràmetres

g = 9.8  # m/s^2
L = 1  # m
m = 1 # kg
t_inicial = 0  # s
t_final = 10  # s
N_punts = 10000

# Concidicions inicials

ci_pendol = [0, 0.01*np.pi] 
'''
Aquestes condicions inicials són d'angles petits. És d'esperar que les figures
del pèndol simple i del pèndol físic siguin pràcticament iguals. Unes
condicions properes a [0, 2*np.pi] mostraran les diferències entre ambdós
pèndols.
'''

# Temps

t = np.linspace(t_inicial, t_final, N_punts)

# Definició de les equacions diferencials dels pèndols

def eqdif1(theta, t):
    '''
    Funció que estableix les equacions del moviment del pèndol simple.
    theta és un vector tal que
    theta[0] = angle
    theta[1] = velocitat
    t és un vector amb els temps on calcularem el moviment
    '''
    return theta[1], -(g/L)*theta[0]

def eqdif2(theta, t):
    '''
    Funció que estableix les equacions del moviment del pèndol físic.
    theta és un vector tal que
    theta[0] = angle
    theta[1] = velocitat
    t és un vector amb els temps on calcularem el moviment
    '''
    return theta[1], -(g/L)*np.sin(theta[0])

def eqdif3(theta, t):
    '''
    Funció que estableix les equacions del moviment del pèndol esmorteït.
    theta és un vector tal que
    theta[0] = angle
    theta[1] = velocitat
    t és un vector amb els temps on calcularem el moviment
    b és un paràmetre que modela la intensitat de l'esmorteïment
    '''
    b = 0.5
    return theta[1], -(g/L)*np.sin(theta[0]) - (1/m*L**2)*(b*theta[1])
                          
def eqdif4(theta, t):
    '''
    Funció que estableix les equacions del moviment del pèndol esmorteït
    forçat.
    theta és un vector tal que
    theta[0] = angle
    theta[1] = velocitat
    t és un vector amb els temps on calcularem el moviment
    b és un paràmetre que modela la intensitat de l'esmorteïment
    A és un paràmetre que modela la intensitat de la força
    omega és un paràmetre que modela la freqüencia de la força
    '''
    b = 0.5
    A = 1.35
    omega = 0.666
    g = 1
    return theta[1], -(g/L)*np.sin(theta[0]) + (1/m*L**2)*(-b*theta[1]
                                                          + A*np.cos(omega*t))

# Resolució dels pèndols
Pendol_simple = spi.odeint(eqdif1, ci_pendol, t)
Pendol_fisic = spi.odeint(eqdif2, ci_pendol, t)
Pendol_esmorteit = spi.odeint(eqdif3, ci_pendol, t)
Pendol_esm_for = spi.odeint(eqdif4, ci_pendol, t)

plt.figure()
# Plot pèndol simple
plt.title('Pèndol simple')
plt.plot(t, Pendol_simple[:, 0], label='Angle')
plt.plot(t, Pendol_simple[:, 1], label='Velocitat (rad/s)')
plt.xlabel('t')
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.figure()
# Plot pèndol físic
plt.title('Pèndol físic')
plt.plot(t, Pendol_fisic[:, 0], label='Angle')
plt.plot(t, Pendol_fisic[:, 1], label='Velocitat (rad/s)')
plt.xlabel('t')
plt.grid(True)
plt.legend()
plt.tight_layout()

# Plot pèndol esmorteït
plt.figure()
plt.title('Pèndol esmorteït')
plt.plot(t, Pendol_esmorteit[:, 0], label='Angle')
plt.plot(t, Pendol_esmorteit[:, 1], label='Velocitat (rad/s)')
plt.xlabel('t')
plt.grid(True)
plt.legend()
plt.tight_layout()

# Plot pèndol esmorteït forçat
plt.figure()
plt.title('Pèndol esmorteït')
plt.plot(t, Pendol_esm_for[:, 0], label='Angle')
plt.plot(t, Pendol_esm_for[:, 1], label='Velocitat (rad/s)')
plt.xlabel('t')
plt.legend()
plt.grid(True)
plt.tight_layout() 
'''
El resultat d'aquest darrer plot indica que el forçament és prou gran per 
provocar gairebé una rotació. El pèndol queda instantàniament a l'inrevés per 
després tornar en sentit oposat i evitar completar una rotació.
'''
# %% Equació de Schrödinger

x_mig = 0
x_final = 5
x_inicial = -5

x1 = np.linspace(x_mig, x_inicial, 10000)
x2 = np.linspace(x_mig, x_final, 10000)
x_sencer = np.linspace(x_inicial, x_final, 20000)

# Representació de les funcions d'Hermite

def Hermite(x, n):
    '''
    Funció que avalua la funció d'Hermite d'ordre n en el punt x.
    '''
    return (((2**n) * special.factorial(n) * np.sqrt(np.pi))**(-1/2)) * \
        np.exp(-(x**2)/2)*special.eval_hermite(n, x)

Hermite_0 = []
Hermite_2 = []
Hermite_4 = []
for x in x_sencer:
    Hermite_0 = np.append(Hermite_0, Hermite(x, 0))
    Hermite_2 = np.append(Hermite_2, Hermite(x, 2))
    Hermite_4 = np.append(Hermite_4, Hermite(x, 4))


plt.figure()
plt.subplot(121)
plt.title('Funcions Hermite')
plt.plot(x_sencer, Hermite_0, label='n=0')
plt.plot(x_sencer, Hermite_2, label='n=2')
plt.plot(x_sencer, Hermite_4, label='n=4')
plt.grid(True)
plt.xlabel('x')
plt.legend()
plt.xlim(-5, 5)
plt.ylim(-0.6, 0.8)
plt.tight_layout()

# Resolució equació de Schrödinger

def eqsch(psi, x):
    '''
    Funció que resol l'equació diferencial de Schrödinger
    determinada. 
    psi[0] = funció d'ona
    psi[1] = derivada
    '''
    return psi[1], (x**2 - (2*n + 1))*psi[0] 

enes = [0, 2, 4]
ci_set = [[0.751126, 0], [-0.531125, 0], [0.459969, 0]]
colors = ['C0', 'C1', 'C2']

plt.subplot(122)
plt.title('Solució equació Schrödinger')
plt.grid(True)
plt.xlabel('x')
plt.xlim(-5, 5)
plt.ylim(-0.6, 0.8)
for n in enes:
    sch_neg = spi.odeint(eqsch, ci_set[int(n/2)], x1)
    sch_pos = spi.odeint(eqsch, ci_set[int(n/2)], x2)

    plt.plot(x1, sch_neg[:, 0], colors[int(n/2)], label = 'n = ' + str(n))
    plt.plot(x2, sch_pos[:, 0], colors[int(n/2)])
plt.legend()
plt.tight_layout()

# %% Oscil·ladors acoblats

# Paràmetres

m_osc = 1 # kg
k1 = 10 # N/m
k2 = 0.5 # N/m
temps = np.linspace(0, 40, 10000)
ci_osc = [1, 0, 0, 0]

# Resolució

def harmosc(x, t):
    '''
    Funció que resol l'equació de dos oscil·ladors acoblats
    x[0] = x1
    x[1] = v1
    x[2] = x2
    x[3] = v2
    '''
    return x[3], -(k1 + k2) * x[2]/m + k2 * x[0]/m, \
        x[1], -(k1 + k2) * x[0]/m + k2 * x[2]/m

oscilador = spi.odeint(harmosc, ci_osc, temps)

plt.figure()
# Plot de les posicions
plt.subplot(121)
plt.title('Posicions dels oscil·ladors')
plt.xlabel('Temps (s)')
plt.ylabel('Posició (m)')
plt.plot(temps, oscilador[:, 0], label = 'x1')
plt.plot(temps, oscilador[:, 2], label = 'x2')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Plot de les velocitats
plt.subplot(122)
plt.title('Velocitats dels oscil·ladors')
plt.xlabel('Temps (s)')
plt.ylabel('Velocitat (m/s)')
plt.plot(temps, oscilador[:, 1], label = 'v1')
plt.plot(temps, oscilador[:, 3], label = 'v2')
plt.legend()
plt.grid(True)
plt.tight_layout()

# %%

# Paràmetres Terra

T_Terra = 365.26*86400 # s
N_punts = 10000
t_Terra = np.linspace(0, T_Terra, N_punts)
au = 1.49598e11 # m
afeli_Terra = 1.0167*au # m
M = 1.9891e30 # kg
G = 6.67e-11 # unitats SI
v_Terra = 29290 # m/s
ci_Terra = [afeli_Terra, 0, 0, v_Terra]

# Paràmetres Halley

T_Halley = 75.32*3.15576e7 # s
afeli_Halley = 35.082*au # m
v_Halley = 0.869e3 # m/s
t_Halley = np.linspace(0, T_Halley, N_punts)

ci_Halley = [afeli_Halley, 0, 0, v_Halley]

# Resolució

def gravitacio(y, t):
    '''
    Funció que resol l'equació del moviment d'un objecte de massa negligible 
    en un camp gravitatori.
    y[0] = x
    y[1] = v_x
    y[2] = y
    y[3] = v_y
    '''
    r_quadrat = y[0]**2+y[2]**2
    return y[1], -G*M*y[0]/(r_quadrat)**(3/2), \
        y[3], -G*M*y[2]/(r_quadrat)**(3/2)

Terra = spi.odeint(gravitacio, ci_Terra, t_Terra)
Halley = spi.odeint(gravitacio, ci_Halley, t_Halley)

# Figura Terra
plt.figure()

# Subplot posició x
plt.subplot(221)
plt.grid(True)
plt.title('x(t)')
plt.xlabel('Temps (s)')
plt.ylabel('Posició (m/s)')
plt.plot(t_Terra, Terra[:, 0])
plt.tight_layout()

# Subplot posició y
plt.subplot(222)
plt.grid(True)
plt.title('y(t)')
plt.xlabel('Temps (s)')
plt.ylabel('Posició (m/s)')
plt.plot(t_Terra, Terra[:, 2])
plt.tight_layout()

# Subplot velocitat
plt.subplot(223)
plt.grid(True)
plt.title('Velocitat tangencial')
plt.xlabel('Temps (s)')
plt.ylabel('Velocitat (m/s)')
plt.plot(t_Terra, np.sqrt(Terra[:, 1]**2 + Terra[:, 3]**2))
plt.ylim(20000, 40000)
plt.tight_layout()

# Subplot òrbita
plt.subplot(224)
plt.grid(True)
plt.title('Òrbita')
plt.xlabel('Posició x (m)')
plt.ylabel('Posició y (m)')
plt.plot(Terra[:, 0], Terra[:, 2])
plt.axis('equal')
plt.tight_layout()

# Figura Halley
plt.figure()

# Subplot posició x
plt.subplot(221)
plt.grid(True)
plt.title('x(t)')
plt.xlabel('Temps (s)')
plt.ylabel('Posició (m/s)')
plt.plot(t_Halley, Halley[:, 0])
plt.tight_layout()

# Subplot posició y
plt.subplot(222)
plt.grid(True)
plt.title('y(t)')
plt.xlabel('Temps (s)')
plt.ylabel('Posició (m/s)')
plt.plot(t_Halley, Halley[:, 2])
plt.tight_layout()

# Subplot velocitat
plt.subplot(223)
plt.grid(True)
plt.title('Velocitat tangencial')
plt.xlabel('Temps (s)')
plt.ylabel('Velocitat (m/s)')
plt.plot(t_Halley, np.sqrt(Halley[:, 1]**2 + Halley[:, 3]**2))
plt.tight_layout()

# Subplot òrbita
plt.subplot(224)
plt.grid(True)
plt.title('Òrbita')
plt.xlabel('Posició x (m)')
plt.ylabel('Posició y (m)')
plt.plot(Halley[:, 0], Halley[:, 2])
plt.axis('equal')
plt.tight_layout()
