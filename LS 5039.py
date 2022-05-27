# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 02:03:43 2021

@author: Sergi
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import scipy.optimize as opt
import scipy.stats as stats

'''
Aquest codi va ser d'ús personal i pot no estar adequadament comentat. Correspon a un codi emprat per al meu TFG on, donades una sèrie de dades,
s'efectua l'anàlisi de diversos models astrofísics per determinar quin és el més plausible. Aquest programa es va utilitzar diverses vegades editant
manualment el codi per extreure anàlisis de diversos conjunts de dades així com era convenient durant el procés d'elaboració del TFG. Això significa
que executar aquest programa no donarà el resultat per a totes les sèries de dades sinó només per a la sèrie seleccionada en el codi.
Degut a la naturalesa de les dades, hi ha paràmetres de l'ajust molt difícils de determinar i poden provocar que l'execució del codi s'allargui més
de l'esperat. Tot i així, el codi es deixava d'executar automàticament passar un cert temps (curt, màxim uns pocs minuts) degut al propi funcionament
de scipy.optimize.
'''

# Dades

Dies = [50856.2, 50883.1, 50913.0, 50946.0, 52540.96, 52545.07, 52546.96,
        52549.04, 52552.03, 52554.98, 52559.01, 52559.93, 52561.03, 52562.02,
        52562.98, 52564.96, 52566.01, 52566.91, 52567.95, 52568.94, 53220.57,
        53377.21, 53391.20, 53392.20, 53644.60, 53720.17, 53809.04, 54519.00,
        54573.81, 54596.78, 54688.63, 56491.56, 56492.60, 56492.75, 56493.56,
        56494.57, 56494.75, 56495.73] # MJD

freqs = np.array([0.235, 0.61, 1.4, 2.3, 4.8, 8.5, 15]) # GHz

S_235 = np.array([27, 18, 28, 26, 31, 19, 25, 34, 24]) # mJy
error_235 = np.array([9, 6, 9, 8, 9, 6, 8, 11, 7])

S_610 = np.array([51, 28, 38, 49, 46, 41, 49, 35, 31.3, 40, 51, 37])
error_610 = np.array([5, 4, 3, 3, 3, 2, 3, 2, 1.8, 2, 3, 2])

S_14 = np.array([27, 38.6, 39.8, 45.5, 36.8, 29.2, 27.5, 33.6, 32, 34.5, 38.5,
                 42, 31.6, 35.0, 36.4, 26.2, 29.1])
error_14 = np.array([4, 1.0, 0.7, 1.0, 1.4, 1.7, 1.3, 1.2, 2, 1.9, 1.1, 2, 1.7,
                     1.8, 1.3, 1.5, 1.6])

S_23 = np.array([30.6, 27.2])
error_23 = np.array([1.0, 0.7])

S_48 = np.array([21.8, 23.7, 22.7, 26.0, 20.5, 24.4, 22.3, 15.6, 18.8, 15.6,
                 20.0, 25.6, 20.25, 25.1, 29.6, 21.6, 23.6, 23.9, 18.9, 24.1,
                 22.42, 18.68])
error_48 = np.array([0.5, 0.1, 0.1, 0.2, 0.5, 0.2, 0.4, 0.5, 0.3, 0.6, 0.4,
                     0.2, 0.19, 0.2, 0.2, 0.2, 0.5, 0.2, 0.2, 0.3, 0.15, 0.11])

S_85 = np.array([17.4, 17.77, 15.85, 18.23, 14.3, 16.5, 14.7, 9.7, 11.3, 6.3,
                 16.7, 18.0, 14.5, 16.5, 20.8, 14.9, 13.0, 16.6, 14.8, 18.6])
error_85 = np.array([0.1, 0.08, 0.07, 0.12, 0.5, 0.3, 0.3, 0.5, 0.4, 0.5, 0.3,
                     0.2, 0.3, 1.0, 0.3, 0.3, 0.6, 0.2, 0.2, 0.3])

S_15 = np.array([12.4, 13.3, 11.8, 13.9, 9.4, 8.4, 6.5, 8.4, 13.3, 12.0, 10.8,
                 10, 15.3, 9.3, 7.9, 3.7, 11.2, 12.3])
error_15 = np.array([0.3, 0.2, 0.18, 0.2, 0.7, 0.7, 0.6, 0.4, 0.4, 0.5, 0.5, 2,
                     0.6, 0.7, 0.7, 0.3, 0.6, 0.4])

diccionari_S = {'0': S_235, '1': S_610, '2': S_14, '3': S_23, '4': S_48,
                '5': S_85, '6': S_15}
diccionari_errors = {'0': error_235, '1': error_610, '2': error_14,
                     '3': error_23, '4': error_48, '5': error_85,
                     '6': error_15}

# Programa

ticksx = [0.2, 0.5, 1, 2, 5, 10, 15]
ticksy = [5, 10, 20, 50]

weights = {}
divisor = {}
mesures = []
desvest = []
errors = []
errors2 = []
mitjana = []
for k, v in diccionari_errors.items():
    weights[k] = 1/(v**2)
    divisor[k] = np.sum(1/(v**2))
    errors = np.append(errors, (1/divisor[k]) *
                       np.sqrt(np.sum((weights[k]*diccionari_errors[k])**2)))

for k, v in diccionari_S.items():
    mesures = np.append(mesures, np.sum(weights[k]*v)/divisor[k])
    mitjana = np.append(mitjana, np.average(diccionari_S[k]))
    desvest = np.append(desvest, np.std(diccionari_S[k]) /
                        np.sqrt(len(diccionari_S[k])))

# Primer tempteig índex espectral

Giga = freqs[4:7]
Reg = mesures[4:7]
Sense = stats.linregress(np.log(Giga), np.log(Reg))
print(Sense[0])
print(Sense[1])
print(Sense[4])

# Funcions per ajustar

def SSA(nu, p, P1, P2):
    return (P1/(4*np.pi)) * nu**(5/2) * (1-np.exp(-P2*nu**(-(p+4)/2)))

def FFA(nu, p, P1, P2):
    return (P1/(4*np.pi)) * nu**((5-p)/2) * (1-np.exp(-P2*nu**(-2)))

def SSAFFA(nu, p, P1, R, K1, K2):
    return (1/(4*np.pi*(K1*nu**(-(p+4)/2) + K2*nu**(-2)))) *\
            P1*nu**(-(p-1)/2) * (1-np.exp(-R(K1*nu**(-(p+4)/2)+K2*nu*(-2))))

def SSARazin(nu, nuR, p, P1, P2):
    return np.exp(-nuR/nu) * SSA(nu, p, P1, P2)

def FFARazin(nu, nuR, p, P1, P2):
    return np.exp(-nuR/nu) * FFA(nu, p, P1, P2)

def SSAFFARazin(nu, nuR, p, P1, R, K1, K2):
    return np.exp(-nuR/nu) * SSAFFA(nu, p, P1, R, K1, K2)

def Razin(nu, nuR, p, P1):
    return np.exp(-nuR/nu) * (P1/(4*np.pi)) * nu**(-(p-1)/2)


# Ajust

mask = [True, True, True, False, True, True, True]

popt, pcov = opt.curve_fit(SSA, freqs[mask], mesures[mask],
                           sigma=desvest[mask], p0=[2, 500, 3])
print(popt)

popt1, pcov1 = opt.curve_fit(FFA, freqs[mask], mesures[mask],
                             sigma=desvest[mask], p0=[2, 500, 3])
print(popt1)

popt2, pcov2 = opt.curve_fit(SSARazin, freqs[mask], mesures[mask],
                             sigma=desvest[mask], p0=[0.4, 2, 500, 3],
                             bounds=((0, -np.inf, -np.inf, -np.inf),
                                     (np.inf, np.inf, np.inf, np.inf)),
                             absolute_sigma=True)
print(popt2)

popt3, pcov3 = opt.curve_fit(FFARazin, freqs[mask], mesures[mask],
                             sigma=desvest[mask], p0=[0.1, 2, 1000000, 0.01],
                             absolute_sigma=True)
print(popt3)

razin, crazin = opt.curve_fit(Razin, freqs[mask], mesures[mask],
                              sigma=desvest[mask], p0=[0.4, 2, 800])
print(razin)
print(crazin)

# Chi quadrat

Chisqu = stats.chisquare(mesures[mask], SSARazin(freqs[mask], popt2[0], popt2[1],
                                                 popt2[2], popt2[3]))

Chisqu2 = stats.chisquare(mesures[mask], FFA(freqs[mask], popt1[0], popt1[1],
                                             popt1[2]))

Chisqu3 = stats.chisquare(mesures[mask], SSA(freqs[mask], popt[0], popt[1],
                                             popt[2]))

Chisqu4 = stats.chisquare(mesures[mask], FFARazin(freqs[mask], popt3[0], popt3[1],
                                                  popt3[2], popt3[3]))

Chisqu5 = stats.chisquare(mesures[mask], Razin(freqs[mask], razin[0], razin[1],
                                               razin[2]))
print(Chisqu)
print(Chisqu2)
print(Chisqu3)
print(Chisqu4)
print(Chisqu5)

# Imatges

x = np.logspace(-1, 1.3, 100)

fig, ax = plt.subplots()
plt.yscale('log')
plt.xscale('log')
plt.tick_params(axis='both', which='both', direction='in', right=True,
                top=True)
ax.set_ylim(2, 70)
ax.set_xlim(0.1, 20)
# plt.plot(freqs, (freqs**(Sense[0])*np.exp(Sense[1])), '--')
ax.errorbar(freqs, mesures, yerr=2*desvest, fmt='s', markersize=5,
            capsize=3, label='Average')
# ax.errorbar(freqs, mitjana, yerr=2*desvest, fmt='s', markersize=5,
            # capsize=3, label='Simple average')
plt.plot(x, SSA(x, popt[0], popt[1], popt[2]), '--', label='SSA', color='k')
plt.plot(x, FFA(x, popt1[0], popt1[1], popt1[2]), ':', label='FFA', color='k')
# plt.plot(x, SSARazin(x, popt2[0], popt2[1], popt2[2], popt2[3]), '-',
         # label='SSA+Razin', color='k')
# plt.plot(x, FFARazin(x, popt3[0], popt3[1], popt3[2], popt3[3]), '--')
plt.plot(x, Razin(x, razin[0], razin[1], razin[2]), label='Razin', color='k')
locs, labels = plt.xticks()
plt.xticks(ticksx, ticksx)
plt.xlabel('$\\nu$ (GHz)')
plt.ylabel('$S_\\nu$ (mJy)')
locs, labels = plt.yticks()
plt.yticks(ticksy, ticksy)
plt.legend()
plt.tight_layout()
