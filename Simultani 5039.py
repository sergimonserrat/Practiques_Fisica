# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 05:16:53 2021

@author: Sergi
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.stats as stats

# Dades

freqs = np.array([0.235, 0.61, 2.3, 4.8])

July19 = np.array([34, 51, 30.6, 22.42])
errors19 = np.array([11, 3, 1, 0.15])

July21 = np.array([24, 37, 27.2, 18.68])
errors21 = np.array([7, 2, 0.7, 0.11])

# Ticks gràfica

ticksx = [0.2, 0.5, 1, 2, 5, 10, 15]
ticksy = [5, 10, 20, 50]

# Funcions ajust

def SSA(nu, p, P1, P2):
    return (P1/(4*np.pi)) * nu**(5/2) * (1-np.exp(-P2*nu**(-(p+4)/2)))


def FFA(nu, p, P1, P2):
    return (P1/(4*np.pi)) * nu**((5-p)/2) * (1-np.exp(-P2*nu**(-2)))

def SSAFFA(nu, p, P1, P2, P3, K1, K2):
    return (1/(4*np.pi*(K1*nu**(-(p+4)/2) + K2*nu**(-2)))) *\
        P1*nu**(-(p-1)/2) * (1-np.exp(-P2*nu**(-(p+4)/2)-P3*nu*(-2)))

def SSARazin(nu, nuR, p, P1, P2):
    return np.exp(-nuR/nu) * SSA(nu, p, P1, P2)

def FFARazin(nu, nuR, p, P1, P2):
    return np.exp(-nuR/nu) * FFA(nu, p, P1, P2)

def SSAFFARazin(nu, nuR, p, P1, P2, P3, K1, K2):
    return np.exp(-nuR/nu) * SSAFFA(nu, p, P1, P2, K1, K2)

def Razin(nu, nuR, p, P1):
    return np.exp(-nuR/nu) * (P1/(4*np.pi)) * nu**(-(p-1)/2)

# 19 Juliol

# Ajustos

opt19, cov19 = opt.curve_fit(SSA, freqs, July19, sigma=errors19,
                             p0=[2, 500, 3])
print(opt19)
print(cov19)

opt119, cov119 = opt.curve_fit(FFA, freqs, July19, sigma=errors19,
                               p0=[2, 500, 3])
print(opt119)

opt219, cov219 = opt.curve_fit(SSARazin, freqs, July19,
                               sigma=errors19, p0=[0.3, 2, 500, 3],
                               bounds=((0, -np.inf, -np.inf, -np.inf),
                                       (np.inf, np.inf, np.inf, np.inf)),
                               absolute_sigma=True)
print(opt219)

opt319, cov319 = opt.curve_fit(FFARazin, freqs, July19,
                               sigma=errors19, p0=[0.3, 7, 700, 100],
                               bounds=((0, -np.inf, -np.inf, -np.inf),
                                       (np.inf, np.inf, np.inf, np.inf)),
                               absolute_sigma=True)
print(opt319)

# Chi quadrat

Chisqu19 = stats.chisquare(July19, SSARazin(freqs, opt219[0], opt219[1],
                                            opt219[2], opt219[3]))

Chisqu219 = stats.chisquare(July19, FFA(freqs, opt119[0], opt119[1],
                                        opt119[2]))

Chisqu319 = stats.chisquare(July19, SSA(freqs, opt19[0], opt19[1], opt19[2]))

Chisqu419 = stats.chisquare(July19, FFARazin(freqs, opt319[0], opt319[1],
                                             opt319[2], opt319[3]))
print('SSARazin 19 Juliol:', Chisqu19)
print('FFA 19 Juliol:', Chisqu219)
print('SSA 19 Juliol:', Chisqu319)
print('FFARazin 19 Juliol:', Chisqu419)

# 21 Juliol

# Ajustos

opt21, cov21 = opt.curve_fit(SSA, freqs, July21, sigma=errors21,
                             p0=[2, 3000, 0.2])
print(opt21)

opt121, cov121 = opt.curve_fit(FFA, freqs, July21, sigma=errors21,
                               p0=[5, 400, 30])
print(opt121)

opt221, cov221 = opt.curve_fit(SSARazin, freqs, July21,
                               sigma=errors21, p0=[0.41, 2.24, 5.59e8, 1.22e-6],
                               bounds=((0, -np.inf, -np.inf, -np.inf),
                                       (np.inf, np.inf, np.inf, np.inf)))
print(opt221)

opt321, cov321 = opt.curve_fit(FFARazin, freqs, July21,
                               sigma=errors21, p0=[0.4, 6, 600, 70],
                               bounds=((0, -np.inf, -np.inf, -np.inf),
                                       (np.inf, np.inf, np.inf, np.inf)))
print(opt321)

razin, crazin= opt.curve_fit(Razin, freqs, July21,
                               sigma=errors21, p0=[0.4, 2, 700])

print(razin)
print(crazin)

# Chi quadrat

Chisqu21 = stats.chisquare(July21, SSARazin(freqs, opt221[0], opt221[1],
                                            opt221[2], opt221[3]))

Chisqu221 = stats.chisquare(July21, FFA(freqs, opt121[0], opt121[1],
                                        opt121[2]))

Chisqu321 = stats.chisquare(July21, SSA(freqs, opt21[0], opt21[1], opt21[2]))

Chisqu421 = stats.chisquare(July21, FFARazin(freqs, opt321[0], opt321[1],
                                             opt321[2], opt321[3]))
Chisqu521 = stats.chisquare(July21, Razin(freqs, razin[0], razin[1], razin[2]))


print('SSARazin 21 Juliol:', Chisqu21)
print('FFA 21 Juliol:', Chisqu221)
print('SSA 21 Juliol:', Chisqu321)
print('FFARazin 21 Juliol:', Chisqu421)
print('Razin 21 Juliol:', Chisqu521)

x = np.logspace(-1, 1.3, 100)

fig, ax = plt.subplots()
plt.yscale('log')
plt.xscale('log')
plt.tick_params(axis='both', which='both', direction='in', right=True,
                top=True)
ax.set_ylim(2, 70)
ax.set_xlim(0.1, 20)
# plt.plot(freqs, (freqs**(Sense[0])*np.exp(Sense[1])), '--')
ax.errorbar(freqs, July19, yerr=errors19, fmt='s', markersize=5,
            capsize=3, color='C0', label='July 19')
ax.errorbar(freqs, July21, yerr=errors21, fmt='s', markersize=5,
            capsize=3, color='C1', label='July 21')
plt.plot(x, SSA(x, opt19[0], opt19[1], opt19[2]), '--', color='k',
         label='July 19 SSA')
# plt.plot(x, FFARazin(x, opt321[0], opt321[1], opt321[2], opt321[3]), '--',
         # color='C2')
plt.plot(x, Razin(x, razin[0], razin[1], razin[2]), ':', color='k',
         label='July 21 Razin')
# plt.plot(x, FFARazin(x, popt3[0], popt3[1], popt3[2], popt3[3]), '--')
locs, labels = plt.xticks()
plt.xticks(ticksx, ticksx)
plt.xlabel('$\\nu$ (GHz)')
plt.ylabel('$S_\\nu$ (mJy)')
locs, labels = plt.yticks()
plt.yticks(ticksy, ticksy)
plt.legend()
plt.tight_layout()
