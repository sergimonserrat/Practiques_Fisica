# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6 17:54:41 2022

@author: Sergi
"""
import numpy as np
import matplotlib.pyplot as plt

mida = np.linspace(0, 2, 101)
w = 1.9
eta = 3
x, y = np.meshgrid(mida, mida)

laplace = np.zeros((101, 101))

distancia = np.sqrt((x-1)**2+(y-1)**2)

exterior = np.abs(distancia) > 0.99
interior = np.abs(distancia) < 1.01

contorn = interior ^ exterior
# Condicions de contorn
for i in range(0, 101):
    for j in range(0, 101):
        if contorn[i][j] == False:
            laplace[i][j] = 1

contorn[50][50] = False

nointegrar = np.abs(distancia) > 1.01

def solver(contorn, laplace, w, nointegrar):
    # Resolució sobrerrelaxació
    for k in range(1, 50): 
        for i in range(0, 100):
            for j in range(0, 100):
                if nointegrar[i][j] == False:
                    if contorn[i][j] == True:
                        laplace[i][j] = (1-w)*laplace[i][j] + w/4 * \
                                        (laplace[i][j+1] + laplace[i][j-1] \
                                         + laplace[i+1][j] + laplace[i-1][j])
    return laplace

lichtenberg = np.zeros((101, 101))
lichtenberg[50][50] = 1
candidats = np.zeros((101, 101))
check = lichtenberg.astype(np.bool)

def candidatures(candidats, lichtenberg, laplace):
    for i in range(1, 100):
        for j in range(1, 100):
            test = lichtenberg[i+1][j] + lichtenberg[i-1][j] \
                   + lichtenberg[i][j-1] + lichtenberg[i][j+1]
            if lichtenberg[i][j] == 0:
                if test > 0:
                    if laplace[i][j] > 0:
                        candidats[i][j] = test*laplace[i][j]
    return candidats

def sorteig(probabilitats):
    x=np.arange(10201).reshape(101,101)
    xy=np.random.choice(x.flatten(), 1, p = probabilitats.flatten())
    index=np.unravel_index(xy.item(), x.shape) # convert back to index
    return index

check = 0
with open('distancia_eta3.txt', 'w') as f:
    while check <= 0.99: 
        solver(contorn, laplace, w, nointegrar)
        candidats = np.zeros((101, 101))
        candidatures(candidats, lichtenberg, laplace)
        probabilitats = (candidats**eta)/np.sum(candidats**eta)
        i, j = sorteig(probabilitats)
        lichtenberg[i][j] = 1
        laplace[i][j] = 0
        contorn[i][j] = False
        check = distancia[i][j]
        f.write(str(check))
        f.write('\n')
        print(check)

# %%
plt.figure()
plt.title('$\eta$=3')
plt.imshow(exterior, cmap='Reds', extent= [0,2,0,2])
plt.imshow(lichtenberg, cmap = 'binary', extent = [0,2,0,2], alpha = 0.9)
plt.tight_layout()

