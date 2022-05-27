# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 22:01:57 2022

@author: Sergi
"""
import numpy as np
import matplotlib.pyplot as plt
import string
import re
from unidecode import unidecode

# Empty list for every word as they appear
x = []
# Empty list for all characters
y = []
with open('1342-0.txt', 'r', encoding='utf8') as inputFile:
    for line in inputFile.readlines():
        line = unidecode(line) #Converting fancy quotation marks
        lineStripped = line.strip() #Removing leading white spaces
        lineSplitted = lineStripped.split() #Splitting into words
        for element in lineSplitted:
            #Removing punctuation marks
            element = element.translate(str.maketrans('', '', string.punctuation))
            #Removing special character double quotes 
            re.sub('"', '', element)
            #Converting to lowercase
            element = element.lower()
            #Adding to the list
            x.append(str(element))

#Empty dictionary
recompte = dict()

for paraula in x:
    if paraula in recompte:
        #Word is already in dictionary
        recompte[paraula] = recompte[paraula] + 1
    else:
        #New key for word not yet in dictionary
        recompte[paraula] = 1

#Empty array for probabilities
probs = np.array([])

for value in recompte:
    valor = recompte.get(value) #Obtaining values
    probs = np.append(probs, valor/len(x)) #Calculating probabilities
    
entropia = np.sum(-probs*np.log2(probs)) #Calculating word entropy

print('Pride and Prejudice word entropy is %.3f bits' %(entropia))

'''
English word entropy found in the literature is 9.8 bits per word
(Grignetti 1964). The lower entropy of Pride and Prejudice points to a higher
degree of order compared with a text of random English words. Given a word it
is easier to predict the next one than a randomly generated English text.
'''

lletres=dict()
ps = np.array([])
for element in x:
    for lletra in element:
        #Running through the text one letter at a time
        if lletra in lletres:
            #Character is already in dictionary
            lletres[lletra] = lletres[lletra] + 1
        else:
            #New key for character not yet in dictionary
            lletres[lletra] = 1

for value in lletres:
    valor = lletres.get(value) #Obtaining values
    ps = np.append(ps, valor) #Calculating probabilities

P = ps/np.sum(ps)
charentropy = np.sum(-P*np.log2(P)) #Calculating character entropy

print('Pride and Prejudice character entropy is %.3f bits' %(charentropy))

'''
English character entropy is 4.03 bits per letter using the first order model
described in Elements of Information Theory by T. M. Cover & J. A. Thomas, 2006
p. 170.
The higher value for Pride and Prejudice might be due to the inclusion of
numbers as characters, as they slightly contribute to the final information
content.
'''