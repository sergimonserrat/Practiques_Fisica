# -*- coding: utf-8 -*-
"""
Created on Wed May 11 12:28:40 2022

@author: Sergi
"""
import numpy as np
import matplotlib.pyplot as plt

#%% Functions

def hextodec(string):
    '''
    Hexadecimal to decimal converter.
    
    Input
    string: str. Hexadecimal string
    
    Output
    int. Decimal integer
    '''
    return int(string, 16)
    
def DictBuilder(llista, empty_dict):
    '''
    Adds elements to a dictionary from a list of two elements. The first
    list element will be used as the key and the second as the value.
    
    Input
    llista: list. List of length 2 containing the pair of desired key and value
    empty_dict: dict. Dictionary we wish to update with a new key-value pair
    
    Output
    codes: dict. Returns the updated dictionary
    '''
    codes = empty_dict[str(llista[0])] = llista[1]
    return codes
    
def MaxRetriever(dictionary):
    '''
    Retrieves the key with the maximum value in a dictionary
    
    Input
    dictionary: dict. Dictionary from which we want to obtain the key with the
    maximum value
    
    Output
    str. The hexadecimal string that acts as the key of the dictionary associated
    to the maximum value
    '''
    return max(dictionary, key=dictionary.get)
#%% Initial values 

codes = {}
order = {}
Solution = ''
#%% Program
with open('DEV-CHALLENGE-1 data', 'r') as inputFile:
    for line in inputFile.readlines():
        lineStripped = line.strip()
        lineSplitted = lineStripped.split(',') #Designating comma as our separator
        #Storing hexadecimal strings with their associated letters in a dictionary
        DictBuilder(lineSplitted, codes)
        #Blank spaces might not be properly saved. We will deal with it later

with open('DEV-CHALLENGE-1 criteria', 'r') as inputFile:
    for line in inputFile.readlines():
        lineStripped = line.strip()
        lineSplitted = lineStripped.split(',') #Designating comma as our separator
        lineSplitted[1] = hextodec(lineSplitted[1]) #Conversion to decimal
        #Storing hexadecimal strings with their associated order numbers in a dictionary
        DictBuilder(lineSplitted, order)
  

for i in range(0, len(order)): #Iterating over the length of the dictionary
    key = MaxRetriever(order) #We retrieve the keys in reverse order
    newletter = codes[str(key)] #Obtaining the letter from the other dictionary
    if len(newletter) == 1:
        Solution = Solution + newletter #Adding the letter to the string
    else:
        Solution = Solution + ' '
        #This takes care of white spaces using the fact len(whitespace) != 1
    order[str(key)] = 0 #Setting value to 0 so we don't repeat letters indefinitely
    '''
    Since the size of a dictionary can't be changed during iteration, we need
    to update our values so they are not the maximum after each addition to the
    string. It's easier to work with the maximum than the minimum since we can
    guarantee that 0 is always a minimum, while multiplying by some factor would
    not guarantee as strongly that we had increased the used value enough for
    it to be a maximum.
    Choosing the maximum has the inconvenient of reversing the solution string
    '''
    Ordered_solution = Solution[::-1] #Reversing the final string
    
print('The enigma solution is:', Ordered_solution)

'''
Coding time: aprox. 2 h
Difficulty: relatively easy. The only problem was having to give up on using
dict.pop to get the values as this changed the size of the dictionary.
'''