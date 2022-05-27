# -*- coding: utf-8 -*-
"""
Created on Thu May 12 14:26:40 2022

@author: Sergi
"""

def getChange(amount, price):
    change = amount-price
    coins = [0.01, 0.05, 0.10, 0.20, 0.50, 1.00] #Available coin values
    coin_hist = [0 for i in range(0,6)] #Histogram of needed coins
    i=-1 #Start at the largest coin denomination
    while change < coins[i]: #Checks the largest denomination we can start from
        i = i-1
    while change != 0: 
        coin_hist[i] += 1 #Adds one coin at a time
        change = change - coins[i]
        change = round(change, 2) #Prevents rounding errors in the updated change
        try:
            while change < coins[i]:
                #Step down to smaller denomination when we can't keep substracting
                i = i-1
        except IndexError:
        #If we let it run it would give an error because the index is out of the list
            break
    return coin_hist

print(getChange(5, 0.99))
print(getChange(3.14, 1.99))
print(getChange(4, 3.14)) 
'''
There is a mistake on the formulation of the problem. The proposed solution
[1, 0, 1, 1, 1, 0] = 0.81 dollars of change while 4-3.14 = 0.86 dollars. The
solution [1, 1, 1, 1, 1, 0] presented here is correct.
'''
'''
Coding time: aprox. 30 min
Difficulty: easy
'''