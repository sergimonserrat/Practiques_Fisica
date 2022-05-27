# -*- coding: utf-8 -*-
"""
Created on Thu May 12 13:32:22 2022

@author: Sergi
"""

def StringInsertion(string, inserted, index):
    '''
    Function that inserts a string at a particular point in another string
    
    Input
    string: str. Original string
    inserted: str. String we want to insert
    index: int. Position where we wish to insert the string
    
    Output
    newstring: str. The new string.
    '''
    newstring = string[:index] + inserted + string[index:]
    return newstring

def findWord(chars):
    '''
    Function that constructs words from the comparison of its constituent
    characters.
    
    Input
    chars: list. List of strings that compares the relative position of the
    individual letter within the word. Only works for '>' symbol.
    
    Output
    Word: str. The word we wanted to find.
    '''
    Word = str(chars[0][0]+chars[0][2]) #We stitch the first two words
    del chars[0] #We delete the info we already used
    while len(chars) != 0: #Loop that keeps going untill we run out of comparisons
        i = 0 #Useful counter to delete used info later
        for comparison in chars:
            if comparison[0] in Word: #Check if first character is already in our word
                #We pull the index so we know where to insert it
                index = Word.index(comparison[0])
                #We update our word
                UpdatedWord = StringInsertion(Word, comparison[2], index+1)
                #We reset the starting string
                Word = UpdatedWord
                #We delete the comparison we just used
                del chars[i]
            elif comparison[2] in Word:
                #Same loop to check the second character
                index = Word.index(comparison[2])
                if index != 0: #Conditional loop to avoid 0-index problems
                    UpdatedWord = StringInsertion(Word, comparison[0], index-1)
                else:
                    UpdatedWord = StringInsertion(Word, comparison[0], 0)
                Word = UpdatedWord
                del chars[i]
            else:
                pass
            i = i+1
    return Word
            
chars1 = ['P>E','E>R','R>U']
chars2 = ["I>N","A>I","P>A","S>P"]
chars3 = ["I>F", "W>I", "S>W", "F>T"]
chars4 = ["R>T", "A>L", "P>O", "O>R", "G>A", "T>U", "U>G"]
chars5 = ["A>S", "C>A", "S>A"]

print(findWord(chars1))
print(findWord(chars2))
print(findWord(chars3))
print(findWord(chars4))
print(findWord(chars5))

'''
Coding time: aprox. 1h
Difficulty: medium. I had to think through carefully how to methodically check
everything.
'''