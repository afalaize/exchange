# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 15:24:06 2017

@author: afalaize
"""

opening = '{'
closing = '}'

string = '{blabla0= \n{blabla2=blabla3 }{blabla4={blabla6=blablablabla7}}}'


def test_bloc(string):
    """
    Return True if the string contains as many closing delimiters as opening
    delimiters.
    """
    return string.count(opening) == string.count(closing)

    
def get_bloc_limits(string):
    """
    Return the indices for opening and closing delimiters of the next bloc.
    """
    start = string.find(opening)
    if start == -1:
            return start
    else:
        end = string.find(closing)
        while not test_bloc(string[start:end]):
            end = end + string[end:].find(closing) + len(closing)
        return start, end

        
def get_bloc(string):
    """
    Return the content of the next bloc in string, without delimiters
    """
    limits = get_bloc_limits(string)
    if limits == -1:
        return limits
    else: 
        return string[limits[0]+1:limits[1]-1].strip()

        
def get_bloc_key_and_value(bloc_string):
    """
    For a given bloc (returned by the get_bloc function), recover the mtc key
    with associated value.
    """
    try:
        key, value = bloc_string.split('=', 1)
    except ValueError:
        key, value = bloc_string.split(':', 1)
    key = key.strip()
    value = value.strip()
    return key, value

    
class Bloc:
    """
    Basically a recusrsive dictionnary
    """
    def __init__(self, string):
        bloc_string = get_bloc(string)
        key, value = get_bloc_key_and_value(bloc_string)
        self.key = key
        value_limits = get_bloc_limits(value)
        if value_limits == -1:
            self.value = value
        else:
            self.value = Bloc(value)

    def string(self):
        if isinstance(self.value, str):
            value = self.value
        else:
            value = self.value.string()
        return opening + ' ' + self.key + '= ' + value + ' ' + closing

        
def get_blocs_recursively(string):
    bloc_string = get_bloc(string)
    key, value = get_bloc_key_and_value(bloc_string)
    limits = get_bloc_limits(value)
    if limits == -1:
        return {key: value}
    else:
        continuation_flag = True
        subblocs = list()
        while continuation_flag:
            limits = get_bloc_limits(value)
            if limits == -1:
                continuation_flag = False
            else:
                subblocs.append(get_blocs_recursively(value[limits[0]:limits[1]]))
                value = value[limits[1]:]
        return {key: subblocs}
      
if __name__ == '__main__':
    dico = get_blocs_recursively(string)
        
        