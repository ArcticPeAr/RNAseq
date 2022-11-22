"""
This program accepts DFs of up and down DEGs as lists for two or more GO terms and compares them
"""
#Import libraries
import pandas as pd

#Eac 
class GoList:
    def __init__(self, term, up, down):
        self.term = term
        self.up = up
        self.down = down
    
    def get_up(self):
        return self.up
    
    def get_down(self):
        return self.down   


def compareUpLists(list1, list2):
    """
    This function compares two lists of up DEGs and returns a list of genes that are in both lists and lists of genes that are in list1 but not in list2 and vice versa
    """
    return list(set(list1).intersection(list2))


def compareDownLists(list1, list2):
    """
    This function compares two lists of down DEGs and returns a list of genes that are in both lists and 
    """
    return list(set(list1).intersection(list2))

def outliersListUp(list1, list2):
    """
    This function returns a list of genes that are in list1 but not in list2
    """
    return list(set(list1).difference(list2))

def outliersListDown(list1, list2):
    """
    This function returns a list of genes that are in list1 but not in list2
    """
    return list(set(list1).difference(list2))

