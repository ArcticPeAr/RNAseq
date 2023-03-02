from matplotlib_venn import venn2

vennDict = {}

def vennForGo(listUP, listDOwn, GOterm):
    '''This function takes two lists (Up and Down) and their respective GO term and returns a dictionary with the GO term as the key and the lists as the values.'''
    vennList = [listUP, listDOwn]
    vennDict[GOterm] = vennList
    return vennDict


def lastTwoGo(goDict):
    '''This function takes the dictionary from the vennForGo function and plots a venn diagram for the last two GO terms.'''
    venn2(subsets = (goDict[goDict.keys()[-2]], goDict[goDict.keys()[-1]]), set_labels = (goDict.keys()[-2], goDict.keys()[-1]))
    plt.title('Venn Diagram for ' + goDict.keys()[-2] + ' and ' + goDict.keys()[-1])
    plt.show()


def vennPlotter(vennDict):
    '''This function takes the dictionary from the vennForGo function and plots the venn diagram for each GO term and save it as Venn.'''
    for key in vennDict:
        venn2(subsets = vennDict[key], set_labels = (key, ''))
        plt.title(key)
        plt.show()

def UpDownVennPlotter(vennDict):
    '''This function accepts a dictionary and for each key, plots a venn diagram for the Up and Down lists and saves the diagram as "venn{Goterm}UpDown.pdf.'''
    for key in vennDict:
        venn2(subsets = (vennDict[key][0], vennDict[key][1]), set_labels = (key, ''))
        plt.title(key)
        plt.savefig('venn' + key + 'UpDown.pdf')


def vennIntersect(vennDict):
    '''This function accepts a dictinonary and produces a venn diagram of the intersection of Up and Down lists between the last two GO terms and saves it as "venn{Goterm1_Goterm2].pdf .'''
    venn2(subsets = (vennDict[vennDict.keys()[-2]][0], vennDict[vennDict.keys()[-1]][1]), set_labels = (vennDict.keys()[-2], vennDict.keys()[-1]))
    plt.title('Venn Diagram for ' + vennDict.keys()[-2] + ' and ' + vennDict.keys()[-1])
    plt.savefig('venn' + vennDict.keys()[-2] + '_' + vennDict.keys()[-1] + '.pdf')
    
