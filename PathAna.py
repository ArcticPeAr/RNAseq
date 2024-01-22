from reactome2py import content, analysis
import pprint
import webbrowser
import itertools
import os
#from datetime import datetime
#now = datetime.now()
##date_time = now.strftime("%d-%m-%Y_%H:%M")

def mergeGenLists(geneList1,geneList2):
    '''This function takes gene lists of both up and down, removes genes appearing more than once and returns a merged list.'''
    mergedList = geneList1 + geneList2
    mergedList = list(set(mergedList))
    return mergedList
        

def getToken(geneList):
    '''This function takes a list of genes and returns a token.'''
    result = analysis.identifiers(ids=geneList)
    token = result['summary']['token']
    return token


def getURL(token):
    '''This function takes a token and returns a URL to reactome.'''
    url = 'https://reactome.org/PathwayBrowser/#/DTAB=AN&ANALYSIS=' + token
    return url


def openBrowser(url):
    '''This function takes a URL and opens it in the default browser.'''
    webbrowser.open(url)


def getReport(token):
    '''This function takes a token and returns a report.'''
    if not os.path.exists('report/'):
        os.makedirs('report/')
    analysis.report(token, path='report/', file='report.pdf', number='25', resource='TOTAL', 
                diagram_profile='Modern', analysis_profile='Standard', fireworks_profile='Barium Lithium', species= 'Homo sapiens', chunk_size=128)
    print('Report saved to report/ folder.')
    


def getFireworks(token):
    '''This function takes a token and returns a fireworks report.'''
    if not os.path.exists('fireworks/'):
        os.makedirs('fireworks/')
    content.export_fireworks(species='9606', ext='jpeg', file='fireworks_report', path='fireworks/', quality='10', 
                         flag=None, flag_interactors=False, sel=[], title=True, margin='15', resource='Total', 
                         diagram_profile='Calcium salts', coverage=False, token=token, exp_column=None)
    print('Fireworks report saved to fireworks/ folder.')


def getDownstreamPathways(token):
    '''This function takes a token and returns a list of downstream pathways.'''
    pathways = token_result['pathways']
    pathways_stId = [p['stId'] for p in pathways]
    fetch_downstream_pathways = [content.pathways_low_diagram(id=stId, species=None, all_forms=False) for stId in pathways_stId]
    downstream_pathways = [low_pathway for low_pathway in fetch_downstream_pathways if low_pathway is not None]
    has_diagram = [p[0]['hasDiagram'] for p in downstream_pathways]
    downstream_pathway_has_diagram = list(itertools.compress(downstream_pathways, has_diagram))
    downstream_pathway_stId = [p[0]['stId'] for p in downstream_pathway_has_diagram]

    if not os.path.exists('diagrams/'):
        os.makedirs('diagrams/')   
        [content.export_diagram(id=stId, ext='png', quality='5', flag_interactors=False, title=True, margin='15',
                            ehld=True, diagram_profile='Modern', resource='Total', analysis_profile='Standard', 
                            token=token, flag=None, sel=[], exp_column=None, file="-".join([stId,'report']), path='diagrams/') 
    for stId in downstream_pathway_stId]

