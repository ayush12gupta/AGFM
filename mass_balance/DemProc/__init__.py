#
# Author: Piyush Agram
# Copyright 2016
#

from .DemProc import *
from .Factories import *

def getFactoriesInfo():
    return  {'TopsProc':
                     {'args':
                           {
                            'procDoc':{'value':None,'type':'Catalog','optional':True}
                            },
                     'factory':'createDemProc'                     
                     }
              
              }

def createDemProc(name=None, procDoc= None):
    from .DemProc import DemProc
    return DemProc(name = name,procDoc = procDoc)
