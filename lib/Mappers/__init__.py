#!/usr/bin/env python


import os
import sys


##updating the library look up path
dir = os.path.dirname(__file__)
sys.path.append(os.path.join(dir,'../lib'));

#load the modules
from BWA import BWA
#from Mappers.BWA import BWA 