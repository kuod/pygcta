import numpy as np
import pdb


class kernel(object):
    """
    kernel class
    """

    def __init__(self, SID = None, K = None):
        
        self.SID = SID
        self.K = K


    def save(self, fn):
        sp.savez(fn, self.K, self.SID)

    def load(self, fn):
        data = sp.load(fn)
        self.K = data[data.keys()[0]]
        self.SID = data[data.keys()[1]]



    
        
