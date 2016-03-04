import numpy as np
import pdb
from pysnptools.snpreader import Pheno

class phenotypes(object):
    """
    class representing phenotypes
    """

    def __init__(self, SID = None, Y = None, fn_phen = None):
        """
        Constructor
        SID: Sample id
        Y: phenotype
        """
        if not fn_phen is None:
            self.read_phen(fn_phen)
        else:
            self.SID = SID
            self.Y = Y

    def read_phen(self,fn_phen = None):
        """
        read phenotype file
        """
        PH = Pheno(fn_phen)
        PHOB = PH.read()
        self.Y = PHOB.val
        self.SID = PHOB.iid[:,1]


    def normalize(type=None):
        """
        normalize Y
        """
        pass

if __name__ == "__main__":
    fn_phen = "../data/test.phen"
    PHEN = phenotypes(fn_phen = fn_phen)

