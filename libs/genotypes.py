import scipy as sp

class genotypes(object):
    """
    Class representing genotypes
    """

    def __init__(self, SID = None, GT = None, POS = None):
        """
        Constructor
        SID: Sample ids
        GT: Genotype matrix
        POS: Positions
        """

        self.SID = SID
        self.GT = GT
        self.POS = POS

    def read_plink(fn_plink = None):
        """
        plink reader
        """

    def getK(self):
        """
        return Kinship matrix
        """
        
    def impute(self):
        """
        mean impute missing values
        """
        
     
