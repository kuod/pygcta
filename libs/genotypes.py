import numpy as np
import pdb
import kernel
from pysnptools.snpreader import Bed

class genotypes(object):
    """
    Class representing genotypes
    """

    def __init__(self, SID = None, GT = None, POS = None, isNormalised = False):
        """
        Constructor
        SID: Sample ids
        GT: Genotype matrix (Position x Sample)
        POS: Positions (Position x 2)
        """

        self.SID = SID
        self.GT = GT
        self.POS = POS
        self.isNormalised = isNormalised

    def read_plink(self, fn_plink = None):
        """
        plink reader
        """
        PL = Bed(fn_plink)
        PLOB = PL.read()
        PLOBST = PLOB.standardize()
        self.GT = PLOBST.val.T
        self.POS = PLOB.pos[:,[0,1]]
        self.SID = PLOB.iid[:,1]
        self.isNormalised = True



    def getK(self):
        """
        return Kinship matrix as kernel object
        """
        if not self.isNormalised:
            self.impute()
        K = np.dot(self.GT, self.GT.T)
        K /= K.diagonal().mean()
        Kobj = kernel.kernel(self.SID, K)
        return Kobj
        
    def impute(self):
        """
        mean impute missing values
        and normalize
        """
        for i in xrange(self.GT.shape[0]):
            iOK = ~np.isnan(self.GT[i,:])
            mean = np.mean(self.GT[i,iOK])
            self.GT[i,:] -= mean
            Inan = np.isnan(self.GT[i,:])
            self.GT[i,Inan] = 0.0
            self.GT[i,iOK] /= np.std(self.GT[i,iOK])


    def filter(self, maf = 0.05, msf = 0.5):
        """
        filter genotypes
        maf = minimal allele frequency
        msf = minimal snp frequency across one genotype
        """
        pass
        
     
if __name__ == "__main__":

    ### Some unit testing
    X = np.random.randn(20*5).reshape(20,5)
    SID = np.array(['a','b','c','d','e'])
    POS = np.array(zip(np.ones(20).tolist(),range(20)))

    

    ### create object
    GEN = genotypes(SID = SID, GT = X, POS = POS)

    ### test read plink file
    fn_in = "../data/test.bed"
    GEN.read_plink(fn_plink = fn_in)


    ### test kinship
    K = GEN.getK()
