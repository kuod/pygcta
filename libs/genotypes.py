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
        GT: Genotype matrix (Sample x Position)
        POS: Positions (Position x 2)
        """
        # Geno matrix was 3925 x 1000 
	# Pheno matrix was 3925 x 1 
	# I think the dimensions are incorrect for the code
	# Suggest that we transpose GT or fix code -J
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
        self.GT = PLOB.val
        self.POS = PLOB.pos[:,[0,1]]
        self.SID = PLOB.iid[:,1]
        self.isNormalised = False

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
        for i in xrange(self.GT.shape[1]):
            iOK = ~np.isnan(self.GT[:,i])
            mean = np.mean(self.GT[iOK,i])
            self.GT[:,i] -= mean
            Inan = np.isnan(self.GT[:,i])
            self.GT[Inan,i] = 0.0
            self.GT[iOK,i] /= np.std(self.GT[iOK,i])


    def filter(self, maf = None, msf = None):
        """
        filter genotypes
        maf = minimal allele frequency
        msf = minimal snp frequency across one genotype
        """

        if maf is None and msf is None:
            print "Please supply maf and/or msf for filtering."
            return
        elif self.isNormalised:
            print "Please filter prior to normalization."
            return

        maf_vec = np.asarray([np.min((af, 1-af)) for af in 0.5*np.nanmean(self.GT, 0)])
        msf_vec = np.mean(np.isnan(self.GT), 0)

        if maf is None:
            print "Excluding %d genotypes based on missing freq > %f" % (sum(msf_vec > msf), msf)
            keep_idx = (msf_vec < msf)
        elif msf is None:
            print "Excluding %d genotypes based on minor allele freq < %f" % (sum(maf_vec <maf), maf)
            keep_idx = (maf_vec > maf)
        else:
            print "Excluding %d genotypes based on minor allele freq < %f" % (sum(maf_vec <maf), maf)
            print "Excluding %d genotypes based on missing freq > %f" % (sum(msf_vec > msf), msf)
            keep_idx = (msf_vec < msf) * (maf_vec > maf)

        self.GT = self.GT[:,keep_idx.ravel()]
        
        
     
if __name__ == "__main__":

    ### Some unit testing
    #X = np.random.randn(20*5).reshape(20,5)
    #SID = np.array(['a','b','c','d','e'])
    #POS = np.array(zip(np.ones(20).tolist(),range(20)))

    ### create object
    GEN = genotypes()

    ### test read plink file
    fn_in = "../data/test.bed"
    GEN.read_plink(fn_plink = fn_in)
    GEN.filter(0.05, 0.1)


    ### test kinship
    #K = GEN.getK()
