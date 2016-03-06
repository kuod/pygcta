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
        PLOBST = PLOB.standardize()
        self.GT = PLOBST.val
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
	n_geno = self.GT.shape[1]
	n_sample = self.GT.shape[0]
	maf_vec = np.empty([n_geno,1])
	msf_vec = np.empty([n_geno,1])
	for i in range(n_geno):
	    na_geno = np.isnan(self.GT[:,i])
	    iOK = ~na_geno
	    na_count = sum(na_geno)
	    # calculate maf based on non nan 
	    msf_vec[i] = na_count / n_sample
	    maf_vec[i] = sum(self.GT[iOK,i]) / (n_sample * 2)
	# Filter missing genotypes first
	print "Filtering %d genotypes based on Missing Fraction\n" % (sum(msf_vec > msf))
	print "Excluding %d genotypes based on minor allele freq\n" % (sum(maf_vec <maf))
	keep_idx = (msf_vec < msf) * (maf_vec > maf)
	# can someone tell me why the ,dtype = boolean messes up the dimensions?
	# I had to add the [0] index because np.where(keep_idx) returned 2 dimensional arrays 
	self.GT = self.GT[:,np.where(keep_idx)[0]]
        
        
     
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
