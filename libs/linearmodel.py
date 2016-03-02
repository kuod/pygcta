import scipy sa sp

class pygcta(object):
    """
    class for pygcta
    """
    
    def __init__(self, Y = None, K = None):
        """
        Consturctor
        Y: Phenotype OBJECT
        K: LIST of kernels
        TODO: add covariates later
        """

        self.Y = Y
        self.K = K



        self.betas = None
        self.pv = None
        self.P = None ### projection matrix
        self.V = None ### covariance matrix
        
    def likelihood():
        """
        implement likelihood function
        """
        pass
    
    def emstep():
        """
        return emstep
        """
        pass

    def optimize():
        """
        MEAT: run optimization
        """
        pass
        
        
