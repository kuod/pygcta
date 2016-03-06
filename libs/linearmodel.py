import numpy as np
import numpy.linalg as la
import pdb

class pygcta(object):
    """
    class for pygcta
    """
    
    def __init__(self, Y = None, K = None, X = None):
        """
        Consturctor
        Y: Phenotype OBJECT
        K: LIST of kernels
        TODO: add covariates later
        """

        self.Y = Y
        self.K = K

        assert np.all(self.Y.SID == self.K[0].SID), 'Tough luck, your samples do not match so complete the matching function!!!!'

        self.betas = None
        self.pv = None
        self.sigma0 = np.array([np.var(self.Y.Y) / (len(self.K) + 1) for x in range(len(self.K) + 1)])### initialize
        if X is None:
            self.X = np.ones((K[0].K.shape[0],1))


    def matching(self, Y, K):
        """
        matching Y and K
        """

        ### sort data Y
        sidx = np.argsort(Y.SID)
        Y.SID = Y.SID[sidx]
        Y.Y = Y.Y[sidx,:]

        ### sort data K
        for kitem in K:
            sidx = np.argsort(kitem.SID)
            kitem.SID = kitem.SID[sidx]
            kitem.K = kitem.K[sidx,:][:,sidx]

        ### find reference set samples
        ### to boring to finish this crap....
        
    def likelihood(self, V, X, P):
        """
        implement likelihood function
        """
	# Why did you leave this commented when it was the correct likelihood equation?
        loglik = -0.5 * (np.log(np.det(V)) + np.log( np.det ( np.dot(np.dot(X.T, np.inv(V)), X) ) ) + np.dot( np.dot( self.Y.T, P), self.Y) )
        return loglik
        pass

    def getV(self, sigmai):
        """
        calculate V
        """
        V = np.zeros(self.K[0].K.shape )

        for i,x in enumerate(self.K):
            V += sigmai[i] * x.K
        V += np.eye(V.shape[0]) * sigmai[-1] ### multiply last one by idendity
        return V


    def getP(self, Vinv):
        """
        calculate P
        P = V^-1 - V^{-1}X(X'V^{-1}X)^{-1}X'V^{-1}
        """
        XVX = np.dot(np.dot(self.X.T, Vinv), self.X)
        P = Vinv - np.dot(np.dot(Vinv, self.X) * (1./XVX), np.dot(self.X.T, Vinv))
        return P
        
    
    def emstep(self, P, A, sigmai):
        """
        return emstep
        """
	# Assumed that n in the equation is sample size
	# Identity matrix should be N x N right?

	N = self.Y.shape()[0]
	# sigma_next = ((sigmai ** 2 ) * np.dot(self.Y.T, np.dot(P, np.dot(A, np.dot(P, self.Y)))) + np.trace(sigmai * np.identity(n = N) - (sigmai ** 2) * np.dot(P,A))) / N
        pass

    def optimize(self, tol = 1E4):
        """
        MEAT: run optimization
        """
        V0 = self.getV(self.sigma0)
        Vinv0 = la.inv(V0)
        P0 = self.getP(Vinv0)

	# while L_new - L_old > 1E4:
	    # continue optimization
        

        
if __name__ == "__main__":
    pass
