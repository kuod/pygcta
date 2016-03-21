import numpy as np
import numpy.linalg as la
import pdb

class pygcta(object):
    """
    class for pygcta
    """
    
    def __init__(self, Y = None, K = None, X = None):
        """
        Constructor
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
        
    def likelihood(self, V, P, Vinv=None, X=None):
        """
        implement likelihood function
        """
        if Vinv is None:
            Vinv = la.inv(V)
        if  X is None:
            X = self.X

        loglik = -0.5 * (np.log(la.det(V)) + np.log( la.det( np.dot(np.dot(X.T, Vinv), X) ) ) + np.dot( np.dot( self.Y.Y.T, P), self.Y.Y) )
        return loglik
        pass

    def getV(self, sigmai):
        """
        calculate V
        """
        V = np.zeros(self.K[0].K.shape )

        for i,x in enumerate(self.K):
            V += sigmai[i] * x.K
        V += np.eye(V.shape[0]) * sigmai[-1] ### multiply last one by identity
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

        Y = self.Y.Y
        N = Y.shape[0]

        sigma_next = ((sigmai ** 2 ) * np.dot(Y.T, np.dot(P, np.dot(A, np.dot(P, Y)))) + np.trace(sigmai * np.eye(N) - (sigmai ** 2) * np.dot(P,A))) / N
        
        return sigma_next.flatten()

    def optimize(self, tol = 1E-4):
        """
        MEAT: run optimization
        """
        N = self.X.shape[0]
        V0 = self.getV(self.sigma0)
        Vinv0 = la.inv(V0)
        P0 = self.getP(Vinv0)

        L_old = self.likelihood(V0, P0, Vinv0)
        print "Likelihood before EM step is: %f" % L_old

        # Update each variance component?? I have no idea what I am doing?
        sigma_next = np.zeros(self.sigma0.shape)
        for i in range(len(self.K)):
            sigma_next[i] = self.emstep(P0, self.K[i].K, self.sigma0[i])

        # How do we update remaining variance component?? Like this or with A=identity?
        sigma_next[-1] = self.emstep(P0, np.eye(N), self.sigma0[-1])

        V_next = self.getV(sigma_next)
        Vinv_next = la.inv(V_next)
        P_next = self.getP(Vinv_next)

        L_new = self.likelihood(V_next, P_next, Vinv_next)
        print "Likelihood after initial EM step is: %f" % L_new
        idx = 2
        while L_new - L_old > tol:
            L_old = L_new
            for i in range(len(self.K)):
                sigma_next[i] = self.emstep(P_next, self.K[i].K, sigma_next[i])
            sigma_next[-1] = self.emstep(P_next, np.eye(N), sigma_next[-1])

            V_next = self.getV(sigma_next)
            Vinv_next = la.inv(V_next)
            P_next = self.getP(Vinv_next)
            L_new = self.likelihood(V_next, P_next, Vinv_next)
            print "Likelihood after EM step %s is: %f (Likelihood Difference: +%e)" % (idx, L_new, L_new-L_old)
            idx += 1

        self.sigma = sigma_next
        print "Final Sigma is %f" % (sigma_next[1])

        
if __name__ == "__main__":
    pass
