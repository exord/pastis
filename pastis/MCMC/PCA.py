###
# FUNCTIONS FOR PRINCIPAL COMPONENT ANALYSIS OF PARAMETERS
###
import numpy as n
import scipy


def get_covariance_matrix(C, bi=0, iteration=0):
    """
    Obtain the covariance matrix based on the paremeters of chain C
    and their eigenvectors, that will be used to convert the parameters
    to principal components.
    """

    if bi >= iteration:
        print('WARNING! BurnIn longer than current chain length. bi changed to 0')
        bi = 0

    # Construct array containing values of jumping parameters
    V = n.zeros((len(C.jumpind), iteration - bi))
    for ii, jj in enumerate(C.jumpind):
        V[ii] = C.values[bi: iteration, jj]

        
    # Subtract mean from each parameter, and trim burn-in period
    meanV = n.mean(V, axis = 1)
    VV = (V.T - meanV).T
    
    # Compute the covariance matrix S and its eigenvalues and vectors
    return n.cov(VV), meanV 


def run_pca(S, iteration = 0):
    """
    Obtain the eigenvalues and eigenvectors of covariance matrix S,
    which will be used to convert the parameters
    to principal components.
    """
    # Compute the covariance matrix S and its eigenvalues and vectors
    w, eigv = scipy.linalg.eigh( S )

    # Sort eigenvector matrix according to eigenvalue values
    ind = n.argsort(w)[::-1]
    coeff = eigv[:, ind]

    # Return transposition of coefficient matrix 
    return coeff.T, n.abs(w[ind])


def estimate_propscale_pca(C, M, S):
    """
    Estimate the size of the proposal scale for the Principal Components of the
    current state X of a given chain C, considering their actual values,
    the principal component coefficients and their covariance.
    """

    # Obtain jumpscales of current state X
    X = C.get_current_state()
    sx = []
    for ind in C.jumpind:
        sx.append(X[ind].propscale)

    sx = n.array(sx)
    sv = n.zeros(sx.shape)

    for i, jj in enumerate(C.jumpind):
        #sv[i] = n.sum((M[i]*sx)**2)
        sv[i] = 0.0
        for j in range(len(S)):
            for k in range(len(S)):
                #if j == k: continue
                sv[i] = sv[i] + M[i, j]*M[i, k]*S[j, k]
        #print n.sqrt(sv[i])

        #X[jj].PCApropscale = n.sqrt(sv[i])
    
    return n.sqrt(sv)
