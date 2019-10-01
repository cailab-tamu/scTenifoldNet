import numpy as np
import scipy as sp
import scipy.linalg
from functools import reduce

def LinearManifold(X, Y, corr, num_dims, Wx, Wy, mu=0.9, eps=1e-8):
    Wxy = corr
    Wxy = mu * (Wx.sum() + Wy.sum()) / (2 * Wxy.sum()) * Wxy
    W = np.asarray(np.bmat(((Wx,Wxy),(Wxy.T,Wy))))
    # Laplacian function
    n_nodes = W.shape[0]
    lap = -np.asarray(W)  # minus sign leads to a copy
    # set diagonal to zero, in case it isn't already
    lap.flat[::n_nodes + 1] = 0
    d = -lap.sum(axis=0)  # re-negate to get positive degrees
    # put the degrees on the diagonal
    lap.flat[::n_nodes + 1] = d
    L = lap
    # LinearDecompose function
    Z = sp.linalg.block_diag(X.T,Y.T)
    u,s,_ = np.linalg.svd(np.dot(Z,Z.T))
    Fplus = np.linalg.pinv(np.dot(u,np.diag(np.sqrt(s))))
    T = reduce(np.dot,(Fplus,Z,L,Z.T,Fplus.T))
    L = 0.5*(T+T.T)
    d1,d2 = X.shape[1],Y.shape[1]
    # ManifoldDecompose function
    vals,vecs = np.linalg.eig(L)
    idx = np.argsort(vals)
    for i in range(len(idx)):
      if vals[idx[i]] >= eps:
        break
    vecs = vecs.real[:,idx[i:]]
    vecs = np.dot(Fplus.T,vecs)
    for i in range(vecs.shape[1]):
      vecs[:,i] /= np.linalg.norm(vecs[:,i])
    return(vecs[:,:num_dims])
