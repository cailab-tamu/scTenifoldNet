import numpy as np
def nonLinearManifold(X,Y,corr,num_dims,Wx,Wy,mu=0.9,eps=1e-8):
  Wxy = corr
  Wxy = mu * (Wx.sum() + Wy.sum()) / (2 * Wxy.sum()) * Wxy
  W = np.asarray(np.bmat(((Wx,Wxy),(Wxy.T,Wy))))
  n_nodes = W.shape[0]
  lap = -np.asarray(W)  # minus sign leads to a copy
  # set diagonal to zero, in case it isn't already
  lap.flat[::n_nodes + 1] = 0
  d = -lap.sum(axis=0)  # re-negate to get positive degrees
  # put the degrees on the diagonal
  lap.flat[::n_nodes + 1] = d
  L = lap
  d1 = X.shape[0]
  d2 = Y.shape[0]
  vals,vecs = np.linalg.eig(L)
  idx = np.argsort(vals)
  for i in range(len(idx)):
    if vals[idx[i]] >= eps:
      break
  vecs = vecs.real[:,idx[i:]]
  for i in range(vecs.shape[1]):
    vecs[:,i] /= np.linalg.norm(vecs[:,i])
  return vecs[:,:num_dims]
