import numpy as np
cimport numpy as np
cimport cython

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def scatter_events(normal, prng, kernel, data_type,
                   int num_det,
                   np.ndarray[np.int64_t, ndim=1] idxs,
                   np.ndarray[np.int64_t, ndim=1] n_ph,
                   np.ndarray[np.float64_t, ndim=2] pos,
                   np.ndarray[np.float64_t, ndim=1] dx,
                   np.ndarray[np.float64_t, ndim=1] x_hat,
                   np.ndarray[np.float64_t, ndim=1] y_hat):

    cdef int num_cells = idxs.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] xsky, ysky, zsky
    cdef int i, j, k, xax, yax, idx

    k = 0

    if isinstance(normal, int):
        
        if normal == 0:
            xax = 1
            yax = 2
        elif normal == 1:
            xax = 2
            yax = 0
        elif normal == 2:
            xax = 0
            yax = 1
    
        if data_type == "cells":
            xsky = prng.uniform(low=-0.5, high=0.5, size=num_det)
            ysky = prng.uniform(low=-0.5, high=0.5, size=num_det)
        elif data_type == "particles":
            if kernel == "gaussian":
                xsky = prng.normal(loc=0.0, scale=1.0, size=num_det)
                ysky = prng.normal(loc=0.0, scale=1.0, size=num_det)
            elif kernel == "top_hat":
                r = prng.uniform(low=0.0, high=1.0, size=num_det)
                theta = 2.0*np.pi*prng.uniform(low=0.0, high=1.0, size=num_det)
                xsky = r*np.cos(theta)
                ysky = r*np.sin(theta)
    
        for i in range(num_det):
            idx = idxs[i]
            for j in range(n_ph[idx]):
                xsky[k] = xsky[k]*dx[idx] + pos[idx, xax]
                ysky[k] = ysky[k]*dx[idx] + pos[idx, yax]
                k += 1

    else:
    
        if data_type == "cells":
            xsky, ysky, zsky = prng.uniform(low=-0.5, high=0.5, 
                                            size=(3, num_det))
            for i in range(num_det):
                idx = idxs[i]
                for j in range(n_ph[idx]):
                    xsky[k] = xsky[k]*dx[idx] + pos[idx, 0]
                    ysky[k] = ysky[k]*dx[idx] + pos[idx, 1]
                    zsky[k] = zsky[k]*dx[idx] + pos[idx, 2]
                    xsky[k] = (xsky[k]*x_hat[0]+
                               ysky[k]*x_hat[1]+
                               zsky[k]*x_hat[2])
                    ysky[k] = (xsky[k]*y_hat[0]+
                               ysky[k]*y_hat[1]+
                               zsky[k]*y_hat[2])
                    k += 1
        elif data_type  == "particles":
            if kernel == "gaussian":
                xsky = prng.normal(loc=0.0, scale=1.0, size=num_det)
                ysky = prng.normal(loc=0.0, scale=1.0, size=num_det)
            elif kernel == "top_hat":
                r = prng.uniform(low=0.0, high=1.0, size=num_det)
                theta = 2.0*np.pi*prng.uniform(low=0.0, high=1.0, size=num_det)
                xsky = r*np.cos(theta)
                ysky = r*np.sin(theta)
            for i in range(num_det):
                idx = idxs[i]
                for j in range(n_ph[idx]):
                    xsky[k] = xsky[k]*dx[idx] + pos[idx, 0]*x_hat[0] + \
                        pos[idx, 1]*x_hat[1] + pos[idx, 2]*x_hat[2]
                    ysky[k] = ysky[k]*dx[idx] + pos[idx, 0]*y_hat[0] + \
                        pos[idx, 1]*y_hat[1] + pos[idx, 2]*y_hat[2]
                    k += 1
    
    return xsky, ysky
