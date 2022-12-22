import numpy as np

cimport cython
cimport numpy as np


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def interp1d_spec(np.ndarray[np.float64_t, ndim=2] table,
                  np.ndarray[np.float64_t, ndim=1] x_vals,
                  np.ndarray[np.float64_t, ndim=1] x_bins,
                  np.ndarray[np.int32_t, ndim=1] x_is):
    cdef double x, xp, xm
    cdef int i, x_i, j
    cdef int nt = x_vals.shape[0]
    cdef int ne = table.shape[1]
    cdef np.ndarray[np.float64_t, ndim=2] output

    output = np.zeros((nt, ne))

    for i in range(nt):
        x_i = x_is[i]
        x = x_vals[i]
        dx_inv = 1.0 / (x_bins[x_i+1] - x_bins[x_i])
        xp = (x - x_bins[x_i]) * dx_inv
        xm = (x_bins[x_i+1] - x) * dx_inv
        for j in range(ne):
            output[i,j] = table[x_i, j] * xm + table[x_i+1, j] * xp
    return output


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def interp1d_var_spec(np.ndarray[np.float64_t, ndim=3] table,
                      np.ndarray[np.float64_t, ndim=1] x_vals,
                      np.ndarray[np.float64_t, ndim=1] x_bins,
                      np.ndarray[np.int32_t, ndim=1] x_is):
    cdef double x, xp, xm
    cdef int i, x_i, j, k
    cdef int ne = table.shape[2]
    cdef int nelem = table.shape[0]
    cdef int nt = x_vals.shape[0]
    cdef np.ndarray[np.float64_t, ndim=3] output

    output = np.zeros((nelem, nt, ne))

    for i in range(nt):
        x_i = x_is[i]
        x = x_vals[i]
        dx_inv = 1.0 / (x_bins[x_i+1] - x_bins[x_i])
        xp = (x - x_bins[x_i]) * dx_inv
        xm = (x_bins[x_i+1] - x) * dx_inv
        for k in range(nelem):
            for j in range(ne):
                output[k,i,j] = table[k, x_i, j] * xm + table[k, x_i+1, j] * xp
    return output
