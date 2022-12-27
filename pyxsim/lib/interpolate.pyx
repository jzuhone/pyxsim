import numpy as np

cimport cython
cimport numpy as np


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def interp1d_spec(np.ndarray[np.float64_t, ndim=2] ctable,
                  np.ndarray[np.float64_t, ndim=2] mtable,
                  np.ndarray[np.float64_t, ndim=3] vtable,
                  np.ndarray[np.float64_t, ndim=1] x_vals,
                  np.ndarray[np.float64_t, ndim=1] x_bins,
                  np.ndarray[np.int32_t, ndim=1] x_is,
                  bint do_var):
    cdef double x, dx_inv
    cdef int i, x_i, j, k, nelem
    cdef int nt = x_vals.shape[0]
    cdef int ne = ctable.shape[1]
    cdef np.ndarray[np.float64_t, ndim=1] xp, xm
    cdef np.ndarray[np.float64_t, ndim=2] coutput, moutput
    cdef np.ndarray[np.float64_t, ndim=3] voutput

    coutput = np.zeros((nt, ne))
    moutput = np.zeros((nt, ne))
    if do_var:
        nelem = <int>vtable.shape[0]
        voutput = np.zeros((nelem, nt, ne))
    else:
        nelem = 0

    xp = np.zeros(nt)
    xm = np.zeros(nt)

    for i in range(nt):
        x_i = x_is[i]
        x = x_vals[i]
        dx_inv = 1.0 / (x_bins[x_i+1] - x_bins[x_i])
        xp[i] = (x - x_bins[x_i]) * dx_inv
        xm[i] = (x_bins[x_i+1] - x) * dx_inv
    for i in range(nt):
        x_i = x_is[i]
        for j in range(ne):
            coutput[i, j] = ctable[x_i, j] * xm[i] + ctable[x_i+1, j] * xp[i]
            moutput[i, j] = mtable[x_i, j] * xm[i] + mtable[x_i + 1, j] * xp[i]
    if do_var:
        for i in range(nt):
            x_i = x_is[i]
            for k in range(nelem):
                for j in range(ne):
                    voutput[k, i, j] = vtable[k, x_i, j] * xm[i] + vtable[k, x_i + 1, j] * xp[i]

    if do_var:
        return coutput, moutput, voutput
    else:
        return coutput, moutput, None
