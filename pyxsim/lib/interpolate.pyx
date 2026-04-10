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
                    voutput[k, i, j] = vtable[k, x_i, j] * xm[i] + \
                                       vtable[k, x_i + 1, j] * xp[i]

    if do_var:
        return coutput, moutput, voutput
    else:
        return coutput, moutput, None


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def interp2d_spec(np.ndarray[np.float64_t, ndim=2] ctable,
                  np.ndarray[np.float64_t, ndim=2] mtable,
                  np.ndarray[np.float64_t, ndim=3] vtable,
                  np.ndarray[np.float64_t, ndim=1] x_vals,
                  np.ndarray[np.float64_t, ndim=1] x_bins,
                  np.ndarray[np.int32_t, ndim=1] x_is,
                  np.ndarray[np.float64_t, ndim=1] y_vals,
                  np.ndarray[np.float64_t, ndim=1] y_bins,
                  np.ndarray[np.int32_t, ndim=1] y_is,
                  bint do_var):
    cdef double x, dx_inv, y, dy_inv
    cdef int i, x_i, j, k, nelem, y_i
    cdef int z1, z2, z3, z4
    cdef int nt = x_vals.shape[0]
    cdef int ne = ctable.shape[1]
    cdef int ntbins = x_bins.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] xp, xm, yp, ym
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
    yp = np.zeros(nt)
    ym = np.zeros(nt)

    for i in range(nt):
        x_i = x_is[i]
        x = x_vals[i]
        dx_inv = 1.0 / (x_bins[x_i+1] - x_bins[x_i])
        xp[i] = (x - x_bins[x_i]) * dx_inv
        xm[i] = (x_bins[x_i+1] - x) * dx_inv
        y_i = y_is[i]
        y = y_vals[i]
        dy_inv = 1.0 / (y_bins[y_i+1] - y_bins[y_i])
        yp[i] = (y - y_bins[y_i]) * dy_inv
        ym[i] = (y_bins[y_i+1] - y) * dy_inv

    for i in range(nt):
        x_i = x_is[i]
        y_i = y_is[i]
        z1 = ntbins*y_i + x_i
        z2 = z1 + ntbins
        z3 = z1 + 1
        z4 = z2 + 1
        for j in range(ne):
            coutput[i, j] = ctable[z1, j] * xm[i] * ym[i] + \
                            ctable[z2, j] * xm[i] * yp[i] + \
                            ctable[z3, j] * xp[i] * ym[i] + \
                            ctable[z4, j] * xp[i] * yp[i]
            moutput[i, j] = mtable[z1, j] * xm[i] * ym[i] + \
                            mtable[z2, j] * xm[i] * yp[i] + \
                            mtable[z3, j] * xp[i] * ym[i] + \
                            mtable[z4, j] * xp[i] * yp[i]
    if do_var:
        for i in range(nt):
            x_i = x_is[i]
            y_i = y_is[i]
            z1 = ntbins * y_i + x_i
            z2 = z1 + ntbins
            z3 = z1 + 1
            z4 = z2 + 1
            for k in range(nelem):
                for j in range(ne):
                    voutput[k, i, j] = vtable[k, z1, j] * xm[i] * ym[i] + \
                                       vtable[k, z2, j] * xm[i] * yp[i] + \
                                       vtable[k, z3, j] * xp[i] * ym[i] + \
                                       vtable[k, z4, j] * xp[i] * yp[i]

    if do_var:
        return coutput, moutput, voutput
    else:
        return coutput, moutput, None
