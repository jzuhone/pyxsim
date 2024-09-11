import numpy as np

cimport cython
cimport numpy as np


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cpow(True)
def power_law_spectrum(
    int num_cells,
    np.ndarray[np.float64_t, ndim=1] emid,
    np.ndarray[np.float64_t, ndim=1] alpha,
    np.ndarray[np.float64_t, ndim=1] K,
    pbar
):
    cdef np.int64_t nbins = emid.size
    cdef np.ndarray[np.float64_t, ndim=1] spec
    cdef int i, j

    spec = np.zeros(nbins)

    for i in range(num_cells):
        for j in range(nbins):
            spec[j] += K[i] * emid[j] ** -alpha[i]
        pbar.update()
    return spec


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def line_spectrum(
    int num_cells,
    np.ndarray[np.float64_t, ndim=1] ee,
    np.ndarray[np.float64_t, ndim=1] sigma,
    np.ndarray[np.float64_t, ndim=1] gx,
    np.ndarray[np.float64_t, ndim=1] gcdf,
    np.ndarray[np.float64_t, ndim=1] N,
    pbar
):
    cdef np.int64_t nbins = ee.size-1
    cdef np.ndarray[np.float64_t, ndim=1] spec, ret, xtmp
    cdef int i, j

    spec = np.zeros(nbins)
    ret = np.zeros(nbins+1)
    xtmp = np.zeros(nbins+1)

    for i in range(num_cells):
        for j in range(nbins):
            xtmp[j] = ee[j] / sigma[i]
        ret = np.interp(xtmp, gx, gcdf)
        for j in range(nbins):
            spec[j] += N[i] * (ret[j+1] - ret[j])
        pbar.update()
    return spec
        pbar.update()
    return spec
