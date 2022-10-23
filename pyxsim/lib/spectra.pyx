import numpy as np
cimport numpy as np
cimport cython


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


def line_spectrum(
    int num_cells,
    np.ndarray[np.float64_t, ndim=1] ee,
    np.ndarray[np.float64_t, ndim=1] sigma,
    np.ndarray[np.float64_t, ndim=1] gx,
    np.ndarray[np.float64_t, ndim=1] gcdf,
    np.ndarray[np.float64_t, ndim=1] N,
    pbar
):
    cdef np.int64_t nbins = ee.size
    cdef np.ndarray[np.float64_t, ndim=1] spec
    cdef int i, j

    spec = np.zeros(nbins)

    for i in range(num_cells):
        xtmp = ee / sigma[i]
        for j in range(nbins):
            ret = np.interp(xtmp, gx, gcdf)
            spec[j] += N[i] * (ret[1:] - ret[:-1])
        pbar.update()
    return spec
