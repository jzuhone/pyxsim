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
    np.ndarray[np.float64_t, ndim=1] shift,
    pbar
):
    cdef np.int64_t nbins = emid.size
    cdef np.ndarray[np.float64_t, ndim=1] spec
    cdef int i, j

    spec = np.zeros(nbins)

    for i in range(num_cells):
        for j in range(nbins):
            spec[j] += shift[i] * shift[i] * K[i] * (emid[j] / shift[i]) ** -alpha[i]
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
    np.ndarray[np.float64_t, ndim=1] shift,
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
            spec[j] += shift[i] * shift[i] * N[i] * (ret[j+1] - ret[j])
        pbar.update()
    return spec


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def shift_spectrum(
    np.ndarray[np.float64_t, ndim=1] ebold,
    np.ndarray[np.float64_t, ndim=1] ebnew,
    np.ndarray[np.float64_t, ndim=2] spec,
    np.ndarray[np.float64_t, ndim=1] shift,
    np.ndarray[np.float64_t, ndim=1] N,
):
    cdef np.int64_t nold = ebold.size-1
    cdef np.int64_t nnew = ebnew.size-1
    cdef np.int64_t num_cells = N.size
    cdef np.ndarray[np.float64_t, ndim=1] nspec, cspec
    cdef int i, j

    nspec = np.zeros(nnew+1)
    cspec = np.zeros(nold+1)
    ospec = np.zeros(nnew)

    for i in range(num_cells):
        for j in range(1, nold+1):
            cspec[j] = spec[i, j-1] + cspec[j-1]
        nspec = np.interp(ebnew/shift[i], ebold, cspec)
        for j in range(nnew):
            ospec[j] += shift[i] * shift[i] * N[i] * (nspec[j+1] - nspec[j])
    return ospec


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def make_band(
    int use_energy,
    double emin,
    double emax,
    np.ndarray[np.float64_t, ndim=1] ebins,
    np.ndarray[np.float64_t, ndim=1] emid,
    np.ndarray[np.float64_t, ndim=2] spec,
    np.ndarray[np.float64_t, ndim=1] shift,
):
    cdef np.int64_t nbins = emid.size
    cdef np.int64_t num_cells = spec.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] ofield
    cdef int i, j
    cdef double spec_sum, shft, ener

    ofield = np.zeros(num_cells)

    for i in range(num_cells):
        spec_sum = 0.0
        for j in range(nbins):
            if use_energy == 1:
                ener = emid[j]
            else:
                ener = 1.0
            if ebins[j] >= emin/shift[i] and ebins[j+1] <= emax/shift[i]:
                spec_sum += ener*spec[i, j]
        shft = shift[i]*shift[i]
        if use_energy == 1:
            shft *= shift[i]
        ofield[i] += shft * spec_sum
    return ofield
