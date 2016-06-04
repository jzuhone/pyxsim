import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport exp

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def broaden_lines(np.ndarray[np.float64_t, ndim=1] E0,
                  np.ndarray[np.float64_t, ndim=1] isigma,
                  np.ndarray[np.float64_t, ndim=1] amp,
                  np.ndarray[np.float64_t, ndim=1] E):

    cdef int i, j, n, m
    cdef double x
    cdef np.ndarray[np.float64_t, ndim=1] lines

    n = E0.shape[0]
    m = E.shape[0]
    lines = np.zeros(m)

    for i in range(n):
        for j in range(m):
            x = (E[j]-E0[i])*isigma[i]
            if x > -6.0 and x < 6.0:
                lines[j] += amp[i]*isigma[i]*exp(-x*x)

    return lines
