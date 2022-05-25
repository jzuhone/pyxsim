import numpy as np
cimport numpy as np
cimport cython

cdef extern from "math.h":
    double sqrt(double x) nogil
    double sin(double x) nogil
    double cos(double x) nogil
    double atan(double x) nogil
    double atan2(double y, double x) nogil
    double asin(double x) nogil
    double acos(double x) nogil

    
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def doppler_shift(np.ndarray[np.float64_t, ndim=1] shift,
                  np.ndarray[np.int64_t, ndim=1] n_ph,
                  np.ndarray[np.float64_t, ndim=1] eobs):

    cdef np.int64_t num_cells = n_ph.shape[0]
    cdef np.float64_t shft
    cdef np.int64_t i, j, k

    k = 0
    for i in range(num_cells):
        shft = sqrt((1.-shift[i])/(1.+shift[i]))
        for j in range(n_ph[i]):
            eobs[k] = eobs[k] * shft
            k += 1

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def scatter_events(normal, prng, kernel, data_type,
                   int num_det,
                   np.ndarray[np.uint8_t, cast=True] det,
                   np.ndarray[np.int64_t, ndim=1] n_ph,
                   np.ndarray[np.float64_t, ndim=1] x,
                   np.ndarray[np.float64_t, ndim=1] y,
                   np.ndarray[np.float64_t, ndim=1] z,
                   np.ndarray[np.float64_t, ndim=1] dx,
                   np.ndarray[np.float64_t, ndim=1] x_hat,
                   np.ndarray[np.float64_t, ndim=1] y_hat):

    cdef np.int64_t num_cells = dx.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] xsky, ysky, zsky, r, theta
    cdef np.int64_t i, j, k, xax, yax, n
    cdef np.float64_t xx, yy

    k = 0
    n = 0

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
    
        for i in range(num_cells):
            for j in range(n_ph[i]):
                if det[n]:
                    xsky[k] *= dx[i]
                    ysky[k] *= dx[i]
                    if normal == 0:
                        xsky[k] += y[i]
                        ysky[k] += z[i]
                    elif normal == 1:
                        xsky[k] += z[i]
                        ysky[k] += x[i]
                    elif normal == 2:
                        xsky[k] += x[i]
                        ysky[k] += y[i]
                    k += 1
                n += 1

    else:
    
        if data_type == "cells":
            xsky, ysky, zsky = prng.uniform(low=-0.5, high=0.5, 
                                            size=(3, num_det))
            for i in range(num_cells):
                for j in range(n_ph[i]):
                    if det[n]:
                        xsky[k] = xsky[k]*dx[i] + x[i]
                        ysky[k] = ysky[k]*dx[i] + y[i]
                        zsky[k] = zsky[k]*dx[i] + z[i]
                        xx = (xsky[k]*x_hat[0]+ ysky[k]*x_hat[1]+
                              zsky[k]*x_hat[2])
                        yy = (xsky[k]*y_hat[0]+ ysky[k]*y_hat[1]+
                              zsky[k]*y_hat[2])
                        xsky[k] = xx
                        ysky[k] = yy
                        k += 1
                    n += 1
        elif data_type  == "particles":
            if kernel == "gaussian":
                xsky = prng.normal(loc=0.0, scale=1.0, size=num_det)
                ysky = prng.normal(loc=0.0, scale=1.0, size=num_det)
            elif kernel == "top_hat":
                r = prng.uniform(low=0.0, high=1.0, size=num_det)
                theta = 2.0*np.pi*prng.uniform(low=0.0, high=1.0, size=num_det)
                xsky = r*np.cos(theta)
                ysky = r*np.sin(theta)
            for i in range(num_cells):
                for j in range(n_ph[i]):
                    if det[n]:
                        xsky[k] = xsky[k]*dx[i] + x[i]*x_hat[0] + \
                            y[i]*x_hat[1] + z[i]*x_hat[2]
                        ysky[k] = ysky[k]*dx[i] + x[i]*y_hat[0] + \
                            y[i]*y_hat[1] + z[i]*y_hat[2]
                        k += 1
                    n += 1
    
    return xsky, ysky


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def pixel_to_cel(np.ndarray[np.float64_t, ndim=1] xsky,
                 np.ndarray[np.float64_t, ndim=1] ysky,
                 np.ndarray[np.float64_t, ndim=1] sky_center):

    cdef int i
    cdef int n = xsky.size
    cdef np.float64_t B, D, cx, cy, sin_cy, cos_cy
    cdef np.float64_t PI = np.pi

    cx = sky_center[0]*PI/180.0
    cy = sky_center[1]*PI/180.0
    sin_cy = sin(cy)
    cos_cy = cos(cy)

    for i in range(n):
        
        D = atan(sqrt(xsky[i]*xsky[i] + ysky[i]*ysky[i]))
        B = atan2(-xsky[i], -ysky[i])

        xsky[i] = sin_cy*sin(D)*cos(B) + cos_cy*cos(D)
        ysky[i] = sin(D)*sin(B)

        xsky[i] = cx + atan2(ysky[i], xsky[i])
        ysky[i] = asin(sin_cy*cos(D) - cos_cy*sin(D)*cos(B))

        xsky[i] *= 180.0/PI
        ysky[i] *= 180.0/PI


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def scatter_events_allsky(data_type, kernel, prng, int num_det,
                          np.ndarray[np.uint8_t, cast=True] det,
                          np.ndarray[np.int64_t, ndim=1] n_ph,
                          np.ndarray[np.float64_t, ndim=1] x,
                          np.ndarray[np.float64_t, ndim=1] y,
                          np.ndarray[np.float64_t, ndim=1] z,
                          np.ndarray[np.float64_t, ndim=1] dx,
                          np.ndarray[np.float64_t, ndim=1] x_hat,
                          np.ndarray[np.float64_t, ndim=1] y_hat,
                          np.ndarray[np.float64_t, ndim=1] z_hat):
    cdef np.int64_t num_cells = x.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] theta, phi, lat, lon
    cdef np.int64_t i, j, k, n
    cdef np.float64_t xx, yy, zz
    cdef np.float64_t PI = np.pi

    k = 0
    n = 0

    lat = np.zeros(num_det)
    lon = np.zeros(num_det)

    if data_type == "cells":
        xsky, ysky, zsky = prng.uniform(low=-0.5, high=0.5,
                                        size=(3, num_det))
    else:
        if kernel == "gaussian":
            xsky = prng.normal(loc=0.0, scale=1.0, size=num_det)
            ysky = prng.normal(loc=0.0, scale=1.0, size=num_det)
            zsky = prng.normal(loc=0.0, scale=1.0, size=num_det)
        elif kernel == "top_hat":
            r = prng.uniform(low=0.0, high=1.0, size=num_det)
            phi = 2.0 * np.pi * prng.uniform(low=0.0, high=1.0, size=num_det)
            theta = np.arccos(prng.uniform(low=-1., high=1., size=num_det))
            xsky = r * np.cos(phi)*np.sin(theta)
            ysky = r * np.sin(phi)*np.sin(theta)
            zsky = r * np.cos(theta)
    for i in range(num_cells):
        if n_ph[i] == 0:
            continue
        for j in range(n_ph[i]):
            if det[n]:
                xsky[k] = xsky[k] * dx[i] + x[i]
                ysky[k] = ysky[k] * dx[i] + y[i]
                zsky[k] = zsky[k] * dx[i] + z[i]
                xx = (xsky[k] * x_hat[0] + ysky[k] * x_hat[1] +
                      zsky[k] * x_hat[2])
                yy = (xsky[k] * y_hat[0] + ysky[k] * y_hat[1] +
                      zsky[k] * y_hat[2])
                zz = (xsky[k] * z_hat[0] + ysky[k] * z_hat[1] +
                      zsky[k] * z_hat[2])
                rr = sqrt(xx * xx + yy * yy + zz * zz)
                lon[k] = atan2(yy, xx)
                if lon[k] < 0.0:
                    lon[k] += 2.0*PI
                lat[k] = 0.5*PI - acos(zz/rr)
                lon[k] *= 180.0 / PI
                lat[k] *= 180.0 / PI
                k += 1
            n += 1

    return lon, lat
