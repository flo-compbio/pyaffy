# Copyright (c) 2016 Florian Wagner
#
# This file is part of pyAffy.
#
# pyAffy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

#cython: profile=False, wraparound=False, boundscheck=False, cdivision=True

cimport cython

import numpy as np
cimport numpy as np

np.import_array()

def medpolish(float[:,:] X, float eps = 0.01, int maxiter = 10, copy = True):

    if copy:
        X = X.copy()
    
    cdef int num_rows = X.shape[0]
    cdef int num_cols = X.shape[1]

    cdef float[::1] row_eff = np.zeros(num_rows, dtype = np.float32)
    cdef float[::1] col_eff = np.zeros(num_cols, dtype = np.float32)
    cdef float global_eff = 0.0
    cdef float sar = 0.0 # sum of absolute residuals
    
    cdef int converged = 0
    cdef int t = 0
    
    cdef float[::1] rdiff, cdiff
    cdef float diff
    cdef float new_sar
    cdef int i, j
    
    while (maxiter == -1 or t < maxiter) and converged == 0:
        
        # sweep rows
        rdiff = np.median(X, axis = 1)

        #X = X - np.tile(rdiff, (num_cols, 1)).T
        for i in range(num_rows):
            # substract median from the row
            for j in range(num_cols):
                X[i,j] -= rdiff[i]
            # add median to the row effect
            row_eff[i] += rdiff[i]
        # also the row containing the column effects
        diff = np.median(col_eff)
        for j in range(num_cols):
            col_eff[j] -= diff
        global_eff += diff

        # sweep columns
        cdiff = np.median(X, axis = 0)

        #X = X - np.tile(cdiff, (num_rows, 1))
        for j in range(num_cols):
            # substract median from the column
            for i in range(num_rows):
                X[i,j] -= cdiff[j]
            # add median to the column effect
            col_eff[j] += cdiff[j]
        # also the column containing the row effects
        diff = np.median(row_eff)
        for i in range(num_rows):
            row_eff[i] -= diff
        global_eff += diff

        # next
        t += 1

         # calculate new sum of absolute residuals
        new_sar = np.sum(np.absolute(X))
        # test for convergence
        if abs(new_sar - sar) < eps * new_sar:
            converged = 1

        sar = new_sar
        
    return X, np.float32(row_eff), np.float32(col_eff), float(global_eff), bool(converged), int(t)


def medpolish_missing(X, double eps = 0.01, int maxiter = 10, int copy = 1):

    if copy != 0:
        X = X.copy()
    
    cdef int num_rows = X.shape[0]
    cdef int num_cols = X.shape[1]

    cdef double[::1] row_eff = np.zeros(num_rows, dtype = np.float64)
    cdef double[::1] col_eff = np.zeros(num_cols, dtype = np.float64)
    cdef double global_eff = 0.0
    cdef double sar = 0.0 # sum of absolute residuals
    
    cdef int converged = 0
    cdef int t = 0
    
    cdef double[::1] rdiff, cdiff
    cdef double diff
    cdef double new_sar
    cdef int i, j
    
    while (maxiter == -1 or t < maxiter) and converged == 0:
        
        # sweep rows
        rdiff = np.ma.median(X, axis = 1).data # allow for masked values

        #X = X - np.tile(rdiff, (num_cols, 1)).T
        for i in range(num_rows):
            # substract median from the row
            for j in range(num_cols):
                X[i,j] -= rdiff[i]
            # add median to the row effect
            row_eff[i] += rdiff[i]
        # also the row containing the column effects
        diff = np.median(col_eff)
        for j in range(num_cols):
            col_eff[j] -= diff
        global_eff += diff

        # sweep columns
        cdiff = np.ma.median(X, axis = 0).data # allow for masked values

        #X = X - np.tile(cdiff, (num_rows, 1))
        for j in range(num_cols):
            # substract median from the column
            for i in range(num_rows):
                X[i,j] -= cdiff[j]
            # add median to the column effect
            col_eff[j] += cdiff[j]
        # also the column containing the row effects
        diff = np.median(row_eff)
        for i in range(num_rows):
            row_eff[i] -= diff
        global_eff += diff

        # next
        t += 1

         # calculate new sum of absolute residuals
        new_sar = np.sum(np.absolute(X))
        # test for convergence
        if abs(new_sar - sar) < eps * new_sar:
            converged = 1

        sar = new_sar
        
    return X, np.float64(row_eff), np.float64(col_eff), float(global_eff), bool(converged), int(t)
