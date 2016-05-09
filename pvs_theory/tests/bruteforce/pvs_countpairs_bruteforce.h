/* File: countpairs.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "function_precision.h" //for definition of DOUBLE
#include <inttypes.h> //for uint64_t


    //define the results structure
    typedef struct{
        uint64_t *npairs;
        DOUBLE *rupp;
        DOUBLE *rpavg;
        DOUBLE *vpavg;
        DOUBLE *spavg;
        DOUBLE *vtavg;
        DOUBLE *stavg;
        int nbin;
    } results_pvs_countpairs;

    results_pvs_countpairs * countpairs_pvs_bruteforce(const int64_t ND1, const DOUBLE * const X1, const DOUBLE * const Y1, const DOUBLE  * const Z1,const DOUBLE *VX1, const DOUBLE *VY1, const DOUBLE *VZ1,
                                    const int64_t ND2, const DOUBLE * const X2, const DOUBLE * const Y2, const DOUBLE  * const Z2,const DOUBLE *VX2, const DOUBLE *VY2, const DOUBLE *VZ2,
                                    const int autocorr,
                                    const char *binfile) __attribute__((warn_unused_result));


    void free_results_pvs(results_pvs_countpairs **results);

#ifdef __cplusplus
}
#endif
