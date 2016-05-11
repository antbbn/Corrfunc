/* File: pvs_countpairs_rp_pi_bruteforce.h */
/*
		This file is part of the PVS extension to the Corrfunc package
		Copyright (C) 2016-- Antonio Bibiano (antbbn@gmail.com)
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

#ifndef DOUBLE_PREC
#define DOUBLE_PREC
#endif

#ifndef OUTPUT_RPAVG
#define OUTPUT_RPAVG
#endif

#ifndef KAHN_SUM
#define KAHN_SUM
#endif
//define the results structure
typedef struct{
	uint64_t *npairs;
	DOUBLE *rupp;
	DOUBLE *rpavg;
        DOUBLE *vpavg;
        DOUBLE *spavg;
        DOUBLE *vtavg;
        DOUBLE *stavg;
	DOUBLE pimax;
	int nbin;
	int npibin;
} results_pvs_countpairs_rp_pi;

results_pvs_countpairs_rp_pi * bruteforce_pvs_countpairs_rp_pi(const int64_t ND1, const DOUBLE *X1, const DOUBLE *Y1, const DOUBLE *Z1,const DOUBLE *VX1, const DOUBLE *VY1, const DOUBLE *VZ1,
											const int64_t ND2, const DOUBLE *X2, const DOUBLE *Y2, const DOUBLE *Z2,const DOUBLE *VX2, const DOUBLE *VY2, const DOUBLE *VZ2,
											const int autocorr,
											const char *binfile,
											const double pimax)  __attribute__((warn_unused_result));

void free_results_pvs_rp_pi(results_pvs_countpairs_rp_pi **results);
void get_max_min(const int64_t ND1, const DOUBLE * restrict X1, const DOUBLE * restrict Y1, const DOUBLE * restrict Z1,
				 DOUBLE *min_x, DOUBLE *min_y, DOUBLE *min_z, DOUBLE *max_x, DOUBLE *max_y, DOUBLE *max_z);

#ifdef __cplusplus
}
#endif
