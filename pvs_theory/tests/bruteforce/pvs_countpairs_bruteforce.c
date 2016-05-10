/* File: countpairs.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "pvs_countpairs_bruteforce.h" //function proto-type
#include "utils.h" //all of the utilities

#ifndef SILENT
#include "progressbar.h" //for the progressbar
#endif

#if defined(USE_AVX) && defined(__AVX__)
#include "avx_calls.h"
#endif


void get_max_min(const int64_t ND1, const DOUBLE * restrict X1, const DOUBLE * restrict Y1, const DOUBLE * restrict Z1,
				 DOUBLE *min_x, DOUBLE *min_y, DOUBLE *min_z, DOUBLE *max_x, DOUBLE *max_y, DOUBLE *max_z)
{
  DOUBLE xmin = *min_x, ymin = *min_y, zmin=*min_z;
  DOUBLE xmax = *max_x, ymax = *max_y, zmax=*max_z;
	
  for(int64_t i=0;i<ND1;i++) {
    if(X1[i] < xmin) xmin=X1[i];
    if(Y1[i] < ymin) ymin=Y1[i];
    if(Z1[i] < zmin) zmin=Z1[i];


    if(X1[i] > xmax) xmax=X1[i];
    if(Y1[i] > ymax) ymax=Y1[i];
    if(Z1[i] > zmax) zmax=Z1[i];
  }
	*min_x=xmin;*min_y=ymin;*min_z=zmin;
	*max_x=xmax;*max_y=ymax;*max_z=zmax;
}

void free_results_pvs(results_pvs_countpairs **results)
{
    if(results == NULL)
        return;
    if(*results == NULL)
        return;

    results_pvs_countpairs *tmp = *results;

    free(tmp->rupp);
    free(tmp->npairs);
    free(tmp->rpavg);
    free(tmp->vpavg);
    free(tmp->spavg);
    free(tmp->vtavg);
    free(tmp->stavg);
    free(tmp);
    tmp = NULL;
}

results_pvs_countpairs * countpairs_pvs_bruteforce(const int64_t ND1, const DOUBLE * const X1, const DOUBLE * const Y1, const DOUBLE  * const Z1,const DOUBLE *VX1, const DOUBLE *VY1, const DOUBLE *VZ1,
                                const int64_t ND2, const DOUBLE * const X2, const DOUBLE * const Y2, const DOUBLE  * const Z2,const DOUBLE *VX2, const DOUBLE *VY2, const DOUBLE *VZ2,
                                const int autocorr,
                                const char *binfile)
{

    /***********************
     *initializing the  bins
     ************************/
    double *rupp;
    int nrpbin ;
    double rpmin,rpmax;
    setup_bins(binfile,&rpmin,&rpmax,&nrpbin,&rupp);
    assert(rpmin > 0.0 && rpmax > 0.0 && rpmin < rpmax && "[rpmin, rpmax] are valid inputs");
    assert(nrpbin > 0 && "Number of rp bins is valid");

    //Find the min/max of the data
    DOUBLE xmin,xmax,ymin,ymax,zmin,zmax;
    xmin=1e10;ymin=1e10;zmin=1e10;
    xmax=0.0;ymax=0.0;zmax=0.0;
    get_max_min(ND1, X1, Y1, Z1, &xmin, &ymin, &zmin, &xmax, &ymax, &zmax);

    if(autocorr==0) {
#ifndef SILENT
        fprintf(stderr,"ND1 = %12"PRId64" [xmin,ymin,zmin] = [%lf,%lf,%lf], [xmax,ymax,zmax] = [%lf,%lf,%lf]\n",ND1,xmin,ymin,zmin,xmax,ymax,zmax);
#endif
        get_max_min(ND2, X2, Y2, Z2, &xmin, &ymin, &zmin, &xmax, &ymax, &zmax);
#ifndef SILENT
        fprintf(stderr,"ND2 = %12"PRId64" [xmin,ymin,zmin] = [%lf,%lf,%lf], [xmax,ymax,zmax] = [%lf,%lf,%lf]\n",ND2,xmin,ymin,zmin,xmax,ymax,zmax);
#endif
    }

#ifndef SILENT
    fprintf(stderr,"Running with [xmin,xmax] = %lf,%lf\n",xmin,xmax);
    fprintf(stderr,"Running with [ymin,ymax] = %lf,%lf\n",ymin,ymax);
    fprintf(stderr,"Running with [zmin,zmax] = %lf,%lf\n",zmin,zmax);
#endif

    /*---Create 3-D lattice--------------------------------------*/
#ifdef PERIODIC
    const DOUBLE xdiff = (xmax-xmin);
    const DOUBLE ydiff = (ymax-ymin);
    const DOUBLE zdiff = (zmax-zmin);
#endif

    uint64_t npairs[nrpbin];
    for(int i=0; i < nrpbin;i++) npairs[i] = 0;

    DOUBLE rupp_sqr[nrpbin];
    for(int i=0; i < nrpbin;i++) {
        rupp_sqr[i] = rupp[i]*rupp[i];
    }


#ifdef OUTPUT_RPAVG
    DOUBLE rpavg[nrpbin];
#endif
    DOUBLE vpavg[nrpbin];
    DOUBLE spavg[nrpbin];
    DOUBLE vtavg[nrpbin];
    DOUBLE stavg[nrpbin];
#ifdef KAHN_SUM
    DOUBLE vpavg_c[nrpbin];
    DOUBLE spavg_c[nrpbin];
    DOUBLE vtavg_c[nrpbin];
    DOUBLE stavg_c[nrpbin];
#endif
    for(int i=0; i < nrpbin;i++) {
#ifdef OUTPUT_RPAVG
        rpavg[i] = 0.0;
#endif
        vpavg[i] = 0.0;
        spavg[i] = 0.0;
        vtavg[i] = 0.0;
        stavg[i] = 0.0;
#ifdef KAHN_SUM
        vpavg_c[i] = 0.0;
        spavg_c[i] = 0.0;
        vtavg_c[i] = 0.0;
        stavg_c[i] = 0.0;
#endif
    }

    DOUBLE sqr_rpmax=rupp_sqr[nrpbin-1];
    DOUBLE sqr_rpmin=rupp_sqr[0];

#ifndef SILENT    
    int interrupted=0;
    int64_t numdone=0;
    init_my_progressbar(ND1,&interrupted);
#endif    

    /*---Loop-over-Data1-particles--------------------*/
    for(int64_t index1=0;index1<ND1;index1++) {

          /* If the silent option is enabled, avoid outputting anything unnecessary*/
#ifndef SILENT          
      my_progressbar(numdone,&interrupted);
      numdone++;
#endif//SILENT

      DOUBLE x1pos = X1[index1];
      DOUBLE y1pos = Y1[index1];
      DOUBLE z1pos = Z1[index1];
      DOUBLE vx1pos = VX1[index1];
      DOUBLE vy1pos = VY1[index1];
      DOUBLE vz1pos = VZ1[index1];


    // *_countpairs_rp_pi double counts
      for (int64_t index2=0;index2 < ND2; index2++)	{ 

        DOUBLE dx = x1pos - X2[index2];
        DOUBLE dy = y1pos - Y2[index2];
        DOUBLE dz = z1pos - Z2[index2];
#ifdef PERIODIC
        if (dx > xdiff/2.) dx = dx - xdiff;
        else if (dx < -1*xdiff/2.) dx = dx + xdiff ;

        if (dy > ydiff/2.) dy = dy - ydiff;
        else if (dy < -1*ydiff/2.) dy = dy + ydiff ;

        if (dz > zdiff/2.) dz = dz - zdiff;
        else if (dz < -1*zdiff/2.) dz = dz + zdiff ;
#endif

        const DOUBLE r2 = (dx*dx + dy*dy + dz*dz);
        if(r2 >= sqr_rpmax || r2 < sqr_rpmin) {
          continue;
        }

        const DOUBLE dvx = vx1pos - VX2[index2];
        const DOUBLE dvy = vy1pos - VY2[index2];
        const DOUBLE dvz = vz1pos - VZ2[index2];

#ifdef OUTPUT_RPAVG
        const DOUBLE r = SQRT(r2);
        const DOUBLE vp = (dvx*dx + dvy*dy + dvz*dz)/r;
#else
        const DOUBLE vp = (dvx*dx + dvy*dy + dvz*dz)/SQRT(r2);
#endif
        const DOUBLE sp = vp*vp;
        const DOUBLE st = (dvx*dvx + dvy*dvy + dvz*dvz)-sp;
        const DOUBLE vt = SQRT(st);

        for(int kbin=nrpbin-1;kbin>=1;kbin--){
          if(r2 >= rupp_sqr[kbin-1]) {
            npairs[kbin]++;
#ifdef OUTPUT_RPAVG
            rpavg[kbin] += r;
#endif
            //vpavg[kbin]+=vp;
#ifndef KAHN_SUM
            vpavg[kbin]+=vp;
#else
            DOUBLE tmp1, tmp2;
            tmp1 = vp - vpavg_c[kbin];
            tmp2 = vpavg[kbin] + tmp1;
            vpavg_c[kbin] = (tmp2 - vpavg[kbin]) - tmp1;
            vpavg[kbin] = tmp2 ;
#endif
            //spavg[kbin]+=sp;
            //vtavg[kbin]+=vt;
            //stavg[kbin]+=st;
#ifndef KAHN_SUM
            spavg[kbin]+=sp;
#else
            tmp1 = sp - spavg_c[kbin];
            tmp2 = spavg[kbin] + tmp1;
            spavg_c[kbin] = (tmp2 - spavg[kbin]) - tmp1;
            spavg[kbin] = tmp2 ;
#endif


#ifndef KAHN_SUM
            vtavg[kbin]+=vt;
#else
            tmp1 = vt - vtavg_c[kbin];
            tmp2 = vtavg[kbin] + tmp1;
            vtavg_c[kbin] = (tmp2 - vtavg[kbin]) - tmp1;
            vtavg[kbin] = tmp2 ;
#endif

#ifndef KAHN_SUM
            stavg[kbin]+=st;
#else
            tmp1 = st - stavg_c[kbin];
            tmp2 = stavg[kbin] + tmp1;
            stavg_c[kbin] = (tmp2 - stavg[kbin]) - tmp1;
            stavg[kbin] = tmp2 ;
#endif
            break;
          }
        }//searching for kbin loop
      }//end of jj loop
    }


#ifndef SILENT    
    finish_myprogressbar(&interrupted);
#endif


    for(int i=0;i<nrpbin;i++) {
        if(npairs[i] > 0) {
#ifdef OUTPUT_RPAVG
            rpavg[i] /= (DOUBLE) npairs[i] ;
#endif
            vpavg[i] /= ((DOUBLE) npairs[i] );
            spavg[i] = SQRT(spavg[i]/((DOUBLE) npairs[i] ));
            vtavg[i] /= ((DOUBLE) npairs[i] );
            stavg[i] = SQRT(stavg[i]/((DOUBLE) npairs[i] ));
        }
    }

    //Pack in the results
    results_pvs_countpairs *results = my_malloc(sizeof(*results), 1);
    results->nbin = nrpbin;
    results->npairs = my_malloc(sizeof(uint64_t), nrpbin);
    results->rupp   = my_malloc(sizeof(DOUBLE)  , nrpbin);
    results->rpavg  = my_malloc(sizeof(DOUBLE)  , nrpbin);
    results->vpavg  = my_malloc(sizeof(DOUBLE)  , nrpbin);
    results->spavg  = my_malloc(sizeof(DOUBLE)  , nrpbin);
    results->vtavg  = my_malloc(sizeof(DOUBLE)  , nrpbin);
    results->stavg  = my_malloc(sizeof(DOUBLE)  , nrpbin);

    for(int i=0;i<nrpbin;i++) {
        results->npairs[i] = npairs[i];
        results->rupp[i] = rupp[i];
#ifdef OUTPUT_RPAVG
        results->rpavg[i] = rpavg[i];
#else
        results->rpavg[i] = 0.0;
#endif
        results->vpavg[i] = vpavg[i];
        results->spavg[i] = spavg[i];
        results->vtavg[i] = vtavg[i];
        results->stavg[i] = stavg[i];
    }

    free(rupp);

    return results;

}
