/* File: pvs_countpairs_rp_pi_bruteforce.c */
/*
		This file is part of the PVS extension to the Corrfunc package
		Copyright (C) 2016-- Antonio Bibiano (antbbn@gmail.com)
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "pvs_countpairs_rp_pi_bruteforce.h" //function proto-type
#include "utils.h" //all of the utilities
#include "progressbar.h" //for the progressbar

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


void free_results_pvs_rp_pi(results_pvs_countpairs_rp_pi **results)
{
	if(results==NULL)
		return;

	if(*results==NULL)
		return;
	
	results_pvs_countpairs_rp_pi *tmp = *results;

	free(tmp->npairs);
	free(tmp->rupp);
	free(tmp->vpavg);
	free(tmp->spavg);
	free(tmp->vtavg);
	free(tmp->stavg);
	free(tmp->rpavg);
	free(tmp);
	tmp = NULL;
}


results_pvs_countpairs_rp_pi * bruteforce_pvs_countpairs_rp_pi(const int64_t ND1, const DOUBLE *X1, const DOUBLE *Y1, const DOUBLE *Z1,const DOUBLE *VX1, const DOUBLE *VY1, const DOUBLE *VZ1,
                                                    const int64_t ND2, const DOUBLE *X2, const DOUBLE *Y2, const DOUBLE *Z2,const DOUBLE *VX2, const DOUBLE *VY2, const DOUBLE *VZ2,
                                                    const int autocorr,
                                                    const char *binfile,
                                                    const double pimax)
{
  const int npibin = (int) pimax;
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


#ifdef PERIODIC
  const DOUBLE xdiff = (xmax-xmin);
  const DOUBLE ydiff = (ymax-ymin);
  const DOUBLE zdiff = (zmax-zmin);
#endif

  /* DOUBLE logrpmax,logrpmin,dlogrp; */
  const DOUBLE dpi = pimax/(DOUBLE)npibin ;
  const DOUBLE inv_dpi = 1.0/dpi;

  DOUBLE rupp_sqr[nrpbin];
  const int64_t totnbins = (npibin+1)*(nrpbin+1);
  for(int i=0; i < nrpbin;i++) {
    rupp_sqr[i] = rupp[i]*rupp[i];
  }

  const DOUBLE sqr_rpmax=rupp_sqr[nrpbin-1];
  const DOUBLE sqr_rpmin=rupp_sqr[0];

  uint64_t npairs[totnbins];
  DOUBLE vpavg[totnbins];
  DOUBLE spavg[totnbins];
  DOUBLE vtavg[totnbins];
  DOUBLE stavg[totnbins];
#ifdef KAHN_SUM
  DOUBLE vpavg_c[totnbins];
  DOUBLE spavg_c[totnbins];
  DOUBLE vtavg_c[totnbins];
  DOUBLE stavg_c[totnbins];
#endif
#ifdef OUTPUT_RPAVG
  DOUBLE rpavg[totnbins];
#endif	
  for(int ibin=0;ibin<totnbins;ibin++) {
    npairs[ibin]=0;
    vpavg[ibin] = 0.0;
    spavg[ibin] = 0.0;
    vtavg[ibin] = 0.0;
    stavg[ibin] = 0.0;
#ifdef KAHN_SUM
    vpavg_c[ibin] = 0.0;
    spavg_c[ibin] = 0.0;
    vtavg_c[ibin] = 0.0;
    stavg_c[ibin] = 0.0;
#endif
#ifdef OUTPUT_RPAVG		
    rpavg[ibin] = 0.0;
#endif		
  }


  int interrupted=0;
  int64_t numdone=0;
  init_my_progressbar(ND1,&interrupted);


  /*---Loop-over-lattice1--------------------*/
  for(int64_t index1=0;index1<ND1;index1++) {

    my_progressbar(numdone,&interrupted);
    numdone++;

    DOUBLE x1pos = X1[index1];
    DOUBLE y1pos = Y1[index1];
    DOUBLE z1pos = Z1[index1];
    DOUBLE vx1pos = VX1[index1]; //name does not make much sense
    DOUBLE vy1pos = VY1[index1];
    DOUBLE vz1pos = VZ1[index1];

    //#ifdef PERIODIC
    //                    x1pos += off_xwrap;
    //                    y1pos += off_ywrap;
    //                    z1pos += off_zwrap;
    //#endif
    
    // *_countpairs_rp_pi double counts
    for (int64_t index2=0;index2 < ND2; index2++)	{ 
      DOUBLE dx =      X2[index2]-x1pos;
      DOUBLE dy =      Y2[index2]-y1pos;
      DOUBLE dz =      Z2[index2]-z1pos;
      const DOUBLE dvx =      VX2[index2]-vx1pos;
      const DOUBLE dvy =      VY2[index2]-vy1pos;
      const DOUBLE dvz =      VZ2[index2]-vz1pos;

#ifdef PERIODIC
      if (dx > xdiff/2.) dx = dx - xdiff;
      else if (dx < -1*xdiff/2.) dx = dx + xdiff ;

      if (dy > ydiff/2.) dy = dy - ydiff;
      else if (dy < -1*ydiff/2.) dy = dy + ydiff ;

      if (dz > zdiff/2.) dz = dz - zdiff;
      else if (dz < -1*zdiff/2.) dz = dz + zdiff ;
#endif
      const DOUBLE dzabs = FABS(dz);

      const DOUBLE r2 = dx*dx + dy*dy;
      if (r2 >= sqr_rpmax || r2 < sqr_rpmin || dzabs >= pimax) {
        continue;
      }

      const DOUBLE r3 = r2 + dz * dz;
      const DOUBLE vp = (dvx*dx + dvy*dy + dvz*dz)/SQRT(r3);
      const DOUBLE sp = vp*vp;
      const DOUBLE st = (dvx*dvx + dvy*dvy + dvz*dvz)-sp;
      const DOUBLE vt = SQRT(st);

#ifdef OUTPUT_RPAVG								
      const DOUBLE r = SQRT(r2);
#endif								
      int pibin = (int) (dzabs*inv_dpi);
      pibin = pibin > npibin ? npibin:pibin;
      
      DOUBLE tmp1, tmp2;
      for(int kbin=nrpbin-1;kbin>=1;kbin--) {
        if(r2 >= rupp_sqr[kbin-1]) {
          const int ibin = kbin*(npibin+1) + pibin;
          npairs[ibin]++;

#ifndef KAHN_SUM
          vpavg[ibin]+=vp;
#else
          tmp1 = vp - vpavg_c[ibin];
          tmp2 = vpavg[ibin] + tmp1;
          vpavg_c[ibin] = (tmp2 - vpavg[ibin]) - tmp1;
          vpavg[ibin] = tmp2 ;
#endif

#ifndef KAHN_SUM
          spavg[ibin]+=sp;
#else
          tmp1 = sp - spavg_c[ibin];
          tmp2 = spavg[ibin] + tmp1;
          spavg_c[ibin] = (tmp2 - spavg[ibin]) - tmp1;
          spavg[ibin] = tmp2 ;
#endif


#ifndef KAHN_SUM
          vtavg[ibin]+=vt;
#else
          tmp1 = vt - vtavg_c[ibin];
          tmp2 = vtavg[ibin] + tmp1;
          vtavg_c[ibin] = (tmp2 - vtavg[ibin]) - tmp1;
          vtavg[ibin] = tmp2 ;
#endif
          
#ifndef KAHN_SUM
          stavg[ibin]+=st;
#else
          tmp1 = st - stavg_c[ibin];
          tmp2 = stavg[ibin] + tmp1;
          stavg_c[ibin] = (tmp2 - stavg[ibin]) - tmp1;
          stavg[ibin] = tmp2 ;
#endif

#ifdef OUTPUT_RPAVG										
          rpavg[ibin]+=r;
#endif		
          break;
        }
      }
    }//end of index2 loop
  }//end of index1 loop
  finish_myprogressbar(&interrupted);



  for(int i=0;i<totnbins;i++){
    if(npairs[i] > 0) {
      vpavg[i] /= ((DOUBLE) npairs[i] );
      spavg[i] = SQRT(spavg[i]/((DOUBLE) npairs[i] ));
      vtavg[i] /= ((DOUBLE) npairs[i] );
      stavg[i] = SQRT(stavg[i]/((DOUBLE) npairs[i] ));
      //spavg[i] /= ((DOUBLE) npairs[i] );
      //stavg[i] /= ((DOUBLE) npairs[i] );
    }
  }

#ifdef OUTPUT_RPAVG	
  for(int i=0;i<totnbins;i++){
    if(npairs[i] > 0) {
      rpavg[i] /= ((DOUBLE) npairs[i] );
    }
  }
#endif


  //Pack in the results
  results_pvs_countpairs_rp_pi *results = my_malloc(sizeof(*results), 1);
  results->nbin   = nrpbin;
  results->npibin = npibin;
  results->pimax  = pimax;
  results->npairs = my_malloc(sizeof(uint64_t), totnbins);
  results->rupp   = my_malloc(sizeof(DOUBLE)  , nrpbin);
  results->vpavg  = my_malloc(sizeof(DOUBLE)  , totnbins);
  results->spavg  = my_malloc(sizeof(DOUBLE)  , totnbins);
  results->vtavg  = my_malloc(sizeof(DOUBLE)  , totnbins);
  results->stavg  = my_malloc(sizeof(DOUBLE)  , totnbins);
  results->rpavg  = my_malloc(sizeof(DOUBLE)  , totnbins);

  for(int i=0;i<nrpbin;i++) {
    results->rupp[i] = rupp[i];
    for(int j=0;j<npibin;j++) {
      int index = i*(npibin+1) + j;
      assert(index < totnbins && "index must be within range");
      results->npairs[index] = npairs[index];
      results->vpavg[index] = vpavg[index];
      results->spavg[index] = spavg[index];
      results->vtavg[index] = vtavg[index];
      results->stavg[index] = stavg[index];
#ifdef OUTPUT_RPAVG
      results->rpavg[index] = rpavg[index];
#else
      results->rpavg[index] = 0.0;
#endif
    }
  }

  free(rupp);

  return results;
}
