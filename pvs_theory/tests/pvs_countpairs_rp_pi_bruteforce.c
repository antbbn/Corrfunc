/* File: pvs_countpairs_rp_pi.c */
/*
		This file is a part of the Corrfunc package
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "pvs_countpairs_rp_pi.h" //function proto-type
#include "gridlink_pvs.h"//function proto-type for gridlink
#include "cellarray_pvs.h" //definition of struct cellarray_pvs
#include "utils.h" //all of the utilities
#include "progressbar.h" //for the progressbar

#if defined(USE_AVX) && defined(__AVX__)
#include "avx_calls.h"
#endif


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


results_pvs_countpairs_rp_pi * pvs_countpairs_rp_pi(const int64_t ND1, const DOUBLE *X1, const DOUBLE *Y1, const DOUBLE *Z1,const DOUBLE *VX1, const DOUBLE *VY1, const DOUBLE *VZ1,
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
#ifdef OUTPUT_RPAVG
  DOUBLE rpavg[totnbins];
#endif	
  for(int ibin=0;ibin<totnbins;ibin++) {
    npairs[ibin]=0;
    vpavg[ibin] = 0.0;
    spavg[ibin] = 0.0;
    vtavg[ibin] = 0.0;
    stavg[ibin] = 0.0;
#ifdef OUTPUT_RPAVG		
    rpavg[ibin] = 0.0;
#endif		
  }



  int interrupted=0;
  int64_t numdone=0;
  init_my_progressbar(totncells,&interrupted);


  /*---Loop-over-lattice1--------------------*/
  for(int64_t index1=0;index1<ND1;index1++) {

    my_progressbar(numdone,&interrupted);
    numdone++;

    DOUBLE x1pos = x1[index1];
    DOUBLE y1pos = y1[index1];
    DOUBLE z1pos = z1[index1];
    DOUBLE vx1pos = vx1[index1]; //name does not make much sense
    DOUBLE vy1pos = vy1[index1];
    DOUBLE vz1pos = vz1[index1];

    //#ifdef PERIODIC
    //                    x1pos += off_xwrap;
    //                    y1pos += off_ywrap;
    //                    z1pos += off_zwrap;
    //#endif

    for (int64_t index2=index1+1;index2 < ND2; index2++)	{
      const DOUBLE dx =      localx2[index2]-x1pos;
      const DOUBLE dy =      localy2[index2]-y1pos;
      const DOUBLE dz =      localz2[index2]-z1pos;
      const DOUBLE dzabs = FABS(dz);
      const DOUBLE dvx =      localvx2[index2]-vx1pos;
      const DOUBLE dvy =      localvy2[index2]-vy1pos;
      const DOUBLE dvz =      localvz2[index2]-vz1pos;

#ifdef PERIODIC
      if (dx > xdiff/2.) dx = xdiff - dx
      else if (dx < -1*xdiff/2.) dx = xdiff + dx

      if (dy > ydiff/2.) dy = ydiff - dy
      else if (dy < -1*ydiff/2.) dy = ydiff + dy

      if (dz > zdiff/2.) dz = zdiff - dz
      else if (dz < -1*zdiff/2.) dz = zdiff + dz
#endif

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
      for(int kbin=nrpbin-1;kbin>=1;kbin--) {
        if(r2 >= rupp_sqr[kbin-1]) {
          const int ibin = kbin*(npibin+1) + pibin;
          npairs[ibin]++;
          vpavg[ibin]+=vp;
          spavg[ibin]+=sp;
          vtavg[ibin]+=vt;
          stavg[ibin]+=st;
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
  for(int64_t i=0;i<totncells;i++) {
    free(lattice1[i].x);
    free(lattice1[i].y);
    free(lattice1[i].z);
    free(lattice1[i].vx);
    free(lattice1[i].vy);
    free(lattice1[i].vz);
    if(autocorr==0) {
      free(lattice2[i].x);
      free(lattice2[i].y);
      free(lattice2[i].z);
      free(lattice2[i].vx);
      free(lattice2[i].vy);
      free(lattice2[i].vz);
    }
  }

  free(lattice1);
  if(autocorr==0) {
    free(lattice2);
  }
	
  return results;
}
