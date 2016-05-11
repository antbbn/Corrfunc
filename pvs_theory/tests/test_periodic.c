/* File: test_periodic.c */
/*
		This file is part of the PVS extension to the Corrfunc package
		Copyright (C) 2016-- Antonio Bibiano (antbbn@gmail.com)
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://github.com/manodeep/Corrfunc/
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>

#ifndef MAXLEN
#define MAXLEN 500
#endif

#ifndef PERIODIC
#define PERIODIC
#endif

#ifndef DOUBLE_PREC
#define DOUBLE_PREC
#endif

#ifndef OUTPUT_RPAVG
#define OUTPUT_RPAVG
#endif

// Kahn summation is needed to make the tests more accurate.
#ifndef KAHN_SUM
#define KAHN_SUM
#endif

#ifndef SILENT
#define SILENT
#endif

#include "function_precision.h"
#include "io.h"
#include "defs.h"
#include "utils.h"

//Including the C files directly
//#include "gridlink.c"
#include "gridlink_pvs.c"
#include "io.c"
#include "ftread.c"
#include "../pvs_of_r/pvs_countpairs.c"
#include "../pvs_rp_pi/pvs_countpairs_rp_pi.c"
//#include "../wp/countpairs_wp.c"
//#include "../vpf/countspheres.c"
//#include "../xi/countpairs_xi.c"

char tmpoutputfile[]="./test_periodic_output.txt";

int test_periodic_DDpvs(const char *correct_outputfile);
int test_periodic_DDrppi_pvs(const char *correct_outputfile);
//int test_wp(const char *correct_outputfile);
//int test_vpf(const char *correct_outputfile);
//int test_xi(const char *correct_outputfile);

void read_data_and_set_globals(const char *firstfilename, const char *firstformat,const char *secondfilename,const char *secondformat);

//Global variables
int ND1;
DOUBLE *X1=NULL,*Y1=NULL,*Z1=NULL;
DOUBLE *VX1=NULL,*VY1=NULL,*VZ1=NULL;

int ND2;
DOUBLE *X2=NULL,*Y2=NULL,*Z2=NULL;
DOUBLE *VX2=NULL,*VY2=NULL,*VZ2=NULL;

char binfile[]="../../xi_theory/tests/bins";
DOUBLE pimax=40.0;
double boxsize=420.0;
#ifdef USE_OMP
const int nthreads=4;
#endif

char current_file1[MAXLEN],current_file2[MAXLEN];

//end of global variables


int test_periodic_DDpvs(const char *correct_outputfile)
{
    int autocorr = (X1==X2) ? 1:0;

    //Do the straight-up DD counts
    results_pvs_countpairs *results = countpairs_pvs(ND1,X1,Y1,Z1,VX1,VY1,VZ1,
                                             ND2,X2,Y2,Z2,VX2,VY2,VZ2,
#ifdef USE_OMP
                                             nthreads,
#endif
                                             autocorr,
                                             binfile);

    DOUBLE rlow=results->rupp[0];
    FILE *fp=NULL;
    fp=my_fopen(tmpoutputfile,"w");
    for(int i=1;i<results->nbin;i++) {
        fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf %20.8lf %20.8lf %20.8lf %20.8lf %20.8lf \n",results->npairs[i],results->vpavg[i],results->spavg[i],results->vtavg[i],results->stavg[i],results->rpavg[i],rlow,results->rupp[i]);
        rlow=results->rupp[i];
    }
    fclose(fp);

    char execstring[MAXLEN];
    my_snprintf(execstring,MAXLEN,"diff -q %s %s",correct_outputfile,tmpoutputfile);
    int ret=system(execstring);

    free_results_pvs(&results);
    return ret;
}



int test_periodic_DDrppi_pvs(const char *correct_outputfile)
{
    int autocorr = (X1==X2) ? 1:0;

    results_pvs_countpairs_rp_pi *results = pvs_countpairs_rp_pi(ND1,X1,Y1,Z1,VX1,VY1,VZ1,
                                                         ND2,X2,Y2,Z2,VX2,VY2,VZ2,
#ifdef USE_OMP
                                                         nthreads,
#endif
                                                         autocorr,
                                                         binfile,
                                                         pimax);

    const int npibin = results->npibin;
    const DOUBLE dpi = pimax/(DOUBLE)results->npibin ;
    FILE *fp=my_fopen(tmpoutputfile,"w");
    for(int i=1;i<results->nbin;i++) {
        const double logrp = LOG10(results->rupp[i]);
        for(int j=0;j<npibin;j++) {
            int index = i*(npibin+1) + j;
            fprintf(fp,"%10"PRIu64" %20.8lf %20.8lf %20.8lf %20.8lf %20.8lf %20.8lf  %20.8lf \n",results->npairs[index],results->vpavg[index],results->spavg[index],results->vtavg[index],results->stavg[index],results->rpavg[index],logrp,(j+1)*dpi);
        }
    }
    fclose(fp);
    char execstring[MAXLEN];
    my_snprintf(execstring,MAXLEN,"diff -q %s %s",correct_outputfile,tmpoutputfile);
    int ret=system(execstring);

    //free the result structure
    free_results_pvs_rp_pi(&results);
    return ret;
}
/*
int test_wp(const char *correct_outputfile)
{
    results_countpairs_wp *results = countpairs_wp(ND1,X1,Y1,Z1,
                                                   boxsize,
#ifdef USE_OMP
                                                   nthreads,
#endif
                                                   binfile,
                                                   pimax);
    DOUBLE rlow=results->rupp[0];
    FILE *fp=my_fopen(tmpoutputfile,"w");
    for(int i=1;i<results->nbin;++i) {
        fprintf(fp,"%e\t%e\t%e\t%e\t%12"PRIu64" \n",results->wp[i],results->rpavg[i],rlow,results->rupp[i],results->npairs[i]);
        rlow=results->rupp[i];
    }
    fclose(fp);
    char execstring[MAXLEN];
    my_snprintf(execstring,MAXLEN,"diff -q %s %s",correct_outputfile,tmpoutputfile);
    int ret=system(execstring);

    //free the result structure
    free_results_wp(&results);
    return ret;
}

int test_vpf(const char *correct_outputfile)
{
    const double rmax = 10.0;
    const int nbin = 10;
    const int nc = 10000;
    const int num_pN=6;
    const unsigned long seed=-1234;
    results_countspheres *results = countspheres(ND1, X1, Y1, Z1,
                                                 rmax, nbin, nc,
                                                 num_pN,
                                                 seed);

    FILE *fp=my_fopen(tmpoutputfile,"w");
    const DOUBLE rstep = rmax/(DOUBLE)nbin ;
    for(int ibin=0;ibin<results->nbin;ibin++) {
        const double r=(ibin+1)*rstep;
        fprintf(fp,"%"DOUBLE_FORMAT" ", r);
        for(int i=0;i<num_pN;i++) {
            fprintf(fp," %10.4e", (results->pN)[ibin][i]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    char execstring[MAXLEN];
    my_snprintf(execstring,MAXLEN,"diff -q %s %s",correct_outputfile,tmpoutputfile);
    int ret=system(execstring);

    //free the result structure
    free_results_countspheres(&results);
    return ret;
}

int test_xi(const char *correct_outputfile)
{

    results_countpairs_xi *results = countpairs_xi(ND1,X1,Y1,Z1,
                                                   boxsize,
#ifdef USE_OMP
                                                   nthreads,
#endif
                                                   binfile);
    DOUBLE rlow=results->rupp[0];
    FILE *fp=my_fopen(tmpoutputfile,"w");
    for(int i=1;i<results->nbin;++i) {
        fprintf(fp,"%e\t%e\t%e\t%e\t%12"PRIu64" \n",results->xi[i],results->rpavg[i],rlow,results->rupp[i],results->npairs[i]);
        rlow=results->rupp[i];
    }
    fclose(fp);
    char execstring[MAXLEN];
    my_snprintf(execstring,MAXLEN,"diff -q %s %s",correct_outputfile,tmpoutputfile);
    int ret=system(execstring);

    //free the result structure
    free_results_xi(&results);
    return ret;
}
*/

void read_data_and_set_globals(const char *firstfilename, const char *firstformat,const char *secondfilename,const char *secondformat)
{
    int free_X2=0;
    if(X2 != NULL && X2 != X1) {
        free_X2=1;
    }


    //Check to see if data has to be read for X1/Y1/Z1
    if (strncmp(current_file1,firstfilename,strlen(current_file1)) != 0) {
        //replace the data-set
        if(X1 != NULL) {
            free(X1);free(Y1);free(Z1);
            free(VX1);free(VY1);free(VZ1);
        }

        //Since X2 was pointing to X1, need to signal that the memory location is no longer valid
        if(free_X2 == 0) {
            X2 = NULL;
            Y2 = NULL;
            Z2 = NULL;
            VX2 = NULL;
            VY2 = NULL;
            VZ2 = NULL;
        }
        ND1 = read_positions(firstfilename,firstformat, sizeof(DOUBLE), 6, &X1, &Y1, &Z1, &VX1, &VY1, &VZ1);
        strncpy(current_file1,firstfilename,MAXLEN);
    }

    //first check if only one unique file is asked for
    if(strncmp(firstfilename,secondfilename,strlen(firstfilename))==0) {
        //But X2 might have read-in a different file->avoid a memory-leak
        if(free_X2 == 1) {
            free(X2);free(Y2);free(Z2);
            free(VX2);free(VY2);free(VZ2);
            free_X2 = 0;//not essential since the code returns after this section
        }
        X2=X1;
        Y2=Y1;
        Z2=Z1;
        VX2=VX1;
        VY2=VY1;
        VZ2=VZ1;
        ND2=ND1;
        strncpy(current_file2,secondfilename,MAXLEN);
        return;
    }


    //Check to see if data has to be read for X2/Y2/Z2
    if (strncmp(current_file2,secondfilename,strlen(current_file2)) != 0 || X2 == NULL) {
        //replace the data-set
        if(free_X2 == 1) {
            free(X2);free(Y2);free(Z2);
            free(VX2);free(VY2);free(VZ2);
        }
        ND2 = read_positions(secondfilename,secondformat, sizeof(DOUBLE), 6, &X2, &Y2, &Z2, &VX2, &VY2, &VZ2);
        strncpy(current_file2,secondfilename,MAXLEN);
    }
}


int main(int argc, char **argv)
{
    struct timeval tstart,t0,t1;
    char file[]="../../xi_theory/tests/data/gals_Mr19.ff";
    char fileformat[]="f";

    gettimeofday(&tstart,NULL);

    //set the globals
    ND1 = read_positions(file,fileformat, sizeof(DOUBLE), 6, &X1, &Y1, &Z1, &VX1, &VY1, &VZ1);
    ND2 = ND1;
    X2 = X1;
    Y2 = Y1;
    Z2 = Z1;
    VX2 = VX1;
    VY2 = VY1;
    VZ2 = VZ1;

    strncpy(current_file1,file,MAXLEN);
    strncpy(current_file2,file,MAXLEN);

    int failed=0;
    int status;

    //const char alltests_names[][MAXLEN] = {"Mr19 DD (periodic)","Mr19 DDrppi (periodic)","Mr19 wp (periodic)","Mr19 vpf (periodic)","Mr19 xi periodic)",
    const char alltests_names[][MAXLEN] = {"Mr19 DDrppi_pvs (periodic)","Mr19 DDpvs (periodic)"};
    const int ntests = sizeof(alltests_names)/(sizeof(char)*MAXLEN);
    const int function_pointer_index[] = {0,1};//0->DDrppi_pvs, 1-> DDpvs

    const char correct_outoutfiles[][MAXLEN] = {"Mr19_DDrppi_pvs_periodic","Mr19_DDpvs_periodic"};
    const char firstfilename[][MAXLEN] = {"../../xi_theory/tests/data/gals_Mr19.ff","../../xi_theory/tests/data/gals_Mr19.ff"};
    const char firstfiletype[][MAXLEN] = {"f","f"};
    const char secondfilename[][MAXLEN] = {"../../xi_theory/tests/data/gals_Mr19.ff","../../xi_theory/tests/data/gals_Mr19.ff"};
    const char secondfiletype[][MAXLEN] = {"f","f"};
    const DOUBLE allpimax[]             = {40.0,40.0};

    int (*allfunctions[]) (const char *) = {test_periodic_DDrppi_pvs, test_periodic_DDpvs};
    const int numfunctions=2;//5 functions total

    int total_tests=0,skipped=0;

    if(argc==1) {
        //nothing was passed at the command-line -> run all tests
        for(int i=0;i<ntests;i++) {
            int function_index = function_pointer_index[i];
            assert(function_index >= 0 && function_index < numfunctions && "Function index is within range");
            const char *testname = alltests_names[i];
            int skip_test=test_all_files_present(2,firstfilename[i],secondfilename[i]);
            if(skip_test != 0) {
                fprintf(stderr,ANSI_COLOR_YELLOW "SKIPPED: " ANSI_COLOR_MAGENTA "%s"  ANSI_COLOR_RESET ". File(s) not found\n", testname);
                skipped++;
                continue;
            }
            read_data_and_set_globals(firstfilename[i],firstfiletype[i],secondfilename[i],secondfiletype[i]);
            pimax=allpimax[i];
            gettimeofday(&t0,NULL);
            status = (*allfunctions[function_index])(correct_outoutfiles[i]);
            gettimeofday(&t1,NULL);
            double pair_time = ADD_DIFF_TIME(t0,t1);
            total_tests++;
            if(status==EXIT_SUCCESS) {
                fprintf(stderr,ANSI_COLOR_GREEN "PASSED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_GREEN ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                char execstring[MAXLEN];
                my_snprintf(execstring,MAXLEN,"rm -f %s",tmpoutputfile);
                run_system_call(execstring);
            } else {
                fprintf(stderr,ANSI_COLOR_RED "FAILED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_RED ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                failed++;
                char execstring[MAXLEN];
                my_snprintf(execstring,MAXLEN,"mv %s %s.%d",tmpoutputfile,tmpoutputfile,i);
                run_system_call(execstring);

            }
        }
    } else {
        //run specific tests
        for(int i=1;i<argc;i++){
            int this_test_num = atoi(argv[i]);
            if(this_test_num >= 0 && this_test_num < ntests) {
                int function_index = function_pointer_index[this_test_num];
                assert(function_index >= 0 && function_index < numfunctions && "Function index is within range");
                const char *testname = alltests_names[this_test_num];
                int skip_test=test_all_files_present(2,firstfilename[this_test_num],secondfilename[this_test_num]);
                if(skip_test != 0) {
                    fprintf(stderr,ANSI_COLOR_YELLOW "SKIPPED: " ANSI_COLOR_MAGENTA "%s"  ANSI_COLOR_RESET ". File(s) not found\n", testname);
                    skipped++;
                    continue;
                }
                total_tests++;
                read_data_and_set_globals(firstfilename[this_test_num],firstfiletype[this_test_num],secondfilename[this_test_num],secondfiletype[this_test_num]);
                pimax=allpimax[this_test_num];
                gettimeofday(&t0,NULL);
                status = (*allfunctions[function_index])(correct_outoutfiles[this_test_num]);
                gettimeofday(&t1,NULL);
                double pair_time = ADD_DIFF_TIME(t0,t1);
                if(status==EXIT_SUCCESS) {
                    fprintf(stderr,ANSI_COLOR_GREEN "PASSED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_GREEN ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                    char execstring[MAXLEN];
                    my_snprintf(execstring,MAXLEN,"rm -f %s",tmpoutputfile);
                    run_system_call(execstring);
                } else {
                    fprintf(stderr,ANSI_COLOR_RED "FAILED: " ANSI_COLOR_MAGENTA "%s" ANSI_COLOR_RED ". Time taken = %8.2lf seconds " ANSI_COLOR_RESET "\n", testname,pair_time);
                    failed++;
                    char execstring[MAXLEN];
                    my_snprintf(execstring,MAXLEN,"mv %s %s.%d",tmpoutputfile,tmpoutputfile,this_test_num);
                    run_system_call(execstring);

                }
            }
        }
    }


    gettimeofday(&t1,NULL);
    double total_time = ADD_DIFF_TIME(tstart,t1);
    if(failed > 0) {
        fprintf(stderr,ANSI_COLOR_RED "FAILED %d out of %d tests. Total time = %8.2lf seconds " ANSI_COLOR_RESET "\n", failed, total_tests, total_time);
    } else {
        fprintf(stderr,ANSI_COLOR_GREEN "PASSED: ALL %d tests. Total time = %8.2lf seconds " ANSI_COLOR_RESET "\n", total_tests, total_time);
    }
    if(skipped > 0) {
        fprintf(stderr,ANSI_COLOR_YELLOW "SKIPPED: %d tests" ANSI_COLOR_RESET "\n", skipped);
        fprintf(stderr,ANSI_COLOR_MAGENTA "Tests are skipped on the PyPI installed code-base. Please use the git repo if you want to run the entire suite of tests"ANSI_COLOR_RESET"\n\n");
    }

    if(X2 != X1) {
        free(X2);free(Y2);free(Z2);
        free(VX2);free(VY2);free(VZ2);
    }
    free(X1);free(Y1);free(Z1);
    free(VX1);free(VY1);free(VZ1);
    return EXIT_SUCCESS;
}