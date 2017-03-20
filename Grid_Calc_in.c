#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <math.h>

#ifdef nompi
#include "mimic_mpi.h"
#else
#include "mpi.h"
#endif

#include "Tools_BandCalc.h"
#include "Inputtools.h"
#include "read_scfout.h"
#include "EigenValue_Problem.h"
#include "GetOrbital.h"


int Grid_Calc()
{
  FILE *fp, *fp1, *fp2;
  char c;
  char fname_FS[256], fname_MP[256], fname_Spin[256];
  char Pdata_s[256];

  int i,j,k,l,m,n,n2, i1,i2, j1,j2, l2;           // loop variable
  int id;
  int Count_data;
  int *S_knum, *E_knum, T_knum, size_Tknum;       // MPI variable (k-point divide)
  int namelen, num_procs, myrank;                 // MPI_variable

  double k1, k2, k3;                              // k-point variable                         
  double **k_xyz, *EIGEN, *EIGEN_Band;                      // Eigen solve array
 
  double d0, d1, d2, d3;

  double TStime, TEtime, Stime, Etime;            // Time variable
  int *S2_knum, *E2_knum, T2_knum;                // MPI variable (k-point divide)

  double line_k, line_k_old, *line_Nk;

  // ### Orbital Data    ###
  double **MulP_Band, ***OrbMulP_Band, ***Orb2MulP_Band;

  int i_vec[20],i_vec2[20];                       // input variable
  char *s_vec[20];                                // input variable
  double r_vec[20];                               // input variable

  double Re11, Re22, Re12, Im12;                  // MulP Calc variable-1
  double Nup[2], Ndw[2], Ntheta[2], Nphi[2];      // MulP Calc variable-2

  // ### MPI_Init ############################################
//  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//  MPI_Get_processor_name(processor_name, &namelen);
//  printf("myrank:%d \n",myrank);

//  dtime(&Stime);

  // ### INPUT_FILE ##########################################

  input_open(fname);
    
  r_vec[0]=-0.2; r_vec[1]=0.1;
  input_doublev("Energy.Range",2,E_Range,r_vec);
  if (E_Range[0]>E_Range[1]){
    d0 = E_Range[0];
    E_Range[0] = E_Range[1];
    E_Range[1] = d0;
  }//if

  r_vec[0]=0.0; r_vec[1]=0.0; r_vec[2]=0.0;
  input_doublev("Search.kCentral",3,k_CNT,r_vec);

  input_int("Eigen.Newton",&switch_Eigen_Newton,0);
  if (switch_Eigen_Newton != 0) switch_Eigen_Newton = 1;
  input_int("Calc.Axis.Grid",&plane_3mesh,1);

  i_vec2[0]=2;  i_vec2[1]=2;
  input_intv("kmesh.Grid.",2,i_vec,i_vec2);
  k1_3mesh = i_vec[0];    k2_3mesh = i_vec[1]; 

  r_vec[0]=0.5; r_vec[1]=0.5;
  input_doublev("kRange.Grid",2,kRange_3mesh,r_vec);
  
  r_vec[0]=1.0; r_vec[1]=1.0; r_vec[2]=1.0;
  input_doublev("MulP.Vec.Scale",3, MulP_VecScale, r_vec);

  input_close();

  // ### Band Total (n2) ###
  k = 1;
  for (i=1; i<=atomnum; i++){ k+= Total_NumOrbs[i]; }
  n = k - 1;    n2 = 2*k + 2;

  // ######################################################
  S_knum = (int*)malloc(sizeof(int)*num_procs);
  E_knum = (int*)malloc(sizeof(int)*num_procs);
   
  // ### (EigenValue Problem) ###
  EIGEN = (double*)malloc(sizeof(double)*n2);
  for (i = 0; i < n2; i++) EIGEN[i] = 0.0;

  // ### (EigenValue Problem) ###
  size_Tknum = (Band_N_perpath[i] + 1);
    
  k_xyz = (double**)malloc(sizeof(double*)*3);
  for(j=0; j<3; j++) k_xyz[j] = (double*)malloc(sizeof(double)*(size_Tknum+1));
  EIGEN_Band = (double*)malloc(sizeof(double)*((size_Tknum+1)*n2));
  for (j = 0; j < ((size_Tknum+1)*n2); j++) EIGEN_Band[j] = 0.0;


}

