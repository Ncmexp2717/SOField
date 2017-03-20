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
#include "GetOrbital.h"
#include "Band_Dispersion.h"
#include "Circular_Search.h"


int main(int argc, char *argv[]) 
{
  FILE *fp;
  int i,j,k, l,n,n2;  
  
  int namelen, num_procs, myrank;                 // MPI_variable
//  char processor_name[MPI_MAX_PROCESSOR_NAME];

  //File*Open
  int i_vec[20],i_vec2[20];                       // input variable
  char *s_vec[20];                                // input variable
  double r_vec[20];                               // input variable
  
  int Calc_Type;

  double TStime, TEtime, Stime, Etime;            // Time variable


// ### MPI_Init ############################################
//  mpi_comm_level1 = MPI_COMM_WORLD;
//  MPI_COMM_WORLD1 = MPI_COMM_WORLD;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//  MPI_Get_processor_name(processor_name, &namelen);
//  printf("myrank:%d \n",myrank);

//  NUMPROCS_MPI_COMM_WORLD = numprocs;
//  MYID_MPI_COMM_WORLD = myid;
//  Num_Procs = numprocs;

  
  //measuring calc time
  dtime(&TStime);

// ### INPUT_FILE ##########################################
  // check argv 
  if (argc==1){
    if (myrank ==0) printf("\nCould not find an input file.\n\n");
    MPI_Finalize();
    exit(0);
  }

  if (myrank ==0) printf("\n");
  sprintf(fname,"%s",argv[1]);
  input_open(fname);

  input_string("Filename.scfout",fname_wf,"default");
  input_string("Filename.outdata",fname_out,"default");

  s_vec[0]="AtomInfo";
  s_vec[1]="CircularSearch";
  s_vec[2]="BandDispersion";
  s_vec[3]="EigenGrid"; 
  s_vec[4]="TriMesh";
  s_vec[5]="MulPonly";
  i_vec[0]=1;  i_vec[1]=2;  i_vec[2]=3;   i_vec[3]=4;   i_vec[4]=5;  i_vec[5]=6;
    
  i = input_string2int("Calc.Type", &Calc_Type, 6, s_vec, i_vec);

  input_close();

  if (i < 0)  exit(0);

  // ### Get Calculation Data ##############################
  if((fp = fopen(fname_wf,"r")) != NULL){
    if (myrank ==0) printf("\nInput filename is \"%s\"  \n\n", fname_wf);
    fclose(fp);
  }else{
    if (myrank ==0) printf("Cannot open *.scfout File. \"%s\" is not found.\n" ,fname_wf);
    exit(0);
  }

  // ### wave functions　###
  if (myrank == 0){
    read_scfout(fname_wf, 1);
  }else {
    read_scfout(fname_wf, 0);
  }

  // #####################################################################
  // #####################################################################
  // #####################################################################

  // ### Total Num Orbs  ###
  TNO_MAX = 0;
  for(i=1;i<=atomnum;i++){
    if(TNO_MAX < Total_NumOrbs[i]) TNO_MAX = Total_NumOrbs[i];
  }
  // ### Classify Orbs   ###
  An2Spe = (char**)malloc(sizeof(char*)*(atomnum+1));
  for (i=0; i<=atomnum; i++){
    An2Spe[i] = (char*)malloc(sizeof(char)*(asize10));
  }
  ClaOrb = (int**)malloc(sizeof(int*)*(atomnum+1));
  for (i=0; i<=atomnum; i++){
    ClaOrb[i] = (int*)malloc(sizeof(int)*(TNO_MAX+1));
    for (j=0; j<=TNO_MAX; j++) ClaOrb[i][j]=0;
  } 
  MPI_Barrier(MPI_COMM_WORLD);

  if (myrank == 0){ 
    Classify_OrbNum(ClaOrb, An2Spe, 1);
  }else{
    Classify_OrbNum(ClaOrb, An2Spe, 0);
  }

  ClaOrb_MAX[1] = 0;
  for(i=1;i<=atomnum;i++){
    for (j=0; j<=TNO_MAX; j++){
      if(ClaOrb_MAX[1] < ClaOrb[i][j]) ClaOrb_MAX[1] = ClaOrb[i][j];
    }
  } ClaOrb_MAX[0] = 0;
  OrbSym[0][0] = 's';  OrbSym[0][1] = '\0';
  OrbName = (char**)malloc(sizeof(char*)*(ClaOrb_MAX[1]+1));
  for (i=0; i<=ClaOrb_MAX[1]; i++){
    OrbName[i] = (char*)malloc(sizeof(char)*3);
    if (i == 0){// 0
      OrbName[i][0]='s';  OrbName[i][1]=(char)(1+48);      OrbName[i][2]='\0';
    }else if(i > 0 && i < 4){// 1,2,3
      if (ClaOrb_MAX[0] < 1){ ClaOrb_MAX[0] = 1; OrbSym[1][0] = 'p';   OrbSym[1][1] = '\0';}
      OrbName[i][0]='p';  OrbName[i][1]=(char)((i+0)+48);  OrbName[i][2]='\0';
    }else if(i > 3 && i < 9){// 4,5,6,7,8
      if (ClaOrb_MAX[0] < 2){ ClaOrb_MAX[0] = 2; OrbSym[2][0] = 'd';   OrbSym[2][1] = '\0';}
      OrbName[i][0]='d';  OrbName[i][1]=(char)((i-3)+48);  OrbName[i][2]='\0';
    }else if(i >8 && i <16){// 9,10,11,12,13,14,15
      if (ClaOrb_MAX[0] < 3){ ClaOrb_MAX[0] = 3; OrbSym[3][0] = 'f';   OrbSym[3][1] = '\0';}
      OrbName[i][0]='f';  OrbName[i][1]=(char)((i-8)+48);  OrbName[i][2]='\0';
    }else{//16
      OrbName[i][0]=(char)((i-15)/10+48);
      OrbName[i][1]=(char)((i-15)%10+48);  OrbName[i][2]='\0';
    } 
//    if (myrank == 0) printf("OrbName:%s\n",OrbName[i]);
  }//i
//  if (myrank == 0) for (i=1; i<=atomnum; i++){  printf("%4d %s\n", i, An2Spe[i]);  }

  // ### Band Total (n2) ###
  k = 1;
  for (i=1; i<=atomnum; i++){ k+= Total_NumOrbs[i]; }
  n = k - 1;    n2 = 2*k + 2;

  if (myrank == 0){
    printf("########### ORBITAL DATA ##################\n");
//    for(i=1;i<=atomnum;i++) printf("%4d:%4d\n", i, Total_NumOrbs[i]);
//    printf("  MAX:%4d\n",TNO_MAX);
    printf("ClaOrb_MAX[0]:%4d\n",ClaOrb_MAX[0]);
    printf("ClaOrb_MAX[1]:%4d\n",ClaOrb_MAX[1]);
    printf("Total Band (2*n):%4d\n",n2-4);
//    printf("Central (%10.6lf %10.6lf %10.6lf)\n",k_CNT1[0],k_CNT1[1],k_CNT1[2]);
    printf("###########################################\n");
  }

  // #####################################################################
  // #####################################################################
  // #####################################################################

  // ### Calclation ###

  if (myrank ==0) printf("\nStart \"%s\" Calculation (%d).  \n", s_vec[Calc_Type-1], Calc_Type); 

  if (Calc_Type == 1){
  } else if (Calc_Type == 2){
    Circular_Search();
  } else if (Calc_Type == 3){
    Band_Dispersion();
  }

  // ### MALLOC FREE #######################################

  for (i=0; i<=ClaOrb_MAX[1]; i++){
    free(OrbName[i]);
  } free(OrbName);
  // ### (MulP Calculation)   ###
  for (i=0; i<=atomnum; i++){
    free(An2Spe[i]);
  } free(An2Spe);
  for (i=0; i<=atomnum; i++){
    free(ClaOrb[i]);
  } free(ClaOrb);

  free_scfout();

  if (myrank ==0) printf("############ CALC TIME ####################\n");
  dtime(&TEtime);
  if (myrank ==0) printf("  Total Calculation Time:%10.6lf (s)\n",TEtime-TStime);
  if (myrank ==0) printf("###########################################\n");

  // ### MPI_Finalize ########################################
  MPI_Finalize();

  return 0;
}

