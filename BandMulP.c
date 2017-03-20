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


int main(int argc, char *argv[])
{
  FILE *fp, *fp1, *fp2;
  char c;
  char fname_BandAll[256], fname_Band[256], fname_MP[256];
  char Pdata_s[256];

  int i,j,k,l,m,n,n2, i1,i2, j1,j2, l2;           // loop variable
  int id;
  int Count_data;
  int l_max, l_min, l_cal;                        // band index
  int *S_knum, *E_knum, T_knum, size_Tknum;       // MPI variable (k-point divide)
  int namelen, num_procs, myrank;                 // MPI_variable

  double k1, k2, k3;                              // k-point variable                         
  double **k_xyz, *EIGEN_Band;                      // Eigen solve array

  double d0, d1, d2, d3;

  double TStime, TEtime, Stime, Etime;            // Time variable
  int *S2_knum, *E2_knum, T2_knum;                // MPI variable (k-point divide)

  double line_k, line_k_old, *line_Nk;

  // ### Orbital Data    ###
  int TNO_MAX ,ClaOrb_MAX[2] ,**ClaOrb;              // Orbital  
  //  char 
  char OrbSym[5][2], **OrbName, **An2Spe;
  double **MulP_Band, ***OrbMulP_Band, ***Orb2MulP_Band;

  int fo_inp = 0;
  int fo_wf  = 0;
  int i_vec[20],i_vec2[20];                       // input variable
  char *s_vec[20];                                // input variable
  double r_vec[20];                               // input variable

  double Re11, Re22, Re12, Im12;                  // MulP Calc variable-1
  double Nup[2], Ndw[2], Ntheta[2], Nphi[2];      // MulP Calc variable-2


  // ### MPI_Init ############################################
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //  MPI_Get_processor_name(processor_name, &namelen);
  //  printf("myrank:%d \n",myrank);
  dtime(&TStime);

  // ### INPUT_FILE ##########################################
  if (myrank ==0) printf("\n");
  sprintf(fname,"%s",argv[1]);
  if((fp = fopen(fname,"r")) != NULL){
    fo_inp = 1;
    if (myrank ==0) printf("open \"%s\" file... \n" ,fname);
    fclose(fp);
  }else{
    fo_inp = 0;
    if (myrank ==0) printf("Cannot open \"%s\" file.\n" ,fname);
  }

  if (fo_inp == 1){
    input_open(fname);
    input_string("Filename.scfout",fname_wf,"default");
    input_string("Filename.outdata",fname_out,"default");

    r_vec[0]=-0.2; r_vec[1]=0.1;
    input_doublev("Energy.Range",2,E_Range,r_vec);
    if (E_Range[0]>E_Range[1]){
      d0 = E_Range[0];
      E_Range[0] = E_Range[1];
      E_Range[1] = d0;
    }//if

    //BAND
    input_int("Band.Nkpath",&Band_Nkpath,0);
    if (Band_Nkpath>0) {

      Band_N_perpath=(int*)malloc(sizeof(int)*(Band_Nkpath+1));
      for (i=0; i<(Band_Nkpath+1); i++) Band_N_perpath[i] = 0;

      Band_kpath = (double***)malloc(sizeof(double**)*(Band_Nkpath+1));
      for (i=0; i<(Band_Nkpath+1); i++){
	Band_kpath[i] = (double**)malloc(sizeof(double*)*3);
	for (j=0; j<3; j++){
	  Band_kpath[i][j] = (double*)malloc(sizeof(double)*4);
	  for (k=0; k<4; k++) Band_kpath[i][j][k] = 0.0;
	}//for(i)
      }//for(j)
      Band_kname = (char***)malloc(sizeof(char**)*(Band_Nkpath+1));
      for (i=0; i<(Band_Nkpath+1); i++){
	Band_kname[i] = (char**)malloc(sizeof(char*)*3);
	for (j=0; j<3; j++){
	  Band_kname[i][j] = (char*)malloc(sizeof(char)*16);
	}//for(i)
      }//for(j)
      if (fp=input_find("<Band.kpath") ) {
	for (i=1; i<=Band_Nkpath; i++) {
	  fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %s %s",
	      &Band_N_perpath[i]  ,
	      &Band_kpath[i][1][1], &Band_kpath[i][1][2], &Band_kpath[i][1][3],
	      &Band_kpath[i][2][1], &Band_kpath[i][2][2], &Band_kpath[i][2][3],
	      Band_kname[i][1], Band_kname[i][2]);
	}//for(Band_Nkpath)
	if ( ! input_last("Band.kpath>") ) { // format error
	  if (myrank ==0) printf("Format error near Band.kpath>\n");
	  fo_inp = 0;
	}
      }//else
      else { // format error
	if (myrank ==0) printf("<Band.kpath is necessary.\n");
	fo_inp = 0;
      }//else
    }//if(Nkpath>0)
    input_close();
  }//fo_inp

  // ### READ_SCFout_FILE ####################################
  if (fo_inp == 1){
    if((fp = fopen(fname_wf,"r")) != NULL){
      fo_wf = 1;
      if (myrank ==0) printf("\nInput filename is \"%s\"  \n\n", fname_wf);
      fclose(fp);
    }else{
      fo_wf = 0;
      if (myrank ==0) printf("Cannot open *.scfout File. \"%s\" is not found.\n" ,fname_wf);
    }
  }//if(fo_inp)

  if ((fo_inp == 1) && (fo_wf == 1)){
    // ### Get Calculation Data ##############################
    // ### wave functions　###
    if (myrank == 0){
      read_scfout(fname_wf, 1);
    }else {
      read_scfout(fname_wf, 0);
    }

    line_Nk=(double*)malloc(sizeof(double)*(Band_Nkpath+1));
    line_Nk[0] = 0;
    for (i=1; i<(Band_Nkpath+1); i++){
      line_Nk[i] = kgrid_dist(Band_kpath[i][1][1], Band_kpath[i][2][1], Band_kpath[i][1][2], Band_kpath[i][2][2], Band_kpath[i][1][3], Band_kpath[i][2][3]) +line_Nk[i-1]; 
      if (myrank ==0) printf("line_Nk[%d]:%10.6lf\n",i,line_Nk[i]);
      if (myrank ==0) printf("%10.6lf-%10.6lf,%10.6lf-%10.6lf,%10.6lf-%10.6lf\n",Band_kpath[i][1][1], Band_kpath[i][2][1], Band_kpath[i][1][2], Band_kpath[i][2][2], Band_kpath[i][1][3], Band_kpath[i][2][3]);
    }
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
    }
    ClaOrb_MAX[0] = 0;
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
      //      if (myrank == 0) printf("OrbName:%s\n",OrbName[i]);
    }//i
    //    if (myrank == 0) for (i=1; i<=atomnum; i++){  printf("%4d %s\n", i, An2Spe[i]);  }
    // ### Band Total (n2) ###
    k = 1;
    for (i=1; i<=atomnum; i++){ k+= Total_NumOrbs[i]; }
    n = k - 1;    n2 = 2*k + 2;

    if (myrank == 0){
      printf("########### ORBITAL DATA ##################\n");
      //      for(i=1;i<=atomnum;i++) printf("%4d:%4d\n", i, Total_NumOrbs[i]);
      //      printf("  MAX:%4d\n",TNO_MAX);
      printf("ClaOrb_MAX[0]:%4d\n",ClaOrb_MAX[0]);
      printf("ClaOrb_MAX[1]:%4d\n",ClaOrb_MAX[1]);
      printf("Total Band (2*n):%4d\n",n2-4);
      //      printf("Central (%10.6lf %10.6lf %10.6lf)\n",k_CNT1[0],k_CNT1[1],k_CNT1[2]);
      printf("###########################################\n");
    }

    // ######################################################
    S_knum = (int*)malloc(sizeof(int)*num_procs);
    E_knum = (int*)malloc(sizeof(int)*num_procs);

    // ### (EigenValue Problem) ###
    EIGEN = (double*)malloc(sizeof(double)*n2);
    for (i = 0; i < n2; i++) EIGEN[i] = 0.0;
    // ### (MulP Calculation)   ###
    Data_MulP = (double****)malloc(sizeof(double***)*4);
    for (i=0; i<4; i++){
      Data_MulP[i] = (double***)malloc(sizeof(double**)*n2);
      for (l=0; l<n2; l++){
	Data_MulP[i][l] = (double**)malloc(sizeof(double*)*(atomnum+1));
	for (j=0; j<=atomnum; j++){
	  Data_MulP[i][l][j] = (double*)malloc(sizeof(double)*(TNO_MAX+1));
	  for (k=0; k<=TNO_MAX; k++) Data_MulP[i][l][j][k] = 0.0;
	}
      }
    }
    if (Band_Nkpath>0) {
      // ################################################## 
      if (myrank ==0) {
	printf("Band.Nkpath: %d\n",Band_Nkpath);
	for (i=1; i<=Band_Nkpath; i++) {
	  printf("%d (%lf %lf %lf) >>> (%lf %lf %lf) %s %s\n",
	      Band_N_perpath[i] ,
	      Band_kpath[i][1][1], Band_kpath[i][1][2],Band_kpath[i][1][3],
	      Band_kpath[i][2][1], Band_kpath[i][2][2],Band_kpath[i][2][3],
	      Band_kname[i][1],Band_kname[i][2]);
	}
      }
      // ################################################## 
      Count_data = 0;
      if (myrank ==0) {
	printf("### BAND_OUTPUT ###########################\n");
	printf("se st d l\n");
	printf("se para\n");
	printf("se xtics (");
	for (i=1; i<(Band_Nkpath+1); i++){
	  printf("\"%s\" %lf, ", Band_kname[i][1], line_Nk[i-1]);
	}printf("\"%s\"  %lf)\n", Band_kname[Band_Nkpath][2], line_Nk[Band_Nkpath]);
	printf("se xra [%10.6lf:%10.6lf]\n",line_Nk[0],line_Nk[Band_Nkpath]);
	printf("se yra [%10.6lf:%10.6lf]\n",E_Range[0],E_Range[1]);
	printf("se tra [%10.6lf:%10.6lf]\n",E_Range[0],E_Range[1]);
	for (i=1; i<Band_Nkpath; i++)  printf("const%d = %lf\n",i,line_Nk[i]);
	for (i=1; i<Band_Nkpath; i++){
	  if (i==1){ printf("pl const%d,t \n",i); }
	  else{ printf("repl const%d,t \n",i); }
	}
	strcpy(fname_BandAll, fname_out);
	strcat(fname_BandAll, ".BAND");
	printf("rep \"%s\" lc rgb \"red\" \n", fname_BandAll);
	fp1 = fopen(fname_BandAll, "w");
	fclose(fp1);
	strcpy(fname_MP, fname_out);
	strcat(fname_MP, ".AtomMulP");
	fp1 = fopen(fname_MP, "w");
	fprintf(fp1,"                                      \n");
	fclose(fp1);
	for (i1=0; i1 <=ClaOrb_MAX[0]; i1++){
	  strcpy(fname_MP,fname_out);
	  strcat(fname_MP,".MulP_");
	  strcat(fname_MP,OrbSym[i1]);
	  fp1 = fopen(fname_MP, "w");
	  fprintf(fp1,"                                      \n");
	  fclose(fp1);
	}
	for (i1=1; i1 <=ClaOrb_MAX[1]; i1++){
	  strcpy(fname_MP,fname_out);
	  strcat(fname_MP,".MulP_");
	  strcat(fname_MP,OrbName[i1]);
	  fp1 = fopen(fname_MP, "w");
	  fprintf(fp1,"                                      \n");
	  fclose(fp1);
	}
      } 
      for (i=1; i<=Band_Nkpath; i++) {
	// ### MALLOC ARRAY ############################### 
	size_Tknum = (Band_N_perpath[i] + 1);

	k_xyz = (double**)malloc(sizeof(double*)*3);
	for(j=0; j<3; j++) k_xyz[j] = (double*)malloc(sizeof(double)*(size_Tknum+1));
	EIGEN_Band = (double*)malloc(sizeof(double)*((size_Tknum+1)*n2));
	for (j = 0; j < ((size_Tknum+1)*n2); j++) EIGEN_Band[j] = 0.0;
	MulP_Band = (double**)malloc(sizeof(double*)*4);
	for (j=0; j<4; j++){
	  MulP_Band[j] = (double*)malloc(sizeof(double)*((atomnum+1)*(size_Tknum+1)*n2));
	  for (k=0; k<((atomnum+1)*(size_Tknum+1)*n2); k++) MulP_Band[j][k] = 0.0;
	}
	OrbMulP_Band = (double***)malloc(sizeof(double**)*4);
	for (j=0; j<4; j++){
	  OrbMulP_Band[j] = (double**)malloc(sizeof(double*)*(ClaOrb_MAX[0]+1));
	  for (j1=0; j1<=ClaOrb_MAX[0]; j1++){
	    OrbMulP_Band[j][j1] = (double*)malloc(sizeof(double)*((atomnum+1)*(size_Tknum+1)*n2));
	    for (k=0; k<((atomnum+1)*(size_Tknum+1)*n2); k++) OrbMulP_Band[j][j1][k] = 0.0;
	  }//j1
	}//j
	Orb2MulP_Band = (double***)malloc(sizeof(double**)*4);
	for (j=0; j<4; j++){
	  Orb2MulP_Band[j] = (double**)malloc(sizeof(double*)*(ClaOrb_MAX[1]+1));
	  for (j1=0; j1<=ClaOrb_MAX[1]; j1++){
	    Orb2MulP_Band[j][j1] = (double*)malloc(sizeof(double)*((atomnum+1)*(size_Tknum+1)*n2));
	    for (k=0; k<((atomnum+1)*(size_Tknum+1)*n2); k++) Orb2MulP_Band[j][j1][k] = 0.0;
	  }//j1
	}//j
	// ### Division CALC_PART #########################
	T_knum = size_Tknum;
	for(j=0; j<num_procs; j++){
	  if (T_knum <= j){
	    S_knum[j] = -10;   E_knum[j] = -100;
	  } else if (T_knum < num_procs) {
	    S_knum[j] = j;     E_knum[j] = i;
	  } else {
	    d0 = (double)T_knum/(double)num_procs;
	    S_knum[j] = (int)((double)j*(d0+0.0001));
	    E_knum[j] = (int)((double)(j+1)*(d0+0.0001)) - 1;
	    if (j==(num_procs-1)) E_knum[j] = T_knum - 1;
	    if (E_knum[j]<0)      E_knum[j] = 0;
	  }
	}
	// ################################################
	for (j = S_knum[myrank]; j <= E_knum[myrank]; j++){
	  id=1;  k1 = Band_kpath[i][1][id]+(Band_kpath[i][2][id]-Band_kpath[i][1][id])*j/Band_N_perpath[i] + Shift_K_Point;
	  id=2;  k2 = Band_kpath[i][1][id]+(Band_kpath[i][2][id]-Band_kpath[i][1][id])*j/Band_N_perpath[i] - Shift_K_Point;
	  id=3;  k3 = Band_kpath[i][1][id]+(Band_kpath[i][2][id]-Band_kpath[i][1][id])*j/Band_N_perpath[i] + 2.0*Shift_K_Point;
	  k_xyz[0][j] = k1;      k_xyz[1][j] = k2;      k_xyz[2][j] = k3;
	  EigenValue_Problem(k1, k2, k3, 1);
	  for(l=1; l<=2*n; l++){
	    EIGEN_Band[j*n2+l] = EIGEN[l];
	    for (i2=0; i2 < atomnum; i2++){
	      for (j1=0; j1<4; j1++){
		//                MulP_Band[j1][j*(atomnum+1)*n2 +i2*n2 +l] = 0.0;
		for (i1=0; i1 < Total_NumOrbs[i2+1]; i1++){
		  Orb2MulP_Band[j1][ClaOrb[i2+1][i1]][j*(atomnum+1)*n2 +i2*n2 +l]+= Data_MulP[j1][l][i2+1][i1];
		  MulP_Band[j1][j*(atomnum+1)*n2 +i2*n2 +l]+= Data_MulP[j1][l][i2+1][i1];
		}//for(i1)
	      }//for(j1)
	    }//for(i2)
	  }//l
	}//j
	// ### MPI part ###################################
	for (j=0; j<num_procs; j++){
	  k = S_knum[j];          l = abs(E_knum[j]-S_knum[j]+1);
	  MPI_Barrier(MPI_COMM_WORLD);
	  MPI_Bcast(&k_xyz[0][k], l, MPI_DOUBLE, j, MPI_COMM_WORLD);
	  MPI_Barrier(MPI_COMM_WORLD);
	  MPI_Bcast(&k_xyz[1][k], l, MPI_DOUBLE, j, MPI_COMM_WORLD);
	  MPI_Barrier(MPI_COMM_WORLD);
	  MPI_Bcast(&k_xyz[2][k], l, MPI_DOUBLE, j, MPI_COMM_WORLD);
	  k = S_knum[j]*n2;       l = abs(E_knum[j]-S_knum[j]+1)*n2;
	  MPI_Barrier(MPI_COMM_WORLD);
	  MPI_Bcast(&EIGEN_Band[k], l, MPI_DOUBLE, j, MPI_COMM_WORLD);
	  k = S_knum[j]*(atomnum+1)*n2;
	  i2 = abs(E_knum[j]-S_knum[j]+1)*(atomnum+1)*n2;
	  for (j1=0; j1<4; j1++){
	    MPI_Barrier(MPI_COMM_WORLD);
	    MPI_Bcast(&MulP_Band[j1][k], i2, MPI_DOUBLE, j, MPI_COMM_WORLD);
	    for (i1=0; i1 <=ClaOrb_MAX[1]; i1++){
	      MPI_Bcast(&Orb2MulP_Band[j1][i1][k], i2, MPI_DOUBLE, j, MPI_COMM_WORLD);
	    }//for(i1)
	  }//for(j1)
	}//j
	// ################################################
	for(l=1; l<=2*n; l++){
	  for(j=0; j<size_Tknum; j++){
	    for (id=0; id<atomnum; id++){
	      for (i1=0; i1 <=ClaOrb_MAX[0]; i1++){
		for (i2=0; i2 <=ClaOrb_MAX[1]; i2++){
		  if (OrbSym[i1][0] == OrbName[i2][0]){
		    for (j1=0; j1<4; j1++){
		      OrbMulP_Band[j1][i1][j*(atomnum+1)*n2 +id*n2 +l]+= Orb2MulP_Band[j1][i2][j*(atomnum+1)*n2 +id*n2 +l];
		    }//j1
		  }//if
		}//i2
	      }//i1
	    }//id
	  }//j
	}//l
	// # SELECT CALC_BAND #############################
	id = 0;
	for(l=1; l<=2*n; l++){
	  if ((E_Range[1]-EIGEN_Band[l])*(E_Range[0]-EIGEN_Band[l])<=0){
	    if (id == 0){ l_min=l; l_max=l; id=1;}
	    else if(id == 1){ l_max=l;}
	  }
	}
	for(l=1; l<=2*n; l++){
	  for(j=0; j<size_Tknum; j++){
	    if ((E_Range[1]-EIGEN_Band[j*n2+l])*(E_Range[0]-EIGEN_Band[j*n2+l])<=0){
	      if (l_min>l){ l_min=l; }
	      else if(l_max<l){ l_max=l;}
	    }
	  } 
	}
	// # Print Band ALL ###############################
	if (myrank == 0){
	  strcpy(fname_BandAll, fname_out);  strcat(fname_BandAll, ".BAND");
	  fp1 = fopen(fname_BandAll, "a");
	  for(l=1; l<=2*n; l++){
	    for(j=0; j<size_Tknum; j++){
	      line_k=line_Nk[i-1]+kgrid_dist(Band_kpath[i][1][1],k_xyz[0][j],Band_kpath[i][1][2],k_xyz[1][j],Band_kpath[i][1][3],k_xyz[2][j]);
	      fprintf(fp1, "%10.6lf %10.6lf\n", line_k, EIGEN_Band[j*n2+l]);
	    } fprintf(fp1, "\n\n");
	  } fclose(fp1);
	  // # Print Band (division) ######################
	  /*
	     for(l=l_min; l<=l_max; l++){
	     strcpy(fname_Band,fname_out);  name_Nband(fname_Band,"_",i);  name_Nband(fname_Band,".Band",l);
	     printf("rep \"%s\" \n", fname_Band);
	     fp1= fopen(fname_Band,"w");
	     for(j=0; j<size_Tknum; j++){
	     line_k=line_Nk[i-1]+kgrid_dist(Band_kpath[i][1][1],k_xyz[0][j],Band_kpath[i][1][2],k_xyz[1][j],Band_kpath[i][1][3],k_xyz[2][j]);
	     if ((E_Range[1]-EIGEN_Band[j*n2+l])*(E_Range[0]-EIGEN_Band[j*n2+l])<=0){
	     fprintf(fp1, "%10.6lf %10.6lf\n", line_k, EIGEN_Band[j*n2+l]);
	     }
	     } fclose(fp1);
	     }//l
	     */
	  // ##############################################
	  strcpy(fname_MP,fname_out);  strcat(fname_MP,".AtomMulP");
	  fp1= fopen(fname_MP,"a");
	  for(l=1; l<=2*n; l++){
	    for(j=0; j<size_Tknum; j++){
	      if ((E_Range[1]-EIGEN_Band[j*n2+l])*(E_Range[0]-EIGEN_Band[j*n2+l])<=0){
		Print_kxyzEig(Pdata_s, k_xyz[0][j], k_xyz[1][j], k_xyz[2][j], l, EIGEN_Band[j*n2+l]);
		fprintf(fp1, "%s", Pdata_s);
		for (k=0; k<atomnum; k++){
		  fprintf(fp1,"%10.6lf %10.6lf ", MulP_Band[0][j*(atomnum+1)*n2 +k*n2 +l], MulP_Band[1][j*(atomnum+1)*n2 +k*n2 +l]);
		  fprintf(fp1,"%10.6lf %10.6lf ", MulP_Band[2][j*(atomnum+1)*n2 +k*n2 +l], MulP_Band[3][j*(atomnum+1)*n2 +k*n2 +l]);
		} fprintf(fp1,"\n"); Count_data++;
	      }//if
	    }//j
	  }//l
	  fclose(fp1);
	  // ##############################################
	  for (i1=0; i1 <=ClaOrb_MAX[0]; i1++){
	    strcpy(fname_MP,fname_out);  strcat(fname_MP,".MulP_");  strcat(fname_MP,OrbSym[i1]);
	    fp1= fopen(fname_MP,"a");
	    for(l=1; l<=2*n; l++){
	      for(j=0; j<size_Tknum; j++){
		if ((E_Range[1]-EIGEN_Band[j*n2+l])*(E_Range[0]-EIGEN_Band[j*n2+l])<=0){
		  Print_kxyzEig(Pdata_s, k_xyz[0][j], k_xyz[1][j], k_xyz[2][j], l, EIGEN_Band[j*n2+l]);
		  fprintf(fp1, "%s", Pdata_s);
		  for (k=0; k<atomnum; k++){
		    fprintf(fp1,"%10.6lf %10.6lf ", OrbMulP_Band[0][i1][j*(atomnum+1)*n2 +k*n2 +l], OrbMulP_Band[1][i1][j*(atomnum+1)*n2 +k*n2 +l]);
		    fprintf(fp1,"%10.6lf %10.6lf ", OrbMulP_Band[2][i1][j*(atomnum+1)*n2 +k*n2 +l], OrbMulP_Band[3][i1][j*(atomnum+1)*n2 +k*n2 +l]);
		  } fprintf(fp1,"\n");
		}//if
	      }//j
	    }//l
	    fclose(fp1);
	  }//i1 
	  // ##############################################
	  for (i1=1; i1 <=ClaOrb_MAX[1]; i1++){
	    strcpy(fname_MP,fname_out);  strcat(fname_MP,".MulP_");  strcat(fname_MP,OrbName[i1]);
	    fp1= fopen(fname_MP,"a");
	    for(l=1; l<=2*n; l++){
	      for(j=0; j<size_Tknum; j++){
		if ((E_Range[1]-EIGEN_Band[j*n2+l])*(E_Range[0]-EIGEN_Band[j*n2+l])<=0){
		  Print_kxyzEig(Pdata_s, k_xyz[0][j], k_xyz[1][j], k_xyz[2][j], l, EIGEN_Band[j*n2+l]);
		  fprintf(fp1, "%s", Pdata_s);
		  for (k=0; k<atomnum; k++){
		    fprintf(fp1,"%10.6lf %10.6lf ", Orb2MulP_Band[0][i1][j*(atomnum+1)*n2 +k*n2 +l], Orb2MulP_Band[1][i1][j*(atomnum+1)*n2 +k*n2 +l]);
		    fprintf(fp1,"%10.6lf %10.6lf ", Orb2MulP_Band[2][i1][j*(atomnum+1)*n2 +k*n2 +l], Orb2MulP_Band[3][i1][j*(atomnum+1)*n2 +k*n2 +l]);
		  } fprintf(fp1,"\n");
		}//if
	      }//j
	    }//l
	    fclose(fp1);
	  }//i1 
	  // ##############################################
	  // kotaka
	}//if(myrank)
	// ################################################

	// ### free ###
	for(j=0; j<3; j++) free(k_xyz[j]);
	free(k_xyz);
	free(EIGEN_Band);
	for (j=0; j<4; j++){
	  free(MulP_Band[j]);
	} free(MulP_Band);
	for (j=0; j<4; j++){
	  for (j1=0; j1<=ClaOrb_MAX[0]; j1++){
	    free(OrbMulP_Band[j][j1]);
	  } free(OrbMulP_Band[j]);
	} free(OrbMulP_Band);
	for (j=0; j<4; j++){
	  for (j1=0; j1<=ClaOrb_MAX[1]; j1++){
	    free(Orb2MulP_Band[j][j1]);
	  } free(Orb2MulP_Band[j]);
	} free(Orb2MulP_Band);
	// ################################################
      }//i
    }//if(Band_Nkpath>0)

    // #######################################################
    if (myrank == 0){
      //### atomnum & data_size ###
      printf("###########################################\n");
      printf("Total MulP data:%4d\n", Count_data); 

      strcpy(fname_MP,fname_out);      strcat(fname_MP,".AtomMulP");
      fp1= fopen(fname_MP,"r+");
      fseek(fp1, 0L, SEEK_SET);        fprintf(fp1,"%6d %4d", Count_data, atomnum);
      fclose(fp1);

      for (i1=0; i1 <=ClaOrb_MAX[0]; i1++){
	strcpy(fname_MP,fname_out);    strcat(fname_MP,".MulP_");    strcat(fname_MP,OrbSym[i1]);
	fp1= fopen(fname_MP,"r+");
	fseek(fp1, 0L, SEEK_SET);      fprintf(fp1,"%6d %4d", Count_data, atomnum);
	fclose(fp1);
      }//i1
      for (i1=1; i1 <=ClaOrb_MAX[1]; i1++){
	strcpy(fname_MP,fname_out);    strcat(fname_MP,".MulP_");    strcat(fname_MP,OrbName[i1]);
	fp1= fopen(fname_MP,"r+");
	fseek(fp1, 0L, SEEK_SET);      fprintf(fp1,"%6d %4d", Count_data, atomnum);
	fclose(fp1);
      }//i1
    }//if(myrank)
    // ### MALLOC FREE #######################################
    for (i=0; i<=ClaOrb_MAX[1]; i++){
      free(OrbName[i]);
    } free(OrbName);
    free(EIGEN); 
    // ### (MulP Calculation)   ###
    for (i=0; i<=atomnum; i++){
      free(An2Spe[i]);
    } free(An2Spe);
    for (i=0; i<=atomnum; i++){
      free(ClaOrb[i]);
    } free(ClaOrb);
    for (i=0; i<4; i++){
      for (l=0; l<n2; l++){
	for (j=0; j<=atomnum; j++){
	  free(Data_MulP[i][l][j]);
	} free(Data_MulP[i][l]);
      } free(Data_MulP[i]);
    } free(Data_MulP);
    // ### (MPI Calculation)    ###
    free(S_knum);
    free(E_knum);
    // ### (Band Calculation)   ###
    free(line_Nk);
    free(Band_N_perpath);
    if (Band_Nkpath>0) {
      for (i=0; i<(Band_Nkpath+1); i++){
	for (j=0; j<3; j++){
	  free(Band_kpath[i][j]);
	}free(Band_kpath[i]);
      }free(Band_kpath);
      for (i=0; i<(Band_Nkpath+1); i++){
	for (j=0; j<3; j++){
	  free(Band_kname[i][j]);
	}free(Band_kname[i]);
      }free(Band_kname);
    }//if(Band_Nkpath>0)
    // #######################################################
    free_scfout();
  }//if(fo_wf)

  // ### MPI_Finalize ########################################
  MPI_Finalize();

  if (myrank == 0){
    printf("############ CALC TIME ####################\n");
    dtime(&TEtime);
    printf("  Total Calculation Time:%10.6lf (s)\n",TEtime-TStime);
    printf("###########################################\n");
  }

  return 0;

}


