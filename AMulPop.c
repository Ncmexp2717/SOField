#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Inputtools.h"
#include "Tools_BandCalc.h"


int main(int argc, char *argv[]){
  FILE *fp,*fp1,*fp2;
  char c;
  int i,j,k,l;
  double d0,d1,d2,d3;

  int Atom_Num, Data_size;
  int l_sum, l_max, l_min, *Data_l;
  double Pxyz[3], MulP_VecScale[3];
  int Num_EA, switch_deg, Num_EB;
  int *EA, *Extract_Atom;
  int Data_Reduction;

  double **Data_k, *Data_kpath, *Data_E, ***Data_atomMulP;
  double **sum_atomMulP;
  double Re11, Re22, Re12, Im12;                  // MulP Calc variable-1
  double Nup[2], Ndw[2], Ntheta[2], Nphi[2];      // MulP Calc variable-2

  int i_vec[20];           // input variable
  int *i_vec2;             // input variable
  char *s_vec[20];         // input variable
  double r_vec[20];        // input variable

  int fl_inp = 0;
  int fl_data  = 0;
  char fname[256], fname_1[256], fname_aMulP[256], fname_out[256];

  //############ INPUT DATA #########################
  sprintf(fname,"%s",argv[1]);
  if((fp = fopen(fname,"r")) != NULL){
    fl_inp = 1;
    printf("open \"%s\" file... \n" ,fname);
    fclose(fp);
  }else{
    fl_inp = 0;
    printf("Cannot open \"%s\" file.\n" ,fname);
  }

  if (fl_inp == 1){
    input_open(fname);
    input_string("Filename.atomMulP",fname_aMulP,"default");
    input_string("Filename.xyzdata",fname_out,"default");
    r_vec[0]=1.0; r_vec[1]=1.0; r_vec[2]=1.0;
    input_doublev("MulP.Vec.Scale",3, MulP_VecScale, r_vec);
    input_int("Num.of.Extract.Atom",&Num_EA,1);
    if (Num_EA > 0){
      i_vec2 = (int*)malloc(sizeof(int)*(Num_EA+1));
      Extract_Atom = (int*)malloc(sizeof(int)*(Num_EA+1));
      for(i=0;i<Num_EA;i++) i_vec2[i]=i+1;
      input_intv("Extract.Atom",Num_EA,Extract_Atom,i_vec2);
      free(i_vec2);
    }
    input_int("Data.Reduction",&Data_Reduction,1);
    input_close();
  }
  
  
  if((fp = fopen(fname_aMulP,"r")) != NULL){ 
    fl_data = 1;
    fscanf(fp,"%d", &Data_size);
    fscanf(fp,"%d", &Atom_Num);
    EA = (int*)malloc(sizeof(int)*(Atom_Num+1));
    for(i=0;i<=Atom_Num;i++) EA[i] = 0;
    if (Num_EA > 0){
      printf("Calculate Atom:");
      for(i=0;i<Num_EA;i++){
        /*kotaka*/
        if (Extract_Atom[i]>0 && Extract_Atom[i]<=Atom_Num){
          EA[Extract_Atom[i]-1] = 1;
          printf("%4d ",Extract_Atom[i]);
        }else{
          printf("%4d is illegal value.\n",Extract_Atom[i]);
        }     
      } printf("\n");
      printf("MulP Scale:(%10.6lf %10.6lf %10.6lf)\n",MulP_VecScale[0],MulP_VecScale[1],MulP_VecScale[2]);  
    }//if
    Data_k = (double**)malloc(sizeof(double*)*(Data_size+1));
    for(j=0;j<=Data_size;j++){
      Data_k[j] = (double*)malloc(sizeof(double)*3);
    }
    Data_l = (int*)malloc(sizeof(int)*(Data_size+1));
    Data_E = (double*)malloc(sizeof(double)*(Data_size+1));
    Data_atomMulP = (double***)malloc(sizeof(double**)*(Data_size+1));
    for(j=0;j<=Data_size;j++){
      Data_atomMulP[j] = (double**)malloc(sizeof(double*)*(Atom_Num));
      for(i=0;i<Atom_Num;i++){
      Data_atomMulP[j][i] = (double*)malloc(sizeof(double)*4);
      }
    }
    for(j=0;j<Data_size;j++){
      for(i=0;i<3;i++) fscanf(fp,"%lf",&Data_k[j][i]);
      fscanf(fp,"%d",&Data_l[j]);
      fscanf(fp,"%lf",&Data_E[j]);
      for(i=0;i<Atom_Num;i++){
        for(k=0;k<4;k++) fscanf(fp,"%lf",&Data_atomMulP[j][i][k]);
      }//i 
    }//j

    fclose(fp);
    l_max = Data_l[0];
    l_min = Data_l[0];
    for(j=0;j<Data_size;j++){
      if (l_min>Data_l[j]) l_min=Data_l[j];
      if (l_max<Data_l[j]) l_max=Data_l[j];
    }
    sum_atomMulP = (double**)malloc(sizeof(double*)*(Data_size+1));
    for(j=0;j<=Data_size;j++){
      sum_atomMulP[j] = (double*)malloc(sizeof(double)*7);
    }
  }


  //############ PRINT PART #########################
  for(j=0;j<Data_size;j++){
    sum_atomMulP[j][0] = 0.0;  sum_atomMulP[j][1] = 0.0;  sum_atomMulP[j][2] = 0.0;
    sum_atomMulP[j][3] = 0.0;  sum_atomMulP[j][4] = 0.0;  sum_atomMulP[j][5] = 0.0;  sum_atomMulP[j][6] = 0.0;
  }  k = 0;

  strcpy(fname_1,fname_out);
  strcat(fname_1,".MulPop");
  
  fp1 = fopen(fname_1,"w");
  fprintf(fp1,"# kx[Ang-1]   ky[Ang-1]   kz[Ang-1]   Eig[eV]     Num_Ele     Spin_x      Spin_y      Spin_z  \n");

  for(l=l_min;l<=l_max;l++){
    strcpy(fname_1,fname_out);
    strcat(fname_1,".MulPop");
    name_Nband(fname_1,"_",l);
    
    fp = fopen(fname_1,"w");
    fprintf(fp,"# kx[Ang-1]   ky[Ang-1]   kz[Ang-1]   Eig[eV]     Num_Ele     Spin_x      Spin_y      Spin_z  \n");

    for(j=0;j<Data_size;j++){
      if (Data_l[j] == l){
        k++;
        sum_atomMulP[j][0]=Data_k[j][0];  sum_atomMulP[j][1]=Data_k[j][1];  sum_atomMulP[j][2]=Data_k[j][2];
        sum_atomMulP[j][3]=0;      sum_atomMulP[j][4]=0;      sum_atomMulP[j][5]=0;      sum_atomMulP[j][6]=0;
        for(i=0;i<Atom_Num;i++){
          if (EA[i]==1){
            sum_atomMulP[j][3]+= Data_atomMulP[j][i][0]; //Re11
            sum_atomMulP[j][4]+= Data_atomMulP[j][i][1];  //Re22
            sum_atomMulP[j][5]+= Data_atomMulP[j][i][2];  //Re12
            sum_atomMulP[j][6]+= Data_atomMulP[j][i][3];  //Im12
          }//if
        }//i
        if (k%Data_Reduction == 0  || Data_Reduction == 0){
          Re11 = sum_atomMulP[j][3];       Re22 = sum_atomMulP[j][4];
          Re12 = sum_atomMulP[j][5];       Im12 = sum_atomMulP[j][6];
          EulerAngle_Spin( 1, Re11, Re22, Re12, Im12, Re12, -Im12, Nup, Ndw, Ntheta, Nphi );
          fprintf(fp, "%10.6lf  " ,sum_atomMulP[j][0]/BohrR);
          fprintf(fp, "%10.6lf  " ,sum_atomMulP[j][1]/BohrR);
          fprintf(fp, "%10.6lf  " ,sum_atomMulP[j][2]/BohrR);
          fprintf(fp, "%10.6lf  ", Data_E[j]);
          fprintf(fp, "%10.6lf  " ,Nup[0] +Ndw[0]);
          fprintf(fp, "%10.6lf  " ,(Nup[0] -Ndw[0]) *sin(Ntheta[0]) *cos(Nphi[0]) *MulP_VecScale[0] );
          fprintf(fp, "%10.6lf  " ,(Nup[0] -Ndw[0]) *sin(Ntheta[0]) *sin(Nphi[0]) *MulP_VecScale[1] );
          fprintf(fp, "%10.6lf\n" ,(Nup[0] -Ndw[0]) *cos(Ntheta[0]) *MulP_VecScale[2] );

          fprintf(fp1, "%10.6lf  " ,sum_atomMulP[j][0]/BohrR);
          fprintf(fp1, "%10.6lf  " ,sum_atomMulP[j][1]/BohrR);
          fprintf(fp1, "%10.6lf  " ,sum_atomMulP[j][2]/BohrR);
          fprintf(fp1, "%10.6lf  ", Data_E[j]);
          fprintf(fp1, "%10.6lf  " ,Nup[0] +Ndw[0]);
          fprintf(fp1, "%10.6lf  " ,(Nup[0] -Ndw[0]) *sin(Ntheta[0]) *cos(Nphi[0]) *MulP_VecScale[0] );
          fprintf(fp1, "%10.6lf  " ,(Nup[0] -Ndw[0]) *sin(Ntheta[0]) *sin(Nphi[0]) *MulP_VecScale[1] );
          fprintf(fp1, "%10.6lf\n" ,(Nup[0] -Ndw[0]) *cos(Ntheta[0]) *MulP_VecScale[2] );
          }//if(Reduction)
      }//if(l)
    }//j
    fclose(fp);
  }//l
  fclose(fp1);

  //############ FREE MALLOC ########################
  if (fl_inp == 1){
    free(Extract_Atom);
  }
  if ((fl_inp == 1) && (fl_data == 1)){
    free(EA);
    for(j=0;j<=Data_size;j++){
      free(sum_atomMulP[j]);
    } free(sum_atomMulP);
    for(j=0;j<=Data_size;j++)
      free(Data_k[j]);
    free(Data_k);
    free(Data_l);
    free(Data_E);
    for(j=0;j<=Data_size;j++){
      for(i=0;i<Atom_Num;i++){
        free(Data_atomMulP[j][i]);
      }free(Data_atomMulP[j]);
    }free(Data_atomMulP);
  }//if
  //############ FREE MALLOC ########################

}



