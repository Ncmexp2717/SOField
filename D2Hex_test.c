#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

char D2Hex(int i){
  if(i<0) return 0;
  else if(i<10) return i+48; 
  else if(i<16) return i+55; 
  else return 0;
}

int main(){

  int i,j,k;
  int iR, iG, iB;
  char ColorCode[7];
  
  for(i=0;i<100;i+=4){
    iR = 128;
    iG = i*2.56;
    iB = 255;
    ColorCode[0] = D2Hex(iR/16);    ColorCode[1] = D2Hex(iR%16);
    ColorCode[2] = D2Hex(iG/16);    ColorCode[3] = D2Hex(iG%16);
    ColorCode[4] = D2Hex(iB/16);    ColorCode[5] = D2Hex(iB%16);
    ColorCode[6] = '\0';
    printf ("%s\n",  ColorCode);
  }
/*
  for(j =0; j<16; j++){
    printf ("%c\n",  D2Hex(j));
  }
*/  
  return 0;
}


