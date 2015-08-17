#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "image_utils.h"
#include "complex.h"

//commmand line test for executable,
// ./julia--forcebyteswap=1

void printusage(){
  fprintf(stderr,"julia <options>\n");
  fprintf(stderr,"options are:\n");
  fprintf(stderr,"--c=<a+bi|#1|#2|#3|#4|#5>    Sets the constant for iteration.  Defaults to -0.4+0.6i. #1,#2,#3,#4,#5 give nice outputs.\n");
  fprintf(stderr,"--minr=a    Sets minimum real coordinate to a. Defaults to -2.\n");
  fprintf(stderr,"--maxr=b    Sets maximum real coordinate to b.  Defaults to 2.\n");
  fprintf(stderr,"--mini=c    Sets minimum imaginary coordinate to c.  Defaults to -2.\n");
  fprintf(stderr,"--maxi=d    Sets maximum imaginary coordinate to d.  Defaults to 2.\n");
  fprintf(stderr,"--numsamples=n  Sets the number of samples in each direction to n.  Results in n x n pixel image.  Defaults to 512.\n");
  fprintf(stderr,"--maxiterations=n  Sets number of iterations to be considered \"surely in Mandelbrot set\" to n.  Defaults to 500 unless monochrome output is selected, in which case it is 250.\n");
  fprintf(stderr,"--maxnorm=M    Sets cutoff on norm of iterations to M.  Defaults to 2.\n");
  fprintf(stderr,"--outputfile=filename  Sets the output file name.  Defaults to mandelbrot. Proper extension is appended.\n");
  fprintf(stderr,"--outputformat=<ppm|monochrome|greyscale|rgb>  Sets the output format.  ppm makes a ppm file.  monochrome makes a monochrome .bmp file.  greyscale makes an 8 bit greyscale bmp.  rgb makes a 24 bit rgb bmp. Default is greyscale\n");
  fprintf(stderr,"--forcebyteswap=<0|1>  Forces output file to be written in with endian-ness reversed.  Defaults to 0.  Debugging purposes only.");
  fprintf(stderr,"-h print this message.\n");
}

void set_opts(int argc, char** argv,double* minr, double* maxr, double* mini, double* maxi, int* numsamples, long int* maxiterations, double* maxnorm, char** outputfile, int* outputformat,complex* c,
	      int* forcebyteswap){
  int i;
  double dummyd;
  double dummyd1;
  long int dummyli;
  int dummyi;
  char dummystr[1000];
  *minr=-2;
  *maxr=2;
  *mini=-2;
  *maxi=2;
  *numsamples=512;
  *maxiterations=500;
  *maxnorm=2;
  *outputfile=malloc(1000*sizeof(char));
  strcpy(*outputfile,"julia");
  *outputformat=2;
  c->real=-0.4;
  c->imag=0.6;
  *forcebyteswap=0;
  
  int maxiteroverride=0;
  
  for (i=1;i<argc;i++){
    if (sscanf(argv[i],"--c=%lf+%lfi",&dummyd,&dummyd1)==2){
      c->real=dummyd;
      c->imag=dummyd1;
    }
    else if(sscanf(argv[i],"--c=%lf-%lfi",&dummyd,&dummyd1)==2){
      c->real=dummyd;
      c->imag=-dummyd1;
    }
    else if(sscanf(argv[i],"--c=%lf",&dummyd)==1 && argv[i][strlen(argv[i])-1]=='i'){
      c->real=0;
      c->imag=dummyd;
    }
    else if(sscanf(argv[i],"--c=%lf",&dummyd)==1){
      c->real=dummyd;
      c->imag=0;
    }
    else if(strcmp(argv[i],"--c=#1")==0){
      c->real=-0.4;
      c->imag=0.6;
    }
    else if(strcmp(argv[i],"--c=#2")==0){
      c->real=-.618034;
      c->imag=0;
    }
    else if(strcmp(argv[i],"--c=#3")==0){
      c->real=-0.123;
      c->imag=0.745;
    }
    else if(strcmp(argv[i],"--c=#4")==0){
      c->real=-1;
      c->imag=0;
    }
    else if(strcmp(argv[i],"--c=#5")==0){
      c->real=-0.8;
      c->imag=0.156;
    }
    
    else if (sscanf(argv[i],"--minr=%lf",&dummyd)==1){
      *minr=dummyd;
    }
    else if (sscanf(argv[i],"--maxr=%lf",&dummyd)==1){
      *maxr=dummyd;
    }
    else if (sscanf(argv[i],"--mini=%lf",&dummyd)==1){
      *mini=dummyd;
    }
    else if (sscanf(argv[i],"--maxi=%lf",&dummyd)==1){
      *maxi=dummyd;
    }
    else if (sscanf(argv[i],"--numsamples=%d",&dummyi)==1){
      *numsamples=dummyi;
    }
    else if (sscanf(argv[i],"--maxiterations=%ld",&dummyli)==1){
      *maxiterations=dummyli;
      maxiteroverride=1;
    }
    else if (sscanf(argv[i],"--maxnorm=%lf",&dummyd)==1){
      *maxnorm=dummyd;
    }
    else if (sscanf(argv[i],"--outputfile=%s",dummystr)==1){
      strcpy(*outputfile,dummystr);
    }
    else if (strcmp(argv[i],"--outputformat=ppm")==0){
      *outputformat=0;
    }
    else if (strcmp(argv[i],"--outputformat=monochrome")==0){
      *outputformat=1;
    }
    else if (strcmp(argv[i],"--outputformat=greyscale")==0){
      *outputformat=2;
    }
    else if (strcmp(argv[i],"--outputformat=rgb")==0){
      *outputformat=3;
    }
    else if (strcmp(argv[i],"-h")==0){
      printusage();
      exit(0);
    }
    else if(sscanf(argv[i],"--forcebyteswap=%d",&dummyi)==1){
      *forcebyteswap=dummyi;
    }
    else{
      fprintf(stderr,"Unrecognized option %s\n",argv[i]);
      printusage();
      exit(1);
    }
  }

  switch (*outputformat){
  case 0:
    strcat(*outputfile,".ppm");
    break;
  default:
    strcat(*outputfile,".bmp");
  }

  if (*outputformat==1 && maxiteroverride==0){
    *maxiterations=250;
  }
  

}

int main(int argc, char** argv){

  char* outputfile;
  int outputformat;
  double rmin;
  double rmax;
  double imin;
  double imax;
  long int maxiter;
  int numsamples;
  double maxnorm;

  complex c;

  int forcebyteswap;
  set_opts(argc,argv,&rmin, &rmax, &imin, &imax, &numsamples, &maxiter, &maxnorm, &outputfile, &outputformat,&c,&forcebyteswap);
  
  double* r,*g,*b;
  short int* m;

  switch(outputformat){
  case 0:
    g=malloc(numsamples*numsamples*sizeof(double));
    break;
  case 1:
    m=malloc(numsamples*numsamples*sizeof(short int));
    break;
  case 2:
    g=malloc(numsamples*numsamples*sizeof(double));
    break;
  case 3:
    r=malloc(numsamples*numsamples*sizeof(double));
    g=malloc(numsamples*numsamples*sizeof(double));
    b=malloc(numsamples*numsamples*sizeof(double));
  }
  
  // escape-time algorithm resource
  // https://books.google.com/books?id=WTQZthj7JQoC&pg=PA175&lpg=PA175&dq=julia+set+escape+time+algorithm+c&source=bl&ots=o5smuboavo&sig=xDBQwtt-LeKuvTqYhCvY70iNtRE&hl=en&sa=X&ei=-XM0VbnIFYrmoAS1hICgBA&ved=0CEAQ6AEwBA#v=onepage&q=julia%20set%20escape%20time%20algorithm%20c&f=false

  int i,j;
  double deltar=(rmax-rmin)/numsamples;
  double deltai=(imax-imin)/numsamples;

  printf("deltai:%f\ndeltar:%f\nmaxiter:%ld\nmaxnorm:%lf\noutput file: %s\nc:%f+%fi\n",deltai,deltar,maxiter,maxnorm,outputfile,c.real,c.imag);

  switch (outputformat){
  case 0:
    for(i=0; i<numsamples; i++)
    {
      for(j=0; j<numsamples; j++)
      {
        int k;
        complex oldComp;
        oldComp.real = rmin + j*deltar;
        oldComp.imag = imin + i*deltai;
        for(k=0; k<maxiter; k++)
        {
          complex newComp = add(mulitply(oldComp, oldComp), c);
          double squared = squaredNorm(newComp);
          if(squared > maxnorm*maxnorm)
          {
            m[i*numsamples+j] = 0;
            break;
          }
          else
          {
            m[i*numsamples+j] = 1;
          }
          oldComp.real = newComp.real;
          oldComp.imag = newComp.imag;
        }
      }
    }
    write_ppm(outputfile,g,numsamples,numsamples,0,maxiter);
    break;
  case 1:
    for(i=0; i<numsamples; i++)
    {
      for(j=0; j<numsamples; j++)
      {
        int k;
        complex oldComp;
        oldComp.real = rmin + j*deltar;
        oldComp.imag = imin + i*deltai;
        for(k=0; k<maxiter; k++)
        {
          complex newComp = add(mulitply(oldComp, oldComp), c);
          double squared = squaredNorm(newComp);
          if(squared > maxnorm*maxnorm)
          {
            m[i*numsamples+j] = 0;
            break;
          }
          else
          {
            m[i*numsamples+j] = 1;
          }
          oldComp.real = newComp.real;
          oldComp.imag = newComp.imag;
        }
      }
    }
    write_monochrome_bmp(outputfile,m,numsamples,numsamples,forcebyteswap);
    break;
  case 2:
    for(i=0; i<numsamples; i++)
    {
      for(j=0; j<numsamples; j++)
      {
        int k;
        complex oldComp;
        oldComp.real = rmin + j*deltar;
        oldComp.imag = imin + i*deltai;
        for(k=0; k<maxiter; k++)
        {
          complex newComp = add(mulitply(oldComp, oldComp), c);
          double squared = squaredNorm(newComp);
          if(squared > maxnorm*maxnorm)
          {
            g[i*numsamples+j] = k;
            break;
          }
          else if(squared == maxnorm*maxnorm)
          {
            g[i*numsamples+j] = 0;
          }
          else
          {
            g[i*numsamples+j] = maxiter;
          }
          oldComp.real = newComp.real;
          oldComp.imag = newComp.imag;
        }
      }
    }
    write_greyscale_bmp(outputfile,g,numsamples,numsamples,0,maxiter,forcebyteswap);
    break;
  case 3:
    for(i=0; i<numsamples; i++)
    {
      for(j=0; j<numsamples; j++)
      {
        int k;
        complex oldComp;
        oldComp.real = rmin + j*deltar;
        oldComp.imag = imin + i*deltai;
        for(k=0; k<maxiter; k++)
        {
          complex newComp = add(mulitply(oldComp, oldComp), c);
          double squared = squaredNorm(newComp);
          if(squared > maxnorm*maxnorm)
          {
            r[i*numsamples+j] = 0;
            g[i*numsamples+j] = k;
            b[i*numsamples+j] = maxiter-k;
            break;
          }
          else if(squared == maxnorm*maxnorm)
          {
            r[i*numsamples+j] = 0;
            g[i*numsamples+j] = 0;
            b[i*numsamples+j] = 0;
          }
          else
          {
            r[i*numsamples+j] = maxiter;
            g[i*numsamples+j] = maxiter;
            b[i*numsamples+j] = maxiter;
          }
          oldComp.real = newComp.real;
          oldComp.imag = newComp.imag;
        }
      }
    }
    write_rgb_bmp(outputfile,r,g,b,numsamples,numsamples,0,maxiter,forcebyteswap);
  }

  free(outputfile);

  switch (outputformat){
  case 0:
    free(g);
    break;
  case 1:
    free(m);
    break;
  case 2:
    free(g);
    break;
  case 3:
    free(g);free(r);free(b);
    break;
  }
  
  
  return 0;
}