#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "extern.h"

double *Jw;
double *Jwpp;
double *Jwxp;
double *Jwpx;
double *Jtheta;
double *Tw;
double *Ttheta;
double *Fw;
double *Ftheta;
double *Fw_prev;
double *Ftheta_prev;
double *fx;
double *fy;
double *p;
double *q;
double *alphax;
double *betax;
double *alphay;
double *betay;
double *psi1;
double *psi2;
double *psi3;
double *tauxtab;
double *tauytab;
int M;

void teascheme(void)
{
  int i,j;
  long k;
  k=-1;
  for(j=0;j<numycells;j++)
    for(i=0;i<numxcells;i++)
    {
      k++;
      if((i!=0) && (i!=numxcells-1) && (j!=0) && (j!=numycells-1))
      {
	Jwpp[k]=0.25*((cells[k+1].psi-cells[k-1].psi)*(cells[k+numxcells].w-cells[k-numxcells].w)-
	    (cells[k+numxcells].psi-cells[k-numxcells].psi)*(cells[k+1].w-cells[k-1].w))/(hx*hy);
      }
    }
  for(j=1;j<numycells-1;j++)
  {
    Jwpp[j*numxcells]=0.25*((cells[j*numxcells+1].psi-cells[(j+1)*numxcells-1].psi)*
	(cells[(j+1)*numxcells].w-cells[(j-1)*numxcells].w)-
      (cells[(j+1)*numxcells].psi-cells[(j-1)*numxcells].psi)*
      (cells[j*numxcells+1].w-cells[(j+1)*numxcells-1].w))/(hx*hy);
    Jwpp[(j+1)*numxcells-1]=0.25*((cells[j*numxcells].psi-cells[(j+1)*numxcells-2].psi)*
	(cells[(j+2)*numxcells-1].w-cells[j*numxcells-1].w)-
      (cells[(j+2)*numxcells-1].psi-cells[j*numxcells-1].psi)*
      (cells[j*numxcells].w-cells[(j+1)*numxcells-2].w))/(hx*hy);
  }
  k=-1;
  for(j=0;j<numycells;j++)
    for(i=0;i<numxcells;i++)
    {
      k++;
      if((i!=0) && (i!=numxcells-1) && (j!=0) && (j!=numycells-1))
      {
	Jwxp[k]=0.25*(cells[k+1].psi*(cells[k+numxcells+1].w-cells[k-numxcells+1].w)-
	    cells[k-1].psi*(cells[k+numxcells-1].w-cells[k-numxcells-1].w)-
	    cells[k+numxcells].psi*(cells[k+numxcells+1].w-cells[k+numxcells-1].w)+
	    cells[k-numxcells].psi*(cells[k-numxcells+1].w-cells[k-numxcells-1].w))/(hx*hy);
      }
    }
  for(j=1;j<numycells-1;j++)
  {
    Jwxp[j*numxcells]=0.25*(cells[j*numxcells+1].psi*(cells[(j+1)*numxcells+1].w-cells[(j-1)*numxcells+1].w)-
      cells[(j+1)*numxcells-1].psi*(cells[(j+2)*numxcells-1].w-cells[j*numxcells-1].w)-
      cells[(j+1)*numxcells].psi*(cells[(j+1)*numxcells+1].w-cells[(j+2)*numxcells-1].w)+
      cells[(j-1)*numxcells].psi*(cells[(j-1)*numxcells+1].w-cells[j*numxcells-1].w))/(hx*hy);
    Jwxp[(j+1)*numxcells-1]=0.25*(cells[j*numxcells].psi*(cells[(j+1)*numxcells].w-cells[(j-1)*numxcells].w)-
      cells[(j+1)*numxcells-2].psi*(cells[(j+2)*numxcells-2].w-cells[j*numxcells-2].w)-
      cells[(j+2)*numxcells-1].psi*(cells[(j+1)*numxcells].w-cells[(j+2)*numxcells-2].w)+
      cells[j*numxcells-1].psi*(cells[(j-1)*numxcells].w-cells[j*numxcells-2].w))/(hx*hy);
  }
  k=-1;
  for(j=0;j<numycells;j++)
    for(i=0;i<numxcells;i++)
    {
      k++;
      if((i!=0) && (i!=numxcells-1) && (j!=0) && (j!=numycells-1))
      {
	Jwpx[k]=0.25*(-cells[k+1].w*(cells[k+numxcells+1].psi-cells[k-numxcells+1].psi)+
	    cells[k-1].w*(cells[k+numxcells-1].psi-cells[k-numxcells-1].psi)+
	    cells[k+numxcells].w*(cells[k+numxcells+1].psi-cells[k+numxcells-1].psi)-
	    cells[k-numxcells].w*(cells[k-numxcells+1].psi-cells[k-numxcells-1].psi))/(hx*hy);
      }
    }
  for(j=1;j<numycells-1;j++)
  {
    Jwpx[j*numxcells]=0.25*(-cells[j*numxcells+1].w*(cells[(j+1)*numxcells+1].psi-cells[(j-1)*numxcells+1].psi)+
	cells[(j+1)*numxcells-1].w*(cells[(j+2)*numxcells-1].psi-cells[j*numxcells-1].psi)+
	cells[(j+1)*numxcells].w*(cells[(j+1)*numxcells+1].psi-cells[(j+2)*numxcells-1].psi)-
	cells[(j-1)*numxcells].w*(cells[(j-1)*numxcells+1].psi-cells[j*numxcells-1].psi))/(hx*hy);
    Jwpx[(j+1)*numxcells-1]=0.25*(-cells[j*numxcells].w*(cells[(j+1)*numxcells].psi-cells[(j-1)*numxcells].psi)+
	cells[(j+1)*numxcells-2].w*(cells[(j+2)*numxcells-2].psi-cells[j*numxcells-2].psi)+
	cells[(j+2)*numxcells-1].w*(cells[(j+1)*numxcells].psi-cells[(j+2)*numxcells-2].psi)-
	cells[j*numxcells-1].w*(cells[(j-1)*numxcells].psi-cells[j*numxcells-2].psi))/(hx*hy);
  }
  k=-1;
  for(j=0;j<numycells;j++)
    for(i=0;i<numxcells;i++)
    {
      k++;
      if((j!=0) && (j!=numycells-1))
      {
//	Jw[k]=Jwpp[k];
	Jw[k]=(Jwpp[k]+Jwxp[k]+Jwpx[k])/3;
      }
    }
  k=-1;
  for(j=0;j<numycells;j++)
    for(i=0;i<numxcells;i++)
    {
      k++;
      if((i!=0) && (i!=numxcells-1) && (j!=0) && (j!=numycells-1))
      {
	Jtheta[k]=-0.5*((cells[k].theta-cells[k-1].theta)*(cells[k-1].psi-cells[k-numxcells-1].psi)+
	    (cells[k+1].theta-cells[k].theta)*(cells[k].psi-cells[k-numxcells].psi)-
	    (cells[k].theta-cells[k-numxcells].theta)*(cells[k-numxcells].psi-cells[k-numxcells-1].psi)-
	    (cells[k+numxcells].theta-cells[k].theta)*(cells[k].psi-cells[k-1].psi))/(hx*hy);
      }
    }
  for(j=1;j<numycells-1;j++)
  {
    Jtheta[j*numxcells]=-0.5*((cells[j*numxcells].theta-cells[(j+1)*numxcells-1].theta)*
	(cells[(j+1)*numxcells-1].psi-cells[j*numxcells-1].psi)+
      (cells[j*numxcells+1].theta-cells[j*numxcells].theta)*
      (cells[j*numxcells].psi-cells[(j-1)*numxcells].psi)-
      (cells[j*numxcells].theta-cells[(j-1)*numxcells].theta)*
      (cells[(j-1)*numxcells].psi-cells[j*numxcells-1].psi)-
      (cells[(j+1)*numxcells].theta-cells[j*numxcells].theta)*
      (cells[j*numxcells].psi-cells[(j+1)*numxcells-1].psi))/(hx*hy);
    Jtheta[(j+1)*numxcells-1]=-0.5*((cells[(j+1)*numxcells-1].theta-cells[(j+1)*numxcells-2].theta)*
	(cells[(j+1)*numxcells-2].psi-cells[j*numxcells-2].psi)+
      (cells[j*numxcells].theta-cells[(j+1)*numxcells-1].theta)*
      (cells[(j+1)*numxcells-1].psi-cells[j*numxcells-1].psi)-
      (cells[(j+1)*numxcells-1].theta-cells[j*numxcells-1].theta)*
      (cells[j*numxcells-1].psi-cells[j*numxcells-2].psi)-
      (cells[(j+2)*numxcells-1].theta-cells[(j+1)*numxcells-1].theta)*
      (cells[(j+1)*numxcells-1].psi-cells[(j+1)*numxcells-2].psi))/(hx*hy);
  }
  k=-1;
  for(j=0;j<numycells;j++)
    for(i=0;i<numxcells;i++)
    {
      k++;
      if((i!=0) && (i!=numxcells-1) && (j!=0) && (j!=numycells-1))
      {
	Tw[k]=K*((cells[k+1].w-2*cells[k].w+cells[k-1].w)/(hx*hx)+
	    (cells[k+numxcells].w-2*cells[k].w+cells[k-numxcells].w)/(hy*hy));
      }
    }
  for(j=1;j<numycells-1;j++)
  {
    Tw[j*numxcells]=K*((cells[j*numxcells+1].w-2*cells[j*numxcells].w+cells[(j+1)*numxcells-1].w)/(hx*hx)+
	(cells[(j+1)*numxcells].w-2*cells[j*numxcells].w+cells[(j-1)*numxcells].w)/(hy*hy));
    Tw[(j+1)*numxcells-1]=K*((cells[j*numxcells].w-2*cells[(j+1)*numxcells-1].w+cells[(j+1)*numxcells-2].w)/(hx*hx)+
	(cells[(j+2)*numxcells-1].w-2*cells[(j+1)*numxcells-1].w+cells[j*numxcells-1].w)/(hy*hy));
  }
  k=-1;
  for(j=0;j<numycells;j++)
    for(i=0;i<numxcells;i++)
    {
      k++;
      if((i!=0) && (i!=numxcells-1) && (j!=0) && (j!=numycells-1))
      {
	Ttheta[k]=K*((cells[k+1].theta-2*cells[k].theta+cells[k-1].theta)/(hx*hx)+
	    (cells[k+numxcells].theta-2*cells[k].theta+cells[k-numxcells].theta)/(hy*hy));
      }
    }
  for(j=1;j<numycells-1;j++)
  {
    Ttheta[j*numxcells]=K*((cells[j*numxcells+1].theta-2*cells[j*numxcells].theta+
	  cells[(j+1)*numxcells-1].theta)/(hx*hx)+
      (cells[(j+1)*numxcells].theta-2*cells[j*numxcells].theta+cells[(j-1)*numxcells].theta)/(hy*hy));
    Ttheta[(j+1)*numxcells-1]=K*((cells[j*numxcells].theta-2*cells[(j+1)*numxcells-1].theta+
	  cells[(j+1)*numxcells-2].theta)/(hx*hx)+
      (cells[(j+2)*numxcells-1].theta-2*cells[(j+1)*numxcells-1].theta+cells[j*numxcells-1].theta)/(hy*hy));
  }
  k=-1;
  for(j=0;j<numycells;j++)
    for(i=0;i<numxcells;i++)
    {
      k++;
      if((i!=numxcells-1) && (j!=0) && (j!=numycells-1))
      {
	Fw[k]=Jw[k]+0.5*beta*(cells[k+1].theta-cells[k].theta+
	    cells[k+numxcells+1].theta-cells[k+numxcells].theta)/hx+Tw[k];
      }
    }
  for(j=1;j<numycells;j++)
  {
    Fw[(j+1)*numxcells-1]=Jw[(j+1)*numxcells-1]+0.5*beta*(cells[j*numxcells].theta-cells[(j+1)*numxcells-1].theta+
      cells[(j+1)*numxcells].theta-cells[(j+2)*numxcells-1].theta)/hx+Tw[k];
  }
  k=-1;
  for(j=0;j<numycells;j++)
    for(i=0;i<numxcells;i++)
    {
      k++;
      if((j!=0) && (j!=numycells-1))
      {
	Ftheta[k]=Jtheta[k]+Ttheta[k];
      }
    }
  if(step!=1)
  {
    k=-1;
    for(j=0;j<numycells;j++)
      for(i=0;i<numxcells;i++)
      {
	k++;
	if((j!=0) && (j!=numycells-1))
	{
	  cells2[k].w=cells[k].w+dt*(1.5*Fw[k]-0.5*Fw_prev[k]);
	  cells2[k].theta=cells[k].theta+dt*(1.5*Ftheta[k]-0.5*Ftheta_prev[k]);
	}
      }
  }
  else
  {
    k=-1;
    for(j=0;j<numycells;j++)
      for(i=0;i<numxcells;i++)
      {
	k++;
	if((j!=0) && (j!=numycells-1))
	{
	  cells2[k].w=cells[k].w+dt*Fw[k];
	  cells2[k].theta=cells[k].theta+dt*Ftheta[k];
	}
      }
  }
  k=-1;
  for(j=0;j<numycells;j++)
    for(i=0;i<numxcells;i++)
    {
      k++;
      if((j!=0) && (j!=numycells-1))
      {
	Fw_prev[k]=Fw[k];
	Ftheta_prev[k]=Ftheta[k];
      }
    }

  // Решение эллиптического уравнения

  double taux;
  double tauy;

  double a;
  int m;

//  temp1=-2*sqrt(((double)numxcells)*((double)numycells))/3.14*log(epsilonpsi);
//  M=(int)(temp1);


  for(m=0;m<M;m++)
  {
//    taux=1/(Lamb1+lamb1);
//    tauy=1/(Lamb2+lamb2);
    taux=tauxtab[m];
    tauy=tauytab[m];

    k=-1;
    for(j=0;j<numycells;j++)
      for(i=0;i<numxcells;i++)
      {
	k++;
	psi1[k]=cells[k].psi;
      }

    a=2+hx*hx/taux;

    // Определение p
    alphax[1]=1/a;
    betax[1]=-1/a;
    for(i=2;i<numxcells-1;i++)
    {
      alphax[i]=-1/(alphax[i-1]+a);
      betax[i]=betax[i-1]/(alphax[i-1]+a);
    }
    p[numxcells-2]=(betax[numxcells-2]-1)/(alphax[numxcells-2]+a);
    for(i=numxcells-3;i>=0;i--)
    {
      p[i]=-alphax[i+1]*p[i+1]+betax[i+1];
    }

    for(j=1;j<numycells-1;j++)
    {
      for(i=1;i<numxcells-1;i++)
      {
	fx[i]=-hx*hx*(cells[j*numxcells+i].w+psi1[j*numxcells+i]/taux+
	    (psi1[j*numxcells+i-1]-2*psi1[j*numxcells+i]+psi1[j*numxcells+i+1])/(hy*hy));
      }
      fx[0]=-hx*hx*(cells[j*numxcells].w+psi1[j*numxcells]/taux+
	  (psi1[(j+1)*numxcells]-2*psi1[j*numxcells]+psi1[(j-1)*numxcells])/(hy*hy));
      fx[numxcells-1]=-hx*hx*(cells[(j+1)*numxcells-1].w+psi1[(j+1)*numxcells-1]/taux+
	  (psi1[(j+2)*numxcells-1]-2*psi1[(j+1)*numxcells-1]+psi1[j*numxcells-1])/(hy*hy));
      // Определение q
      betax[1]=-fx[0]/a;
      for(i=2;i<numxcells-1;i++)
      {
	betax[i]=(betax[i-1]-fx[i-1])/(alphax[i-1]+a);
      }
      q[numxcells-2]=(betax[numxcells-2]-fx[numxcells-2])/(alphax[numxcells-2]+a);
      for(i=numxcells-3;i>=0;i--)
      {
	q[i]=-alphax[i+1]*q[i+1]+betax[i+1];
      }
      psi2[(j+1)*numxcells-1]=(q[0]+q[numxcells-2]-fx[numxcells-1])/(a+p[0]+p[numxcells-2]);
      for(i=0;i<numxcells-1;i++)
      {
	psi2[j*numxcells+i]=-p[i]*psi2[(j+1)*numxcells-1]+q[i];
      }
    }

    a=2+hy*hy/tauy;
    alphay[1]=0;
    betay[1]=0;
    for(j=2;j<numycells;j++)
    {
      alphay[j]=-1/(alphay[j-1]+a);
    }
    for(i=0;i<numxcells;i++)
    {
      for(j=1;j<numycells-1;j++)
      {
	if(i==0)
	{
	  fy[j]=-hy*hy*(cells[j*numxcells].w+psi2[j*numxcells]/tauy+
	      (psi2[j*numxcells+1]-2*psi2[j*numxcells]+psi2[(j+1)*numxcells-1])/(hx*hx));
	}
	else if (i==numxcells-1)
	{
	  fy[j]=-hy*hy*(cells[(j+1)*numxcells-1].w+psi2[(j+1)*numxcells-1]/tauy+
	      (psi2[j*numxcells]-2*psi2[(j+1)*numxcells-1]+psi2[(j+1)*numxcells-2])/(hx*hx));
	}
	else
	{
	  fy[j]=-hy*hy*(cells[j*numxcells+i].w+psi2[j*numxcells+i]/tauy+
	      (psi2[j*numxcells+i+1]-2*psi2[j*numxcells+i]+psi2[j*numxcells+i-1])/(hx*hx));
	}
      }
      for(j=2;j<numycells;j++)
      {
	betay[j]=(betay[j-1]-fy[j-1])/(alphay[j-1]+a);
      }
      psi3[(numycells-1)*numxcells+i]=0;
      for(j=numycells-1;j>0;j--)
      {
	psi3[(j-1)*numxcells+i]=-alphay[j]*psi3[j*numxcells+i]+betay[j];
      }
    }
  }
  k=-1;
  for(j=0;j<numycells;j++)
    for(i=0;i<numxcells;i++)
    {
      k++;
      cells2[k].psi=psi3[k];
    }
}

void teascheme_init(void)
{
  Jw=(double*)malloc(numxcells*numycells*sizeof(double));
  if(Jw==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  Jwpp=(double*)malloc(numxcells*numycells*sizeof(double));
  if(Jwpp==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  Jwxp=(double*)malloc(numxcells*numycells*sizeof(double));
  if(Jwxp==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  Jwpx=(double*)malloc(numxcells*numycells*sizeof(double));
  if(Jwpx==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  Jtheta=(double*)malloc(numxcells*numycells*sizeof(double));
  if(Jtheta==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  Tw=(double*)malloc(numxcells*numycells*sizeof(double));
  if(Tw==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  Ttheta=(double*)malloc(numxcells*numycells*sizeof(double));
  if(Ttheta==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  Fw=(double*)malloc(numxcells*numycells*sizeof(double));
  if(Fw==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  Ftheta=(double*)malloc(numxcells*numycells*sizeof(double));
  if(Ftheta==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  Fw_prev=(double*)malloc(numxcells*numycells*sizeof(double));
  if(Fw_prev==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  Ftheta_prev=(double*)malloc(numxcells*numycells*sizeof(double));
  if(Ftheta_prev==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  psi1=(double*)malloc(numxcells*numycells*sizeof(double));
  if(psi1==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  psi2=(double*)malloc(numxcells*numycells*sizeof(double));
  if(psi2==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  psi3=(double*)malloc(numxcells*numycells*sizeof(double));
  if(psi3==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  fx=(double*)malloc(numxcells*sizeof(double));
  if(fx==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  fy=(double*)malloc(numycells*sizeof(double));
  if(fy==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  p=(double*)malloc(numxcells*sizeof(double));
  if(p==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  q=(double*)malloc(numxcells*sizeof(double));
  if(q==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  alphax=(double*)malloc(numxcells*sizeof(double));
  if(alphax==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  betax=(double*)malloc(numxcells*sizeof(double));
  if(betax==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  alphay=(double*)malloc(numycells*sizeof(double));
  if(alphay==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  betay=(double*)malloc(numycells*sizeof(double));
  if(betay==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }

  double lamb1,lamb2;
  double Lamb1,Lamb2;
  double epsilonpsi=1e-5;
  double a,nu;

  lamb1=4*sin(pi/numxcells)*sin(pi/numxcells)/(hx*hx); 
  Lamb1=4*sin(0.5*pi*(numxcells-1)/numxcells)*sin(0.5*pi*(numxcells-1)/numxcells)/(hx*hx);
  lamb2=4*sin(0.5*pi/numycells)*sin(0.5*pi/numycells)/(hy*hy); 
  Lamb2=4*sin(0.5*pi*(numycells-1)/numycells)*sin(0.5*pi*(numycells-1)/numycells)/(hy*hy);
  a=sqrt((Lamb1-lamb1)*(Lamb2-lamb2)/(Lamb1+lamb2)/(Lamb2+lamb1));
  nu=(1-a)/(1+a);

  printf("%f %f %f %f\n",Lamb1,lamb1,Lamb2,lamb2);
  
  double temp1;

//  temp1=(log(4/nu)*log(4/epsilonpsi)/pi/pi);

  temp1=-2*sqrt(((double)numxcells)*((double)numycells))/3.14*log(epsilonpsi);

//  temp1=-log(epsilonpsi)/(2*(sqrt(lamb1/Lamb1)+sqrt(lamb2/Lamb2)));

  M=(int)(temp1)+1;

  tauxtab=(double*)malloc(M*sizeof(double));
  if(tauxtab==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  tauytab=(double*)malloc(M*sizeof(double));
  if(tauytab==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }

  double *sigma,*xi;
  sigma=(double*)malloc(M*sizeof(double));
  if(sigma==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  xi=(double*)malloc(M*sizeof(double));
  if(xi==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }

  
  double q;
  double b,t,r,s;
  
  q=nu*nu*(1+nu*nu/2)/16;
  b=a*(Lamb2+lamb1)/(Lamb1-lamb1);
  t=(1-b)/(1+b);
  r=(Lamb2+b*Lamb1)/(1+b);
  s=(Lamb2-b*Lamb1)/(1+b);  

  int m;

  for(m=M/2;m<M;m++)
  {
    sigma[m]=(2*((double)m)+1)/(2*((double)M));   //  + is right!
    xi[m]=sqrt(nu)*exp(0.25*(2*sigma[m]-1)*log(q))*(1+exp((1-sigma[m])*log(q))+exp((1+sigma[m])*log(q)))/
      (1+exp(sigma[m]*log(q))+exp((2-sigma[m])*log(q)));
    tauxtab[m]=(1+t*xi[m])/(r*xi[m]+s);
    tauytab[m]=(1-t*xi[m])/(r*xi[m]-s);
  }
  for(m=0;m<M/2;m++)
  {
    xi[m]=nu/xi[M-1-m];
    tauxtab[m]=(1+t*xi[m])/(r*xi[m]+s);
    tauytab[m]=(1-t*xi[m])/(r*xi[m]-s);
  }

  for(m=0;m<M;m++)
  {
    tauxtab[m]=1/(Lamb1+lamb1);
    tauytab[m]=1/(Lamb2+lamb2);
  }

  FILE *out;
  out=fopen("tau.txt","w");
  fprintf(out,"M=%d\n",M);
  for(m=0;m<M;m++)
    fprintf(out,"taux=%f\ttauy=%f\n",tauxtab[m],tauytab[m]);
  for(m=0;m<M;m++)
    fprintf(out,"xi=%f\n",xi[m]);
  fprintf(out,"nu=%f\n",nu);
  fclose(out);

  free(sigma);
  free(xi);
}

void teascheme_finish(void)
{
  free(Jw);
  free(Jwpp);
  free(Jwxp);
  free(Jwpx);
  free(Jtheta);
  free(Tw);
  free(Ttheta);
  free(Fw);
  free(Ftheta);
  free(Fw_prev);
  free(Ftheta_prev);
  free(psi1);
  free(psi2);
  free(psi3);
  free(fx);
  free(fy);
  free(p);
  free(q);
  free(alphax);
  free(betax);
  free(alphay);
  free(betay);
  free(tauxtab);
  free(tauytab);
}
