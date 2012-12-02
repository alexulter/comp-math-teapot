#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "types.h"
#include "extern.h"
#include "init.h"
#include "shared.h"

int filenom=0;

double min(double x, double y)
{
  if(x<y) 
    return x;
  return y;
}

double absd(double x)
{
  if(x<0.0)
    return -x;
  return x;
}

void save(void)
{
  printf("Saving %.4d T=%e\n",filenom,T);
  char filename[30];
  FILE *out;

  sprintf(filename,"out/time.dat");
  out=fopen(filename,"a");
  fprintf(out,"%e\n",T);
  fclose(out);

  sprintf(filename,"out/w%.4d.dat",filenom);
  out=fopen(filename,"w");
  int i,j,k;
  k=-1;
  for(j=0;j<numycells;j++)
  {
    for(i=0;i<numxcells;i++)
    {
      k++;
      fprintf(out,"%e ",cells[k].w);
    }
    fprintf(out,"\n");
  }
  fclose(out);
  
  sprintf(filename,"out/theta%.4d.dat",filenom);
  out=fopen(filename,"w");
  k=-1;
  for(j=0;j<numycells;j++)
  {
    for(i=0;i<numxcells;i++)
    {
      k++;
      fprintf(out,"%e ",cells[k].theta);
    }
    fprintf(out,"\n");
  }
  fclose(out);
  
  sprintf(filename,"out/psi%.4d.dat",filenom);
  out=fopen(filename,"w");
  k=-1;
  for(j=0;j<numycells;j++)
  {
    for(i=0;i<numxcells;i++)
    {
      k++;
      fprintf(out,"%e ",cells[k].psi);
    }
    fprintf(out,"\n");
  }
  fclose(out);

  filenom++;
}

void copycells(struct field** space_array)
{
  int i,j,k;
  k=-1;
  for(j=0;j<numycells;j++)
  {
    for(i=0;i<numxcells;i++)
    {
      k++;
      if((j>0) && (j<numycells-1))
      {
	cells[k].w=cells2[k].w;
	cells[k].theta=cells2[k].theta;
	cells[k].psi=cells2[k].psi;
	
	space_array[i+1][j+1].psi = cells2[k].psi*1e+10;
	
	space_array[i+1][j+1].omega = cells2[k].w*1e+8;
	space_array[i+1][j+1].teta = cells2[k].theta;

      }
    }
  }
printf("%e \n", space_array[150][20].psi);
}

void calcdt(void)
{
  dt=dt0;
}

void calcnextevent(void)
{
  int i;
  nextevent=events[0].T;  // Исправить для отработанных событий
  for(i=1;i<numevents;i++)
  {
    if(T<events[i].T)
      nextevent=min(nextevent,events[i].T);
  }
}

void doevents(void)
{
  int i;
  double epsilon=1e-2*dt;
  for(i=0;i<numevents;i++)
  {
    if(events[i].type==2)
    {
      if(absd(T-events[i].T)<epsilon)
      {
	save();
	events[i].T=events[i].T+events[i].param;
      }
    }
    if(events[i].type==1)
    {
      if(absd(T-events[i].T)<epsilon)
      {
	save();
	Finish();
	exit(0);
      }
    }
  }
  calcnextevent();
}

/*void save_calc(void)
{
  char filename[30];
  FILE *out;
  out=fopen("inf.dat","w");
  fprintf(out,"%d %d\n",numxcells, numycells);
  fprintf(out,"%e %e\n",T, dt0);
  fclose(out);

  out=fopen("data.dat","w");
  int i,j;
  long k;
  k=-1;
  for(j=0;j<numycells;j++)
    for(i=0;i<numxcells;i++)
    {
      k++;
      fprintf(out,"%e %e %e %e",cells[k].rho,cells[k].x,cells[k].y,cells[k].vx,cells[k].vy,
	  cells[k].p,cells[k].e);
    }
  fclose(out);
}*/

