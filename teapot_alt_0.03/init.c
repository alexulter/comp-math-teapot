#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "types.h"
#include "extern.h"
#include "bounds.h"
#include "func.h"
#include "scheme.h"

void Init(void)
{
  printf("Init\n");
  T=0.0;
  dt0=0.1;
  hx=0.2;
  hy=0.2;
  K=1e-2;
  beta=4e-3;
  numxcells = 200;
  numycells = 50;
  cells=(struct Cell*)malloc(numxcells*numycells*sizeof(struct Cell));
  if(cells==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }
  cells2=(struct Cell*)malloc(numxcells*numycells*sizeof(struct Cell));
  if(cells2==NULL)
  {
    fprintf(stderr,"Can't allocate memory\n");
    exit(1);
  }

  int i,j;
  long k;
  k=-1;
  for(j=0;j<numycells;j++)
  {
    for(i=0;i<numxcells;i++)
    {
      k++;
      cells[k].w=0.0;
      cells[k].theta=0.0;
      cells[k].psi=0.0;
    }
  }
  
  teascheme_init();
  scheme=&teascheme;
  schemefinish=&teascheme_finish;
  bounds=&BoundTea;

  numevents=2;
  events=(struct Event*)malloc(numevents*sizeof(struct Event));
  events[0].type=1;
  events[0].T=5e2;
  events[1].type=2;
  events[1].param=5;
  events[1].T=events[1].param;
  calcnextevent();

  FILE *out;
  char filename[30];
  sprintf(filename,"out/time.dat");
  out=fopen(filename,"w");
  fclose(out);

}

void Finish(void)
{
  schemefinish();
  free(cells);
  free(cells2);
}
