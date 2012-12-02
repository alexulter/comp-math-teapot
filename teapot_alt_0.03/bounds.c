#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "types.h"
#include "extern.h"

void BoundTea(void)
{
  int i;

  srand( (unsigned)time(NULL) );
  int rndnom;
//  double gausrnd;

  double theta0=10.0;
  double dtheta=1.0;

  double curtheta;

  for(i=0;i<numxcells;i++)
  {
    cells[i].w=0;
    cells[(numycells-1)*numxcells+i].w=0;
    cells[i].psi=0;
    cells[(numycells-1)*numxcells+i].psi=0;

//    Old version
    rndnom=rand()%2000-1000;
    curtheta=theta0+dtheta*((double)rndnom)/1000;

/*    rndnom=0;
    for(j=0;j<5;j++)
      rndnom=rndnom+rand()%1000;

    gausrnd=((double)rndnom)/5000;

    curtheta=cells[i].theta+dtheta*gausrnd;*/

    cells[i].theta=curtheta;
    cells[(numycells-1)*numxcells+i].theta=0;
  }

/*  for(i=numxcells/2-5;i<(numxcells/2)+5;i++)
  {
    cells[i].theta=10.0;
  }*/
}

