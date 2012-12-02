#include <stdio.h>
#include "types.h"
#include "extern.h"
#include "func.h"
#include "shared.h"

void Calculate(struct field **space_array)
{
  save();
  //step=0;
  //for(;;)
  //{
    //step++;
    dt = dt0;
    //printf("step=%ld T=%f dt=%e\n",step,T,dt);
	printf ("Calculation...\n");
    scheme();
    T=T+dt;
    copycells(space_array);
    bounds();
    doevents();
    //save();
  //}
}
