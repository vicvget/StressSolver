// Submodel.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <stdlib.h>
#include "time.h"

#include "rezrsub.h"

int main()
{

	Rezrsub *my;
	if (my->testCUDA())
		my = new RezrsubGPU();
	else
		my = new RezrsubCPU();

	delete my;

	return 0;
}

