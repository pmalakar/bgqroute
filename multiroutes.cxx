/*
 * Author: Preeti Malakar
 *
 * Test code for testing libroute
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "route.h"

int main (int argc, char **argv) {

	MPI_Init( &argc, &argv);

	int myrank, commsize, source, destination;
	
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
	MPI_Comm_size (MPI_COMM_WORLD, &commsize);

//  source = atoi(argv[1]);
//	destination = atoi(argv[2]);

	FILE *fp = fopen(argv[1], "r");
  while(fscanf(fp, "%d %d", &source, &destination) > 0)
		if (myrank == source)
		{
			printf ("Route from %d - %d\n", source, destination);
			puts("************");
		
			char path[64];
			getRoute (source, destination, path);
			puts(path);
	
			puts("************");
			fflush(stdout);
		}

	fclose(fp);

	MPI_Finalize();

	return 0;
}
