/*
 *  Created by Preeti Malakar
 *  Argonne National Laboratory
 *
 *  For BG/Q: Determine routing order, and find the route
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <mpi.h>
#include <mpix.h>

#include "route.h" 

/*
 * Read the DCR registers to get the routing order 
 * The routing order depends on the partition in which the job runs
 * [Thanks to Phil Heidelberger]
 *
 * getRoutingOrder() populates parameter ro with the routing order, viz.  
 * A=0 B=1 C=2 D=3 E=4
 *
 * Example, for 512 nodes on Mira, ro = {0, 1, 2, 3, 4}
 *
 */

/* 
 * Structure containing BGQ parameters 
 */
MPIX_Hardware_t hw; 

/*
 * Ranks per node, Core ID (0...15) 
 */ 
int ppn, coreID, nodeID;

/*
 * Size of each dimension 
 */ 
int dimSize[MPIX_TORUS_MAX_DIMS];	// torus dimension size

/*
 * Whether dimension is torus or not	//TODO bool 
 */ 
int isTorus[MPIX_TORUS_MAX_DIMS];	// torus wrap = 0/1 

/*
 *  Routing order : varies based on the partition, node count etc. 
 */
int *routingOrder;

inline int min (int x, int y) {
  return x<y? x: y;
}

int getRoutingOrder (int *ro) {

	uint64_t dcr_det_order = DCRReadUser(ND_500_DCR(CTRL_DET_ORDER));

	int A = ND_500_DCR__CTRL_DET_ORDER__MASK0_get(dcr_det_order);
	int B = ND_500_DCR__CTRL_DET_ORDER__MASK1_get(dcr_det_order);
	int C = ND_500_DCR__CTRL_DET_ORDER__MASK2_get(dcr_det_order);
	int D = ND_500_DCR__CTRL_DET_ORDER__MASK3_get(dcr_det_order);
	int E = ND_500_DCR__CTRL_DET_ORDER__MASK4_get(dcr_det_order);

	int torus[5] = {A, B, C, D, E};

	int index = 0;
  for (int i=0; i<5 ; i++) { 
	 	if (torus[i] == 1) 					index = 4;
	 	else if (torus[i] == 2)			index = 3;
	 	else if (torus[i] == 4)			index = 2;
	 	else if (torus[i] == 8)			index = 1;
	 	else if (torus[i] == 16)		index = 0;
    ro[i] = index;
	}

	return 0;
}

void getRoute(int srcRank, int destRank, char *path) {

	//local variables
  int i;
	int unitHop = 1;
  int srcCoords[6], destCoords[6];
	int childCoords[6], intmdtCoords[6];
	char buf[64];

	routingOrder = new int[MPIX_TORUS_MAX_DIMS];
	getRoutingOrder(routingOrder);

	MPIX_Hardware(&hw);
	ppn = hw.ppn;
  coreID = hw.coreID;	
  nodeID = srcRank/ppn;

//	printf("MPIX_TORUS_MAX_DIMS = %d\n", MPIX_TORUS_MAX_DIMS);
	for (i=0; i<MPIX_TORUS_MAX_DIMS ; i++) {
		isTorus[i] = hw.isTorus[i];
		dimSize[i] = hw.Size[i];
	}

	MPIX_Rank2torus (srcRank, srcCoords);
	MPIX_Rank2torus (destRank, destCoords);

#ifdef DEBUG

	printf("Torus dimensions = (%u,%u,%u,%u,%u) Routing order = (%d,%d,%d,%d,%d)\n", hw.Size[0], hw.Size[1], hw.Size[2], hw.Size[3], hw.Size[4], routingOrder[0], routingOrder[1], routingOrder[2], routingOrder[3], routingOrder[4]);

	printf("Torus wraps? %u,%u,%u,%u,%u\n", hw.isTorus[0], hw.isTorus[1], hw.isTorus[2], hw.isTorus[3], hw.isTorus[4]);

	printf("Source rank: %d Node: %d Torus coords = (%u,%u,%u,%u,%u,%u) link ID: %d\n", srcRank, srcRank/ppn, srcCoords[0], srcCoords[1], srcCoords[2], srcCoords[3], srcCoords[4], srcCoords[5], destRank);

#endif

	//Initialize intermediate nodes in original path to the destination node
	for (int dim=0; dim < MPIX_TORUS_MAX_DIMS; dim++) 
		intmdtCoords[dim] = srcCoords[dim];
	   
	intmdtCoords[MPIX_TORUS_MAX_DIMS] = destCoords[MPIX_TORUS_MAX_DIMS];	//T

	int hopnum = 0;
	int hopDiff, intmdt_rank, child, parent;

	child = srcRank;

	for (int dim=0; dim<MPIX_TORUS_MAX_DIMS; dim++) {

		int dimID = routingOrder[dim];
		hopDiff = abs(destCoords[dimID] - srcCoords[dimID]);

		if (hw.isTorus[dimID] == 1) 
			hopDiff = min (hopDiff, hw.Size[dimID] - hopDiff) ;

#ifdef DEBUG
		printf("%d to %d difference in dim %d = %d\n", srcRank, destRank, dimID, hopDiff);
#endif

		for(int diff=0; diff<hopDiff ;diff++) {
			if (hw.isTorus[dimID] == 0) {
				if(destCoords[dimID] < srcCoords[dimID]) 
					intmdtCoords[dimID] -= unitHop;  
				else intmdtCoords[dimID] += unitHop;
			}
			else {		// torus
				if (abs(destCoords[dimID] - srcCoords[dimID])*2 > hw.Size[dimID]) 
				{
					if (destCoords[dimID] > srcCoords[dimID]) 
						intmdtCoords[dimID] = ((intmdtCoords[dimID] - unitHop) + hw.Size[dimID]) % hw.Size[dimID];  
					else 
						intmdtCoords[dimID] = (intmdtCoords[dimID] + unitHop) % hw.Size[dimID];
				}
				else if (abs(destCoords[dimID] - srcCoords[dimID])*2 < hw.Size[dimID]) 
				{
					if (destCoords[dimID] < srcCoords[dimID]) 
						intmdtCoords[dimID] = ((intmdtCoords[dimID] - unitHop) + hw.Size[dimID]) % hw.Size[dimID];  
					else 
						intmdtCoords[dimID] = (intmdtCoords[dimID] + unitHop) % hw.Size[dimID];
				}
				else 
				{
					//if source coord is even, plus direction
					if (srcCoords[dimID]%2 == 0)
						intmdtCoords[dimID] = (intmdtCoords[dimID] + unitHop) % hw.Size[dimID]; 		
					//even source coord: traverse in plus direction  
					else 
						intmdtCoords[dimID] = ((intmdtCoords[dimID] - unitHop) + hw.Size[dimID]) % hw.Size[dimID];  
				}
			}

			++hopnum;

			//get the rank
			MPIX_Torus2rank (intmdtCoords, &intmdt_rank);
			parent = intmdt_rank;

#ifdef DEBUG
			printf ("Enroute %d (%d %d %d %d %d %d) to %d (%d %d %d %d %d %d) Hop %d: in dimension %d Child %d to Parent %d (%d %d %d %d %d %d)\n", \
			srcRank, srcCoords[0], srcCoords[1], srcCoords[2], srcCoords[3], srcCoords[4], srcCoords[5], \
			destRank, destCoords[0], destCoords[1], destCoords[2], destCoords[3], destCoords[4], destCoords[5], \
			hopnum, dimID, child, intmdt_rank, intmdtCoords[0], intmdtCoords[1], intmdtCoords[2], intmdtCoords[3], intmdtCoords[4], intmdtCoords[5]);
#endif

	//		sprintf(buf, "%d ", child);
//			strcat(path, buf); 

#ifdef DEBUG
			//printf ("Route %d to %d Hop %d\n", srcRank, intmdt_rank, hopnum);
			printf ("Route %d to %d Hop %d\n", child, parent, hopnum);
			printf ("%d->%d;\n", hopnum, child, parent);
#endif
			MPIX_Rank2torus (child, childCoords);

			printf ("Hop %d: [%d-%d] %d (%d %d %d %d %d %d) -> %d (%d %d %d %d %d %d)\n", hopnum, srcRank, destRank, child, childCoords[0], childCoords[1], childCoords[2], childCoords[3], childCoords[4], childCoords[5], parent, intmdtCoords[0], intmdtCoords[1], intmdtCoords[2], intmdtCoords[3], intmdtCoords[4], intmdtCoords[5]); 

			child = parent;
		}
	}   

//	sprintf(buf, "%d", destRank);
//	strcat(path, buf); 


//#ifdef DEBUG
//	printf ("path: %s\n", path);
//#endif

	return;
			
}

