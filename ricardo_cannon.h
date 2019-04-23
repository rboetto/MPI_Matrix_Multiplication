#ifndef RICARDO_CANNON
#define RICARDO_CANNON
#include <stdio.h>
#include <mpi.h>

#define BLOCK_LOW(i,p,n) ((i)*(n)/(p))
#define BLOCK_HIGH(i,p,n) (BLOCK_LOW((i)+1,p,n)-1)
#define BLOCK_SIZE(i,p,n) (BLOCK_LOW((i)+1,p,n)-BLOCK_LOW(i,p,n))


#endif
