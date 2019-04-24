#ifndef RICARDO_CANNON
#define RICARDO_CANNON
#include <stdio.h>
#include <mpi.h>

#ifndef BLOCK_LOW
#define BLOCK_LOW(i,p,n) ((i)*(n)/(p))
#endif
#ifndef BLOCK_HIGH
#define BLOCK_HIGH(i,p,n) (BLOCK_LOW((i)+1,p,n)-1)
#endif
#ifndef BLOCK_SIZE
#define BLOCK_SIZE(i,p,n) (BLOCK_LOW((i)+1,p,n)-BLOCK_LOW(i,p,n))
#endif
#define MIN(a,b) ((a)<(b)?(a):(b))

#endif
