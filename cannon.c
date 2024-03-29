/*
 * Code template for implementing ...
 *   matrix multiplication using Cannon's Algorithm in MPI
 *
 * The program takes three command-line arguments: fileA, fileB, and
 * fileC. The first two files contain matrix A and B as the input. The
 * third file is used to store the result matrix C as the output. The
 * program compute: C = A x B. The program assumes the matrices A, B,
 * and C are n x n matrices, the number of processors p is square, and
 * n is evenly divisible by sqrt(p).
 *
 * The files containing the matrices are all binary files and have the
 * following format. The matrix is stored in row-wise order and
 * preceded with two integers that specify the dimensions of the
 * matrix. The matrix elements are of integer types.
 *
 */

#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* set this parameter to reflect the cache line size of the particular
   machine you're running this program */
#define CACHE_SIZE 1024

/* in case later we decide to use another data type */
#define mpitype MPI_DOUBLE
typedef double datatype;

/* block decomposition macros */
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(j,p,n) (((p)*((j)+1)-1)/(n))

/* print out error message and exit the program */
void my_abort(const char* fmt, ...)
{
  int id;     /* process rank */
  va_list ap; /* argument list */

  va_start(ap, fmt);

  /* only process 0 reports */
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  if(!id) vprintf(fmt, ap);

  va_end(ap);

  /* all MPI processes exit at this point */
  exit(1);
}

/* return the data size in bytes */
int get_size(MPI_Datatype t)
{
   if(t == MPI_BYTE) return sizeof(char);
   else if(t == MPI_DOUBLE) return sizeof(double);
   else if(t == MPI_FLOAT) return sizeof(float);
   else if(t == MPI_INT) return sizeof(int);
   else {
     printf ("Error: Unrecognized argument to 'get_size'\n");
     fflush (stdout);
     MPI_Abort (MPI_COMM_WORLD, -3);
   }
   return 0;
}


/* allocate memory from heap */
void *my_malloc(int id, int bytes)
{
   void *buffer;
   if ((buffer = malloc ((size_t) bytes)) == NULL) {
      printf ("Error: Malloc failed for process %d\n", id);
      fflush (stdout);
      MPI_Abort (MPI_COMM_WORLD, -2);
   }
   return buffer;
}

/* Read a matrix from a file. */
void read_checkerboard_matrix (
   char *s,              /* IN - File name */
   void ***subs,         /* OUT - 2D array */
   void **storage,       /* OUT - Array elements */
   MPI_Datatype dtype,   /* IN - Element type */
   int *m,               /* OUT - Array rows */
   int *n,               /* OUT - Array cols */
   MPI_Comm grid_comm)   /* IN - Communicator */
{
   void      *buffer;         /* File buffer */
   int        coords[2];      /* Coords of proc receiving
                                 next row of matrix */
   int        datum_size;     /* Bytes per elements */
   int        dest_id;        /* Rank of receiving proc */
   int        grid_coord[2];  /* Process coords */
   int        grid_id;        /* Process rank */
   int        grid_period[2]; /* Wraparound */
   int        grid_size[2];   /* Dimensions of grid */
   int        i, j, k;
   FILE      *infileptr;      /* Input file pointer */
   void      *laddr;          /* Used when proc 0 gets row */
   int        local_cols;     /* Matrix cols on this proc */
   int        local_rows;     /* Matrix rows on this proc */
   void     **lptr;           /* Pointer into 'subs' */
   int        p;              /* Number of processes */
   void      *raddr;          /* Address of first element
                                 to send */
   void      *rptr;           /* Pointer into 'storage' */
   MPI_Status status;         /* Results of read */

   MPI_Comm_rank (grid_comm, &grid_id);
   MPI_Comm_size (grid_comm, &p);
   datum_size = get_size (dtype);

   /* Process 0 opens file, gets number of rows and
      number of cols, and broadcasts this information
      to the other processes. */

   if (grid_id == 0) {
      infileptr = fopen (s, "r");
      if (infileptr == NULL) *m = 0;
      else {
         fread (m, sizeof(int), 1, infileptr);
         fread (n, sizeof(int), 1, infileptr);
      }
   }
   MPI_Bcast (m, 1, MPI_INT, 0, grid_comm);

   if (!(*m)) MPI_Abort (MPI_COMM_WORLD, -1);

   MPI_Bcast (n, 1, MPI_INT, 0, grid_comm);

   /* Each process determines the size of the submatrix
      it is responsible for. */

   MPI_Cart_get (grid_comm, 2, grid_size, grid_period,
      grid_coord);
   local_rows = BLOCK_SIZE(grid_coord[0],grid_size[0],*m);
   local_cols = BLOCK_SIZE(grid_coord[1],grid_size[1],*n);

   /* Dynamically allocate two-dimensional matrix 'subs' */

   *storage = my_malloc (grid_id,
      local_rows * local_cols * datum_size);
   *subs = (void **) my_malloc (grid_id,local_rows*sizeof(void*));
   lptr = (void *) *subs;
   rptr = (void *) *storage;
   for (i = 0; i < local_rows; i++) {
      *(lptr++) = (void *) rptr;
      rptr += local_cols * datum_size;
   }

   /* Grid process 0 reads in the matrix one row at a time
      and distributes each row among the MPI processes. */

   if (grid_id == 0)
      buffer = my_malloc (grid_id, *n * datum_size);

   /* For each row of processes in the process grid... */
   for (i = 0; i < grid_size[0]; i++) {
      coords[0] = i;

      /* For each matrix row controlled by this proc row...*/
      for (j = 0; j < BLOCK_SIZE(i,grid_size[0],*m); j++) {

         /* Read in a row of the matrix */

         if (grid_id == 0) {
            fread (buffer, datum_size, *n, infileptr);
         }

         /* Distribute it among process in the grid row */

         for (k = 0; k < grid_size[1]; k++) {
            coords[1] = k;

            /* Find address of first element to send */
            raddr = buffer +
               BLOCK_LOW(k,grid_size[1],*n) * datum_size;

            /* Determine the grid ID of the process getting
               the subrow */
            MPI_Cart_rank (grid_comm, coords, &dest_id);

            /* Process 0 is responsible for sending...*/
            if (grid_id == 0) {

               /* It is sending (copying) to itself */
               if (dest_id == 0) {
                  laddr = (*subs)[j];
                  memcpy (laddr, raddr,
                     local_cols * datum_size);

               /* It is sending to another process */
               } else {
                  MPI_Send (raddr,
                     BLOCK_SIZE(k,grid_size[1],*n), dtype,
                  dest_id, 0, grid_comm);
               }

            /* Process 'dest_id' is responsible for
               receiving... */
            } else if (grid_id == dest_id) {
               MPI_Recv ((*subs)[j], local_cols, dtype, 0,
                  0, grid_comm,&status);
            }
         }
      }
   }
   if (grid_id == 0) free (buffer);
}

/*
 * Write a matrix distributed in checkerboard fashion to a file.
 */
void write_checkerboard_matrix (
   char        *s,	      /* IN -File name */
   void       **a,            /* IN -2D matrix */
   MPI_Datatype dtype,        /* IN -Matrix element type */
   int          m,            /* IN -Matrix rows */
   int          n,            /* IN -Matrix columns */
   MPI_Comm     grid_comm)    /* IN -Communicator */
{
   void      *buffer;         /* Room to hold 1 matrix row */
   int        coords[2];      /* Grid coords of process
                                 sending elements */
   int        datum_size;     /* Bytes per matrix element */
   int        els;            /* Elements received */
   int        grid_coords[2]; /* Coords of this process */
   int        grid_id;        /* Process rank in grid */
   int        grid_period[2]; /* Wraparound */
   int        grid_size[2];   /* Dims of process grid */
   int        i, j, k;
   void      *laddr;          /* Where to put subrow */
   int        local_cols;     /* Matrix cols on this proc */
   int        p;              /* Number of processes */
   int        src;            /* ID of proc with subrow */
   MPI_Status status;         /* Result of receive */
   FILE      *outfileptr;     /* Output file */

   MPI_Comm_rank (grid_comm, &grid_id);
   MPI_Comm_size (grid_comm, &p);
   datum_size = get_size (dtype);

   if(grid_id == 0) {
     outfileptr = fopen (s, "w");
     if(!outfileptr ||
	fwrite(&m, sizeof(int), 1, outfileptr) != 1 ||
	fwrite(&n, sizeof(int), 1, outfileptr) != 1)
       MPI_Abort (MPI_COMM_WORLD, -1);
   }

   MPI_Cart_get (grid_comm, 2, grid_size, grid_period,
      grid_coords);
   local_cols = BLOCK_SIZE(grid_coords[1],grid_size[1],n);

   if (!grid_id)
      buffer = my_malloc (grid_id, n*datum_size);

   /* For each row of the process grid */
   for (i = 0; i < grid_size[0]; i++) {
      coords[0] = i;

      /* For each matrix row controlled by the process row */
      for (j = 0; j < BLOCK_SIZE(i,grid_size[0],m); j++) {

         /* Collect the matrix row on grid process 0 and
            print it */
         if (!grid_id) {
            for (k = 0; k < grid_size[1]; k++) {
               coords[1] = k;
               MPI_Cart_rank (grid_comm, coords, &src);
               els = BLOCK_SIZE(k,grid_size[1],n);
               laddr = buffer +
                  BLOCK_LOW(k,grid_size[1],n) * datum_size;
               if (src == 0) {
                  memcpy (laddr, a[j], els * datum_size);
               } else {
                  MPI_Recv(laddr, els, dtype, src, 0,
                     grid_comm, &status);
               }
            }

	    if(fwrite(buffer, datum_size, n, outfileptr) != n)
	      MPI_Abort (MPI_COMM_WORLD, -1);
         } else if (grid_coords[0] == i) {
            MPI_Send (a[j], local_cols, dtype, 0, 0,
               grid_comm);
         }
      }
   }
   if (!grid_id) {
      free (buffer);
      fclose(outfileptr);
   }
}

/* recursive, block-oriented, sequential matrix multiplication */
void my_matmul(int crow, int ccol, /* corner of C block */
	       int arow, int acol, /* corner of A block */
	       int brow, int bcol, /* corner of B block */
	       int l, int m, int n, /* block A is l*m, block B is m*n, block C is l*n */
	       int N, /* matrices are N*N */
	       datatype** a, datatype** b, datatype** c) { /* 2D matrices */
  int i, j, k; 
  int lhalf[3], mhalf[3], nhalf[3]; /* quadrant sizes */
  datatype *aptr, *bptr, *cptr;

  if(m*n*sizeof(datatype) > CACHE_SIZE) { /* block B doesn't fit in cache */
    lhalf[0] = 0; lhalf[1] = l/2; lhalf[2] = l-l/2;
    mhalf[0] = 0; mhalf[1] = m/2; mhalf[2] = m-m/2;
    nhalf[0] = 0; nhalf[1] = n/2; nhalf[2] = n-n/2;
    for(i=0; i<2; i++)
      for(j=0; j<2; j++)
	for(k=0; k<2; k++)
	  my_matmul(crow+lhalf[i], ccol+nhalf[j],
		    arow+lhalf[i], acol+mhalf[k],
		    brow+mhalf[k], bcol+nhalf[j],
		    lhalf[i+1], mhalf[k+1], nhalf[j+1],
		    N, a, b, c);
  } else { /* block B fits in cache */
    for(i=0; i<l; i++) {
      for(j=0; j<n; j++) {
	cptr = &c[crow+i][ccol+j];
	aptr = &a[arow+i][acol];
	bptr = &b[brow][bcol+j];
	for(k=0; k<m; k++) {
	  *cptr += *(aptr++) * (*bptr);
	  bptr += N;
	}
      }
    }
  }
}
