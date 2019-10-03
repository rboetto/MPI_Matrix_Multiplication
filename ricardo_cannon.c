#include "ricardo_cannon.h"

int fast_square_root(int n);
void shift_up (double * mat, int n, int horiz, int vert, double * new_row);
void shift_left (double * mat, int n, int horiz, int vert, double * new_col);

int main(int argc, char *argv[]) {

	MPI_Init(&argc, &argv);

	int id, p;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// Decompose communicator into greatest square < p
	int root_p = fast_square_root(p);
	int g_sqr = root_p * root_p;

	FILE * fpa = fopen(argv[1], "r");
	FILE * fpb = fopen(argv[2], "r");

	// Validate matrices, get size of N
	int n;
	if (!id) {
		int m, o, p;
		fread(&n, sizeof(int), 1, fpa);
		fread(&m, sizeof(int), 1, fpa);
		fread(&o, sizeof(int), 1, fpb);
		fread(&p, sizeof(int), 1, fpb);

		if ((m != n)||(o != p)) {
			perror("ERROR: Files must contain NxN matrices\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		if (m != o) {
			perror("ERROR: Matrices must have same dimensions\n");
			MPI_Abort(MPI_COMM_WORLD, 2);
		}
	} else {
		fseek(fpa, sizeof(int)*2, SEEK_SET);
		fseek(fpb, sizeof(int)*2, SEEK_SET);
	}

	MPI_Bcast(&n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
	fpos_t top_a, top_b;
	fgetpos(fpa, &top_a);
	fgetpos(fpb, &top_b);

	// Clamp # of processes to number of elements in matrices
	g_sqr = MIN(g_sqr,n*n);
	root_p = MIN(root_p,n);

	int sub_size = n*n/g_sqr;
	int z = fast_square_root(sub_size);
	if (z*z != sub_size) {
		perror("ERROR: Submatrix dimensions must be square numbers. View README for instructions.\n");
		MPI_Abort(MPI_COMM_WORLD, 3);
	}

	// Create mew communicator with usable processes
	MPI_Comm cart_comm, sub_comm;
	MPI_Comm_split(MPI_COMM_WORLD, id < g_sqr, id, &sub_comm);

	// Only processes within n / greatest square within num proc. do anything
	if (id < g_sqr) {
		int periods[2] = {root_p,root_p};
		MPI_Cart_create(sub_comm, 2, periods, periods, 0, &cart_comm);
	} else {
		MPI_Finalize();
		return 0;
	}

	// Determine Cartesian ID for process
	int coords[2];
	MPI_Cart_coords(cart_comm, id, 2, coords);
	int source, destination;
	MPI_Cart_shift(cart_comm, 1, root_p-1, &source, &destination);


	int horiz = BLOCK_SIZE(coords[0], root_p, n);
	int vert = BLOCK_SIZE(coords[1], root_p, n);
	// For communication reasons, each submatrix must be NxN as well

	double mat_a[z][z];
	double mat_b[z][z];
	double mat_c[z][z];

	int x_low = BLOCK_LOW(coords[0],root_p,n);
	int x_high = BLOCK_HIGH(coords[0],root_p,n);
	int y_low = BLOCK_LOW(coords[1],root_p,n);
	int y_high = BLOCK_HIGH(coords[1],root_p,n);


	// Read files
	// Note: Can be optimized by reading contiguous memory blocks as a single fread.
	int i, j, offset;
	int k = 0;
	for (i = x_low; i <= x_high; i++) {
		for (j = y_low; j <= y_high; j++) {
			offset = i+((n-j-1)*n);
			printf("Process %d will read the %dth (%d,%d) element\n", id, offset, i, n-j-1);
			fsetpos(fpa, &top_a);
			fseek(fpa, sizeof(double)*offset,SEEK_CUR);
			fread((mat_a+k), sizeof(double), 1, fpa);
			offset = ((i+j)%n)*n+j;
			fsetpos(fpb, &top_b);
			fseek(fpb, sizeof(double)*offset,SEEK_CUR);
			fread((mat_b+k), sizeof(double), 1, fpb);
			k++;
		}
	}

	MPI_Barrier(sub_comm);

	// Communication metadata
	int src, dst;
	MPI_Cart_shift(cart_comm, 0, -1, &src, &dst);

	// Initial skew of matrix A
	for (i = 0; i < vert; i++) {
		for (j = i; j > 0; j--) {

		}
	}
	/*
	int p_blw, p_abv, p_lft, p_rgt;
	MPI_Cart_shift(cart_comm, 0, 1, &p_lft, &p_rgt);
	MPI_Cart_shift(cart_comm, 1, 1, &p_blw, &p_abv);

	MPI_Status not_used;

	MPI_Barrier(sub_comm);
	printf("%i --> %i(%i/%i/%i) --> %i\n", p_blw, id, horiz, vert, area, p_abv);
	MPI_Barrier(sub_comm);

	double a, b;
	for (k = 1; k < n; k++) {
		for (i = 0; i <= area; i++) {
			*(mat_c+i) = *(mat_a+i) * *(mat_b+i);
		}

		for (i = 0; i <= horiz; i++) *(row_send+i) = *(mat_a+(a_top*n)+i);
		for (i = 0; i < vert; i++) *(col_send+i) = *(mat_b+(i*n)+b_left);

		if (!coords[1]) {
			MPI_Send(row_send, horiz, MPI_DOUBLE, p_abv, 0, sub_comm);
			MPI_Recv(row_recv, horiz, MPI_DOUBLE, p_blw, 0, sub_comm, &not_used);

		} else {
			MPI_Recv(row_recv, horiz, MPI_DOUBLE, p_blw, 0, sub_comm, &not_used);
			MPI_Send(row_send, horiz, MPI_DOUBLE, p_abv, 0, sub_comm);
		}

		if (!coords[0]) {
			MPI_Send(col_send, vert, MPI_DOUBLE, p_rgt, 0, sub_comm);
			MPI_Recv(col_recv, vert, MPI_DOUBLE, p_lft, 0, sub_comm, &not_used);
		} else {
			MPI_Recv(col_recv, vert, MPI_DOUBLE, p_rgt, 0, sub_comm, &not_used);
			MPI_Send(col_send, vert, MPI_DOUBLE, p_lft, 0, sub_comm);
		}

		printf("flag");
		MPI_Barrier(cart_comm);

		shift_up(mat_a, n, horiz, vert, row_recv);
		shift_left(mat_b, n, horiz, vert, col_recv);


	}
	*/
	// write_checkerboard_matrix(argv[3], (void**)&mat_c, MPI_DOUBLE, n, n, cart_comm);


	// Shift


	MPI_Finalize();
	return 0;
	// X = (i/(upper_bound-lower_bound+1))+lower_bound
	// Y = (i%(upper_bound-lower_bound))+lower_bound 

}

void shift_up (double * mat, int n, int horiz, int vert, double * new_row) {

	int i, j;
	for (i = 0; i < vert - 1; i++) {
		for (j = 0; j < horiz; j++) {
			*(mat + i*n + j) = *(mat + (i+1)*n + j);

		}
	}

	for (i = 0; i < horiz; i++) {
		*(mat + (vert-1)*n + i) = *(new_row+i);

	}


}

void shift_left (double * mat, int n, int horiz, int vert, double * new_col) {

	int i, j;
	for (i = 0; i < vert; i++) {
		for (j = 0; j < horiz - 1; j++) {
			*(mat + i*n + j) = *(mat + i*n + j+1);
		}
	}

	for (i = 0; i < vert; i++) {
		*(mat + i*n + horiz - 1) = *(new_col+i);
	}


}


// This algorithm is basically black magic
// https://en.wikipedia.org/wiki/Fast_inverse_square_root
int fast_square_root(int num) {
	const float x2 = num * 0.5f;
	const float threehalfs = 1.5f;

	union {
		float f;
		int i;
	} conv = {(float)num};

	conv.i = 0x5f3759df - ( conv.i >> 1);
	conv.f *= (threehalfs - (x2 * conv.f * conv.f));
	return (int)(1/conv.f);
}
