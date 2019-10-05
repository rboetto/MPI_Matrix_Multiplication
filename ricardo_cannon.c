#include "ricardo_cannon.h"

int fast_square_root (int n);
void shift_up (double * mat, int n, int horiz, int vert, double * new_row);
void shift_left (double * mat, int n, int horiz, int vert, double * new_col);
void shiftRow (double send[], int n, int amt, MPI_Comm comm);

int main(int argc, char *argv[]) {

	MPI_Init (&argc, &argv);

	int id, p;
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &p);

	FILE * fpa = fopen (argv[1], "r");
	FILE * fpb = fopen (argv[2], "r");

	if (!id) printf("id=%d p=%d\n", id, p);

	// Validate matrices, get size of N
	int n;
	if (!id) {
		int m, o, v;
		fread(&n, sizeof(int), 1, fpa);
		fread(&m, sizeof(int), 1, fpa);
		fread(&o, sizeof(int), 1, fpb);
		fread(&v, sizeof(int), 1, fpb);

		if ((m != n)||(o != v)) {
			perror ("ERROR: Files must contain NxN matrices\n");
			MPI_Abort (MPI_COMM_WORLD, 1);
		}
		if (m != o) {
			perror("ERROR: Matrices must have same dimensions\n");
			MPI_Abort(MPI_COMM_WORLD, 2);
		}
		if (p < n*n) {
			printf("%d %d\n", n*n, p);
			perror("ERROR: # of elements > # of processes\n");
			MPI_Abort(MPI_COMM_WORLD, 3);
		}
	} else {
		fread(&n, sizeof(int), 1, fpa);
		fseek(fpa, sizeof(int)*2, SEEK_SET);
		fseek(fpb, sizeof(int)*2, SEEK_SET);
	}

	// Get position of matrix beginning (no longer used)
	// fpos_t top_a, top_b;
	// fgetpos (fpa, &top_a);
	// fgetpos (fpb, &top_b);

	// Create separate communicator with only number of needed processes
	MPI_Comm cart_comm, sub_comm;
	MPI_Comm_split(MPI_COMM_WORLD, id < n*n, id, &sub_comm);
	p = n*n;

	// Create cartesian communicator for active processes, finalize other processes
	if (id < n*n) {
		int periods[2] = {n,n};
		MPI_Cart_create(sub_comm, 2, periods, periods, 0, &cart_comm);
	} else {
		MPI_Finalize();
		return 0;
	}

	// Determine Cartesian ID for process
	int coords[2];
	MPI_Cart_coords(cart_comm, id, 2, coords);

	// Calc. math coordinates from process coordinates
	int x = n-coords[1]-1;
	int y = coords[0];

	// Get init. numbers
	// Matrix alignment integrated
	double a, b, c;
	int offset = x*n+((y+x)%n);
	fseek(fpa, sizeof(double)*offset, SEEK_CUR);
	fread(&a, sizeof(double), 1, fpa);
	offset = ((x+y)%n)*n+y;
	fseek(fpb, sizeof(double)*offset, SEEK_CUR);
	fread(&b, sizeof(double), 1, fpb);
	c = 0;

	fclose(fpa);
	fclose(fpb);

	printf("%d: a=%f b=%f\n", id, a, b);

	// Determine adjacent processes
	int lft, rgt, top, btm;
	MPI_Cart_shift(cart_comm, 0, -1, &rgt, &lft);
	MPI_Cart_shift(cart_comm, 1, 1, &btm, &top);

	// Execute cannon's algorithm
	// Receive a from right, send to left
	// Receive b from bottom, send to top
	double a2, b2;
	MPI_Status status;

	int i;
	for (i = 0; i < n; i++) {
		c += a*b;
		if (!coords[0]) {
			MPI_Send(&a, 1, MPI_DOUBLE, lft, 0, cart_comm);
			MPI_Recv(&a2, 1, MPI_DOUBLE, rgt, 0, cart_comm, &status);
		} else {
			MPI_Recv(&a2, 1, MPI_DOUBLE, rgt, 0, cart_comm, &status);
			MPI_Send(&a, 1, MPI_DOUBLE, lft, 0, cart_comm);
		}
		if (!coords[1]) {
			MPI_Send(&b, 1, MPI_DOUBLE, top, 0, cart_comm);
			MPI_Recv(&b2, 1, MPI_DOUBLE, btm, 0, cart_comm, &status);
		} else {
			MPI_Recv(&b2, 1, MPI_DOUBLE, btm, 0, cart_comm, &status);
			MPI_Send(&b, 1, MPI_DOUBLE, top, 0, cart_comm);
		}
		a = a2;
		b = b2;
	}

	MPI_Barrier(cart_comm);
	printf("Result for cell (%d,%d) at position %d= %f\n", x, y, x*n+y, c);

	MPI_Barrier(cart_comm);


	int ord[p];
	double bucket[p];
	offset = x*n+y;
	MPI_Gather(&offset,1,MPI_INTEGER,ord,1,MPI_INTEGER,0,cart_comm);
	MPI_Gather(&c,1,MPI_DOUBLE,bucket,1,MPI_DOUBLE,0,cart_comm);

	if (!id) {
		FILE * fpc = fopen(argv[3],"w");
		fwrite(&n, sizeof(int), 1, fpc);
		fwrite(&n, sizeof(int), 1, fpc);
		int i;
		double result[p];
		for (i = 0; i < p; i++) {
			result[ord[i]] = bucket[i];
		}
		fwrite(result,sizeof(double),p,fpc);
		fclose(fpc);
	}


MPI_Finalize();
return 0;

}


// Everything below this line is my original code attempting to
// have each process operate on a variable number of elements
// before I took the 1-element-per-process approach

	/*
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

	// Initial skew of matrix A
	double buff[z];
	for (i = 0; i < z; i++) {
		for (j = 0; j < z; j++) {
			buff[j] = mat_a[j][i];
		}
		shiftRow(buff, z, i, cart_comm);
		for (j = 0; j < z; j++) {
		printf("%d: %f --> %f\n", id, mat_a[j][i], buff[j]);
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

void shiftRow (double send[], int n, int amt, MPI_Comm comm) {

int src, dst;
double recv[n];

// Assume a perfect displacement
int disp = (amt/n)+(amt%n!=0);
MPI_Cart_shift(comm, 0, -disp, &src, &dst);
int id;
MPI_Comm_rank(comm, &id);
// printf("%d: Receives from %d and sends to %d\n", id, src, dst);


MPI_Status status;

// Something is wrong with my usage of send/recv
MPI_Sendrecv(send, n, MPI_DOUBLE, dst, 0,
	recv, n, MPI_DOUBLE, src, 0, comm, &status);

int k;
for (k = 0; k < n; k++) {
	printf("%d: %f -> %f\n", id, send[k], recv[k]);
}

// If displacement fits evenly, wrap up
// Else, correct
if (amt%n==0) {
	int i;
	for (i = 0; i < n; i++) {
		send[i]=recv[i];
	}
} else {
	int d = amt%n;
	double send_p[d], recv_p[d];
	int i;
	for (i = 0; i < d; i++) send_p[i]=recv[d+i];
	MPI_Cart_shift(comm, 0, 1, &src, &dst);
	MPI_Sendrecv(send_p, d, MPI_DOUBLE, dst, 0, recv_p, d, MPI_DOUBLE, src, 0, comm, &status);
	for (i = 0; i < n; i++) {
		send[i]= i<d? recv_p[i] : recv[i-d];
	}
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
*/
