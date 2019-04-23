#include "ricardo_cannon.h"

int fast_square_root(int n);
int main_cont(MPI_Comm comm, int argc, char ** argv, int root);

int main(int argc, char * argv[]) {

	MPI_Init(&argc, &argv);

	int id, p;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// Decompose communicator into greatest square < p
	int d_sz = fast_square_root(p);
	int g_sqr = d_sz * d_sz;

	MPI_Comm sqr_comm;
	MPI_Comm_split(MPI_COMM_WORLD, id < g_sqr, id, &sqr_comm);

	// Only processes within square do anything
	if (id < g_sqr) {
		MPI_Comm cart_comm;
		int periods[2] = {d_sz,d_sz};
		MPI_Cart_create(sqr_comm, 2, periods, periods, 0, &cart_comm);
		main_cont(cart_comm, argc, argv, d_sz);
	}

	MPI_Finalize();
	return 0;
}

int main_cont(MPI_Comm comm, int argc, char * argv[], int root) {

	int id, p;
	MPI_Comm_rank(comm, &id);
	MPI_Comm_size(comm, &p);

	int coords[2];
	MPI_Cart_coords(comm, id, 2, coords);

	FILE * fpa = fopen(argv[1], "r");
	FILE * fpb = fopen(argv[2], "r");

	int n;

	if (!id) {
		int m, o, p;
		fread(&m, sizeof(int), 1, fpa);
		fread(&n, sizeof(int), 1, fpa);
		fread(&o, sizeof(int), 1, fpb);
		fread(&p, sizeof(int), 1, fpb);
		if ((m != n)||(o != p)) {
			perror("ERROR: Input files must be NxN matrices");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		if (m != o) {
			perror("Error: Input files must be matrices with same dimensions");
			MPI_Abort(MPI_COMM_WORLD, 2);
		}
	} else {
		fseek(fpa, sizeof(int)*2, SEEK_SET);
		fseek(fpb, sizeof(int)*2, SEEK_SET);
	}

	MPI_Bcast(&n, 1, MPI_INTEGER, 0, comm);

	int x_low = BLOCK_LOW(coords[0], root, n);
	int x_high = BLOCK_HIGH(coords[0], root, n);
	int y_low = BLOCK_LOW(coords[1], root, n);
	int y_high = BLOCK_HIGH(coords[1], root, n);

	int area = (x_high-x_low)*(y_high-y_low);
	double * mat_a = (double*)malloc(sizeof(double)*area);
	doulbe * mat_b = (double*)malloc(sizeof(double)*area);

	// X = (i/(upper_bound-lower_bound+1))+lower_bound
	// Y = (i%(upper_bound-lower_bound))+lower_bound 

	int x_off, y_off;
	for (int i = 0; i < area; ++i) {
		x_off = (i/(x_high-x_low+1))+x_low;
		y_off = (i%(x_high-x_low+1))+y_low; // ???
	}

}


// This algorithm is basically black magic
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
