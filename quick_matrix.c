#include <stdio.h>
#include <stdlib.h>

int main (int argc, char * argv[]) {

FILE * fp = fopen("new_matrix", "w");

int n = 4;
fwrite(&n, sizeof(int), 1, fp);
fwrite(&n, sizeof(int), 1, fp);

double d;
int i;
for (i = 0; i < 16; i++) {
	d = 10*(double)rand()/RAND_MAX;
	fwrite(&d, sizeof(double), 1, fp);
}


}
