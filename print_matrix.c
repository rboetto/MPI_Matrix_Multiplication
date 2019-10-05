#include <stdio.h>

int main (int argc, char * argv[]) {

FILE * fp = fopen(argv[1], "r");
int n[2];
double f;
fread(&n, sizeof(int), 2, fp);
printf("%d %d\n", n[0], n[1]);

int i;
for (i = 0; i < 16; i++) {
	fread(&f, sizeof(double), 1, fp);
	printf("%f ", f);
	if (i == 3 || i == 7 || i == 11 || i == 15) {
	printf("\n");
	}

}


}
