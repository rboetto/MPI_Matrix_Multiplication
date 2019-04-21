double * read_matrix (File * fp, int * m, int * n) {
	
	fread(m, sizeof(int), 1, fp);
	fread(n, sizeof(int), 1, fp);
	
	int m_size = (*m)*(*n);
	double * matrix = (double *)calloc(m_size, sizeof(double));
	
	fread(matrix, sizeof(double), m_size, fp);

	return matrix;
	
}

int main(int argc, int * argv[]) {
	
	File fp_a = fopen(argv[1], "r");
	File fp_b = fopen(argv[2], "r");
	
	int m_a, n_a, m_b, n_b;
	double * matrix_a, * matrix_b;
	
	matrix_a = read_matrix(fp_a, &m_a, &n_a);
	matrix_b = read_matrix(fp_b, &m_b, &n_b);
	
	fclose(fp_a);
	fclose(fp_b);
	
	
	
	free(matrix_a);
	free(matrix_b);
	return 0;
}