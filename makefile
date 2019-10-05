cannon:
	mpicc ricardo_cannon.c -o mult


run:
	mpirun -np 40 ./mult input1 input2 output


testbench:
	gcc quick_matrix.c -o quick_matrix
	gcc print_matrix.c -o print_matrix
	./quick_matrix
	mv new_matrix input1
	./quick_matrix
	mv new_matrix input2

clean:
	rm mult
	rm print_matrix
	rm quick_matrix
	rm input1
	rm input2
