cannon:
	mpicc ricardo_cannon.c

run:
	mpirun -np 4 ./a.out input input2 output

clean:
	rm *.out
