# MPI_Ass_5
Matrix multiply

Uses message-passing interface to multiply two N by N matrices

$ make
to compile multiplication program

$ make testbench
to compile auxilliary programs and create two matrix files "input1" and "input2"

$ make run
to run with input files named "input1" and "input2". Output file will be "output".

$ make clean
to executables as well as "input1", "input2", and "output".

Matrix file specification:
1 int N describing first dimension |
1 int M describing second dimension |
N*M doubles |

Note that the multiplier only handles two N by N matrices.
