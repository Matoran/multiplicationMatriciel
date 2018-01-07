all: fox dns cannon
fox: multmat
	mpirun -np 9 ./multmat 1 0 3 > results.dat
	octave-cli check_multmat.m
cannon: multmat
	mpirun -np 9 ./multmat 2 0 3 > results.dat
	octave-cli check_multmat.m
dns: multmat
	mpirun -np 27 ./multmat 3 0 3 > results.dat
	octave-cli check_multmat.m

multmat: multmat.cc
	mpiCC multmat.cc -o multmat
