CC=mpiCC
#CC=mpiCC -O3
freedos.out: grid.o multigrid.o mesh.o main.o
	$(CC) main.o grid.o multigrid.o mesh.o -o freedos.out -lboost_program_options

grid.o: grid.cpp grid.h
	$(CC) -c grid.cpp 
multigrid.o: multigrid.cpp multigrid.h
	$(CC) -c multigrid.cpp 
mesh.o: mesh.cpp mesh.h
	$(CC) -c mesh.cpp 
main.o: main.cpp multigrid.o mesh.o grid.o 
	$(CC) -c main.cpp 
clean:
	rm *.o *.out

