#include<iostream>
#include<mpi.h>
using namespace std;

int main(int argc, char* argv[]) {
	MPI_Comm grid_comm;
	int dims[2];
	int periodic[2];
	int reorder = 1, q = 2, ndims = 2, maxdims = 2;
	int coordinates[2];
	int my_grid_rank;
	int coords[2];
	MPI_Comm row_comm;
	dims[0] = dims[1] = q;
	periodic[0] = periodic[1] = 1;
	coords[0] = 0; coords[1] = 1;
	MPI_Init(&argc, &argv);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic, reorder, &grid_comm);
	MPI_Comm_rank(grid_comm, &my_grid_rank);
	MPI_Cart_coords(grid_comm, my_grid_rank, maxdims, coordinates);
	printf("Process rank %i has coordinates %i %i\n", my_grid_rank, coordinates[0], coordinates[1]);
	MPI_Finalize();
	return 0;
}