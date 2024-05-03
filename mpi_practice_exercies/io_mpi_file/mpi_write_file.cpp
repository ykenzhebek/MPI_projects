#include<iostream>
#include<mpi.h>
#define n 10
using namespace std;
int main(int argc, char** argv) {

	int rank, size, i;
	char a[n];
	MPI_File fh;
	MPI_Status st;
	MPI_Request req;
	MPI_Offset displace;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	for (i = 0; i < n; i++)
		a[i] = (char)('0' + rank);

	displace = rank * n * sizeof(char);
	MPI_File_open(MPI_COMM_WORLD, "input.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_File_set_view(fh, displace, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
	MPI_File_write(fh, a, n, MPI_CHAR, &st);
	MPI_File_close(&fh);
	MPI_Finalize();
	return 0;
}
