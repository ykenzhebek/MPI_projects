#include <mpi.h>
#include <iostream>
#include <fstream>

#define FILENAME "readfile.txt"

int main(int argc, char** argv) {
    int rank, size, ierr;
    MPI_File fh;
    MPI_Offset filesize;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Open the file
    MPI_File_open(MPI_COMM_WORLD, FILENAME, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    // Determine the file size
    MPI_File_get_size(fh, &filesize);

    // Allocate buffer for reading data
    char* buf = new char[filesize + 1];

    // Read the data
    MPI_File_read(fh, buf, filesize, MPI_CHAR, MPI_STATUS_IGNORE);
    buf[filesize] = '\0';

    // Print the data
    std::cout << "Process " << rank << " read: " << buf << std::endl;

    // Close the file
    MPI_File_close(&fh);

    delete[] buf;
    MPI_Finalize();
    return 0;
}