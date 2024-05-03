#include<iostream>
#include<mpi.h>
using namespace std;

/*
	This code implements a ring communication pattern among MPI processes. 
	Each process sends a message to its right neighbor, receives a message from its left neighbor, 
	and updates the received value by adding its own rank to it. The updated value is then passed to the next process in the ring.

	Here's how the code works:

		Each process determines its left and right neighbors in the ring topology based on its rank.
		Inside the loop, each process sends its current sum (initialized to 0) to its right neighbor using MPI_Send.
		Simultaneously, each process starts a non-blocking receive operation (MPI_Irecv) to receive a message from its left neighbor.
		After completing the receive operation (MPI_Wait), each process updates the received value by adding its own rank to it.
		The updated value becomes the new sum for the current process.
		This process continues for a total of "size" iterations, where "size" is the number of MPI processes.
		Finally, each process prints its rank and the computed sum.

	The purpose of this code is to demonstrate a simple ring communication pattern using MPI_Send, MPI_Irecv, MPI_Wait,
	and non-blocking communication to achieve communication between neighboring processes in a ring topology.
*/

int main(int argc,char *argv[])
{
	int rank,size,recvbuf,i,sum;
	int left_neighbor, right_neighbor;
	int tag=12;
	MPI_Status statSend, statRecv;
	MPI_Request reqSend, reqRecv;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	left_neighbor = (rank-1 + size)%size;
	right_neighbor = (rank+1)%size;
	
	sum=0;
	for(i=1;i<=size;i++){
		
			MPI_Send(&sum,1,MPI_INT,right_neighbor,tag,MPI_COMM_WORLD);
			MPI_Irecv(&recvbuf,1,MPI_INT,left_neighbor,tag,MPI_COMM_WORLD,&reqRecv);

			MPI_Wait(&reqRecv,&statRecv);

			recvbuf=recvbuf+rank;
			sum=recvbuf;
		}
			printf("Among %d processes, process %d summa: %d\n", size, rank, sum);
		
	
	MPI_Finalize();
	return 0;
}