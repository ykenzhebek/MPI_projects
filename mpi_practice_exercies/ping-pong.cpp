#include<iostream>
#include<mpi.h>
using namespace std;

/*

	This code performs a similar "ping-pong" task as the previous one, but using MPI_Ssend and MPI_Recv functions for message passing between two MPI processes. 
	Each process takes turns sending and receiving messages in a loop.

	The task for this code is to demonstrate the usage of MPI_Ssend and MPI_Recv functions for synchronous message passing between two MPI processes. 
	Each process sends its rank to the other process and receives the rank from it, and then prints out the ranks along with the iteration number and tag.

	The purpose of this task is to compare the behavior and performance of synchronous message passing using MPI_Ssend and MPI_Recv functions with the previous implementation using MPI_Sendrecv.
*/

int main(int argc,char *argv[])
{
	int rank,size,recvrank,i;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	if(size!=2){
		cout<<"error!";
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	else
	for(i=1;i<=100;i++){
		if(rank==0 )
		{    
			MPI_Ssend(&rank,1,MPI_INT,1,0,MPI_COMM_WORLD);
			MPI_Recv(&recvrank,1,MPI_INT,1,1,MPI_COMM_WORLD,&Status);
			cout<<"i"<<" "<<i<<" "<<"I am rank"<<" "<<rank<<" "<<"which Recv rank "<<" "<<recvrank<<" "<<"tag 1"<<endl;
		}
		else  {
	
			MPI_Recv(&recvrank,1,MPI_INT,0,0,MPI_COMM_WORLD,&Status);
			cout<<"i"<<" "<<i<<" "<<"I am rank"<<" "<<rank<<" "<<"which Recv rank "<<" "<<recvrank<<" "<<"tag 0"<<endl;
			MPI_Ssend(&rank,1,MPI_INT,0,1,MPI_COMM_WORLD);
	
		}
	}
	MPI_Finalize();
	return 0;
}