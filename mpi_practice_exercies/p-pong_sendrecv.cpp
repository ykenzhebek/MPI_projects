#include<iostream>
#include<mpi.h>
using namespace std;

/*
	This code performs a "ping-pong" task, where MPI processes exchange messages in a loop.
	One process sends a message to the other, then receives a response message, and so on.
	This is called a "ping-pong game," where each process takes turns sending and receiving messages.

	The task for this code is to demonstrate the MPI_Sendrecv function, which simultaneously sendsand receives messages between two MPI processes.
	Each process sends its rank to the other processand receives the rank from it.The processes continue to exchange messages in a loop.

	The purpose of this task is to demonstrate the usage of the MPI_Sendrecv function for message exchange between MPI processesand measure the execution time of this operation.
*/

int main(int argc,char *argv[])
{
	int rank,size,recvrank,i,sd,tag=12;
	double start,end;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Status status;

	if(size!=2){
		cout<<"error!";
		MPI_Abort(MPI_COMM_WORLD,1);
	}
		
	if(rank==0)
		{sd=1;}
	else sd=0;

	start = MPI_Wtime();
	for(i=1;i<=100;i++)
	{
		MPI_Sendrecv(&rank, 1, MPI_INT, sd, tag, &recvrank, 1, MPI_INT, sd, tag, MPI_COMM_WORLD, &status);
		
		if(rank==1)
			cout<<"i "<<i<<" rank="<<rank<<" receive "<<recvrank<<endl;
		else 
			cout<<"i "<<i<<" rank="<<rank<<" receive "<<recvrank<<endl;
	}
	end=MPI_Wtime();
	
	if(rank==0){cout<<"time: "<<end-start;}
	
	MPI_Finalize();
	return 0;
}