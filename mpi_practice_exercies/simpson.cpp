#include <iostream>
#include <math.h>
#include<mpi.h>
using namespace std;

/*
	This file implements Simpson's rule for numerical integration using MPI (Message Passing Interface) for parallelization.

	The purpose of this task is to demonstrate parallel computation of numerical integration using Simpson's rule with MPI for distributing the workload among multiple processes.
*/

double f(double x)
{
	return sin(x);
}

int main(int argc, char** argv) {
	
	int n,k,rank,size,i;
	double a,b,y_part=0,y,h,In,S1_part=0,S2_part=0,S1,S2,starttime,endtime;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	a=0;
	b=3.14;
	n=100000;
	k=n/size;
	h=(b-a)/n;

	starttime=MPI_Wtime();
	
	for (i=rank*k; i<rank*k+k; i++) {
		if(i%2==0){ S1_part=S1_part+f(a+h*i);}
		else S2_part=S2_part+f(a+h*i);
		
	}
	
	if(rank==0) S1_part-=f(a);
	MPI_Reduce(&S1_part,&S1,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&S2_part,&S2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	if(rank==0) {
		In=h/3*(f(a)+f(b)+2*S1+4*S2);
		printf("%f\n",In);
	}
	
	endtime=MPI_Wtime();
	
	if(rank==0){cout<<"time: "<<endtime-starttime;}
	
	MPI_Finalize();
	
}