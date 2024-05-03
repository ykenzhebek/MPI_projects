#include <iostream>
#include <math.h>
#include<mpi.h>

using namespace std;

/*
	This file implements the trapezoidal rule for numerical integration using MPI (Message Passing Interface) for parallelization.

	The task for this code is to approximate the definite integral of a given function 𝑓(𝑥) over the interval [𝑎,𝑏] using the trapezoidal rule with 𝑛 subintervals. 
	The code divides the interval [𝑎,𝑏] into 𝑛 subintervals and calculates the integral contribution for each process using the trapezoidal rule formula. 
	The results are then combined using MPI_Reduce to obtain the final integral value.

	The purpose of this task is to demonstrate parallel computation of numerical integration using the trapezoidal rule with MPI for distributing the workload among multiple processes.
*/

double f(double x)
{
	return sin(x);
}

int main(int argc, char** argv) {
	int n,k,rank,size,i;
	double a,b,y_part=0,y,h,In,starttime,endtime;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	a=0;
	b=3.14;
	n=1000;
	k=n/size;
	h=(b-a)/n;
	
	starttime=MPI_Wtime();
	
	if(rank==0) y_part+=f(a)+f(b);
	
	for (i=rank*k; i<rank*k+k; i++) {
		y_part+=2*(f(a+h*i));
	}
	
	if(rank==0) y_part-=2*f(a);
	
	MPI_Reduce(&y_part,&y,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	if(rank==0) {
		In=(b-a)/(2*n)*y;
		printf("%f\n",In);
	}
	endtime=MPI_Wtime();

	if(rank==0){cout<<"time: "<<endtime-starttime;}
	MPI_Finalize();
	
}