#include <time.h>
#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 1024
#define EPS 0.0000001

/*
	1. u_old, u_new --> array
	2. f_old --> array
	3. eps --> double
	4. N --> double
	5. initial conditions
	6. boundary conditions
	...*/


int main(int argc, char **argv) {

	int iteration=0,i,j,size,rank;
	printf("test");	
	MPI_Init(&argc, &argv);
	MPI_Status st;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	int nthreads = 8;
	if (argc>1) sscanf(argv[1],"%d",&nthreads);
	MPI_Bcast(&nthreads,1,MPI_INT,0,MPI_COMM_WORLD);
	omp_set_num_threads(nthreads);

	double *u1,*u2,*u;

        double starttime,endtime;
	double *u_old_part=(double * ) malloc(((N/size)+2)*N * sizeof(double));
	double *f_old_part=(double * ) malloc(((N/size)+2)*N * sizeof(double));
	double *u_new_part=(double * ) malloc(((N/size)+2)*N * sizeof(double));
	double *u_old=(double * ) malloc(N*N * sizeof(double));
	double *f_old=(double * ) malloc(N*N * sizeof(double));
	double *u_new=(double * ) malloc(N*N * sizeof(double));
	
	double *tempup=(double * ) malloc(N * sizeof(double));
	double *tempdown=(double * ) malloc(N * sizeof(double));
	double *temp1=(double * ) malloc(N * sizeof(double));
	double *temp2=(double * ) malloc(N * sizeof(double));

	double  dmax,max,up,down;
	
	int *sendcount = (int * ) malloc(size * sizeof(int));
	int *displs = (int * ) malloc(size * sizeof(int));
	
	double h = 1.0 / (N-1);
	
	if(rank==0){
	 
		for (i = 0; i < N; ++i) {
			for (j = 0; j < N; ++j) {
				double x = i * h;
				double y = j * h;
				u_old[i*N + j] = 0.0;
				u_new[i*N + j] = 0.0;
				f_old[i*N + j ] =2.0*((x*(1-x))+(y*(1-y)));
			}		
	    }
	}
	
	for(i=1;i<size-1;i++){
		sendcount[i]=((N/size)+2)*N;
	}
	
	sendcount[0]=((N/size)+1)*N;
	sendcount[size-1]=((N/size)+1)*N;
	
	if(size==1){
		sendcount[0]=N*N;
	}

	displs[0]=0;
	for(i=1;i<size;i++){
		displs[i]=(displs[i-1]+sendcount[i-1])-2*N;
	}
	
	for(i=0;i<(N/size)+2;i++){
		for(j=0;j<N;j++) {
			u_new_part[i*N+j]=0.0;
		}
	}

	MPI_Scatterv(u_old,sendcount, displs, MPI_DOUBLE, u_old_part, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(f_old,sendcount, displs, MPI_DOUBLE, f_old_part, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if(rank==0){
		for(i=0;i<N; ++i)
		u_new_part[i]=0;
	}

	if(rank==size-1){
		for(i=0;i<N; ++i)
		u_new_part[(sendcount[rank]/N-1)*N+i]=0;
	}
	
	for(i=0;i<sendcount[rank]/N;i++){
		u_new_part[i*N]=0;
		u_new_part[i*N+N-1]=0;
	}

	#pragma omp parallel
	#pragma omp single
        printf("potok %d\n",omp_get_num_threads());
	
	starttime=MPI_Wtime();
	do {
		++iteration;
		dmax = 0;
		
        #pragma omp parallel for private(i,j) reduction(max:dmax)
		for ( i = 1; i <sendcount[rank]/N-1; ++i) {
			for ( j = 1; j < N-1; ++j) {				
				u_new_part[i*N + j ] = ( u_old_part[(i+1)*N + j ] + u_old_part[(i-1)*N + j] + 
								     u_old_part[i*N + (j+1)] + u_old_part[i*N + (j-1)] 
								     + f_old_part[i*N + j ] * h * h ) / 4.0;								
			if(fabs(u_old_part[i *N+ j ] - u_new_part[i *N+ j ]) >= dmax) {
					dmax = fabs(u_old_part[i *N+ j ] - u_new_part[i *N+ j ]);
				}	

			}		
		}
	
		
		if(rank!=0) up=rank-1;//для процесса с рангом 0 нет верхнего соседа
		else up=MPI_PROC_NULL;//верхний сосед
		if(rank!=size-1) down=rank+1;//для последнего процесса нет нижнего соседа
		else down=MPI_PROC_NULL;//нижний сосед

		if(rank!=0){
			for(i=0;i<N;i++){
				tempup[i]=u_new_part[i+N];
			}
		}
		
		if(rank!=size-1){
			for(i=0;i<N;i++){
				tempdown[i]=u_new_part[i+sendcount[rank]-2*N];
			}
		}
		
		MPI_Sendrecv(tempup,N,MPI_DOUBLE,up,12,temp1,N,MPI_DOUBLE,up,12,MPI_COMM_WORLD,&st);
		MPI_Sendrecv(tempdown,N,MPI_DOUBLE,down,12,temp2,N,MPI_DOUBLE,down,12,MPI_COMM_WORLD,&st);
		
		if(rank!=0){
			for(i=0;i<N;i++){
				u_new_part[i]=temp1[i];
			}
		}

		if(rank!=size-1){
			for(i=0;i<N;i++){
				u_new_part[i+sendcount[rank]-1*N]=temp2[i];
			}
		}
		
		
				for (i = 1; i < sendcount[rank]/N-1; i++) {
			for (j = 1; j < N-1; j++) {				
							
			}		
		}
		


		for (i = 0; i < sendcount[rank]/N; ++i) {
			for (j = 0; j < N; ++j) {						
					u_old_part[i *N + j ] = u_new_part[i *N + j ];				
			}		
		}
		MPI_Allreduce(&dmax,&max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);	
				
	} while(max > EPS);
	if(rank==0){
	endtime=MPI_Wtime();	
	}
	int *recvcount = (int * ) malloc(size * sizeof(int));
	int *recvdispls = (int * ) malloc(size * sizeof(int));
	
	for(i=0;i<size;i++)
		recvcount[i]=sendcount[i];
		recvdispls[0]=0;
	for(i=1;i<size;i++)
		recvdispls[i]=recvdispls[i-1]+recvcount[i-1]-2*N;
	MPI_Gatherv(u_new_part,sendcount[rank],MPI_DOUBLE,u_new,recvcount,recvdispls,MPI_DOUBLE,0,MPI_COMM_WORLD);	
	if(rank==0){
		//clock_t t2 = clock(); 
		printf("Time: %lf",(endtime - starttime) );  
	}

	if(rank==0) {
		printf("\nIteration : %d", iteration);
	}
		
	MPI_Finalize();
	
	free(u_old_part);
	free(f_old_part);
	free(u_new_part);

	free(u_old);
	free(f_old);
	free(u_new);
	
	free(recvcount);
	free(recvdispls);
	
	free(sendcount);
	free(displs);

	return 0;
}




