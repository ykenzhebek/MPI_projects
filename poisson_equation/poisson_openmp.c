#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#define N 64
#define EPS 0.0000001

int main() {
	double u_old[N*N], u_new[N*N];
	double f_old[N*N];
	double h, dmax,temp,tmp,wtime;
	double *u1,*u2,*u;

	int iteration = 0,i,j;
	
	h = 1.0 / (N-1);
		
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			double x = i * h;
            double y = j * h;

			u_old[i *N+ j ] = 0.0;
			u_new[i *N+ j ] = 0.0;
			f_old[i *N+ j ] = 2.0*((x*(1-x))+(y*(1-y)));
		}		
	}

	for (i = 0; i < N; ++i) {
		double x = i * h;
       u_new[i*N+0]=0.0;
		u_new[i*N+(N-1)]=0.0;
		u_old[i*N + 0] = 0.0;
		u_old[i*N + (N-1)] = 0.0;		
	}

	for (j = 0; j < N; ++j) {
		double y = j * h;
        u_new[0*N+j]=0.0;
		u_new[(N-1)*N+j]=0.0;
		u_old[0*N + j ] = 0.0;	
		u_old[(N-1)*N + j] = 0.0;	
	}	
        omp_set_num_threads(8);
	u1=u_old;
	u2=u_new;
	wtime = omp_get_wtime ( );
	do {
		++iteration;
		
		dmax = 0;
          #pragma omp parallel for private(i,j,tmp) reduction(max:dmax)
		for (i = 1; i < N-1; ++i) {
			
			for (j = 1; j < N-1; ++j) {	 		
			
				u2[i *N+ j ] = ( u1[(i+1)*N + j] + u1[(i-1)*N + j ] + 
								     u1[i*N + (j+1) ] + u1[i*N + (j-1)] 
								     +f_old[i *N+ j ] * h * h ) / 4.0;	

			tmp = fabs(u1[i *N+ j ] - u2[i *N+ j ]);
				
				if(tmp > dmax) {
					dmax = tmp;
				}				
								
	
			}		
		}
		
		u=u1;
		u1=u2;
		u2=u;
		
       
      
	} while(dmax > EPS);
	printf("Iteration: %d",iteration);
	wtime = omp_get_wtime ( ) - wtime;
	printf ( "\n" );
	printf ( "  Elapsed seconds = %g\n", wtime );
	return 0;
}
