#include<math.h>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<cmath>
#include<time.h>


#define N 2048
#define M 2048
#define time 100
#define Pi 3.14
using namespace std;
/*
string convertInt(int number)
{
	if (number == 0)
		return "0";
	string temp = "";
	string returnvalue = "";
	while (number>0) {
		temp += number % 10 + 48;
		number /= 10;
	}
	for (int i = 0; i<temp.length(); i++)
		returnvalue += temp[temp.length() - i - 1];
	return returnvalue;
}*/


int main() {
	clock_t t1 = clock();
	int iteration = 0;
	int i, j, k, l, s, t;
	double h = 1.0 / (N - 1);
	double tk;
	double T = 1.0;
	double dt = T / time;
	double *U = new double[(N + 1)*(M + 1)];
	double *Upr = new double[(N + 1)*(M + 1)];
	double *U1 = new double[(N + 1)*(M + 1)];
	double *Utoch = new double[(N + 1)*(M + 1)];
	double f = -3;
	double *F = new double[M + 2];
	double *alphax = new double[N + 2];
	double *alphay = new double[M + 2];
	double *betax = new double[N + 2];
	double *betay = new double[M + 2];
	double *A = new double[M + 2];
	double *B = new double[M + 2];
	double *C = new double[M + 2];
	double *x = new double[M + 2];
	double *y = new double[M + 2];
	double v = dt / (h*h);

	for (i = 1; i <= N; i++) {
		A[i] = 0.5*v;
		B[i] = 1 + v;
		C[i] = 0.5*v;
	}

	//---set x, y and initial conditions-----------------//
	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++) {
			x[i] = i * h;
			y[j] = j * h;
	}
	}
	//---initial conditions --//
	for (i = 1; i < N-1; i++) {
		for (j = 1; j < M-1; j++) {
			U[i*N + j] = x[i]*x[i] + y[j]*y[j];
		}
	}
	
	/*cout << "Nachalnie U:=" << endl;
	for (i = 1; i < N-1; i++) {
		for (j = 1; j < M-1; j++) {
			printf("%.12f\t", U[i*N + j]);//cout << Upr[i][j] << " ";
		}
		cout << endl;
	}*/

	float timer1 = 0, timer2 = 0;
		for (t = 1; t <= time; t++) {

			//---boundary conditions---//
			tk = t*dt;
			
		for (j = 0; j < M; j++) {
				U[0 * N + j] = x[j] * x[j] + tk;
				U1[0 * N + j] = x[j] * x[j] + tk;
				Upr[0 * N + j] = x[j] * x[j] + tk;
				//Utoch[0 * N + j] = x[j] * x[j] + tk; 

				U[(M - 1)*N + j] = 1 + x[j] * x[j] + tk; 
				U1[(M - 1)*N + j] = 1 + x[j] * x[j] + tk; 
				Upr[(M - 1)*N + j] = 1 + x[j] * x[j] + tk;
				//Utoch[(M - 1)*N + j] = 1 + x[j] * x[j] + tk; 
			}
			//-left and right --//
			for (i = 0; i < N; i++) {
				U[i*N + 0] = y[i] * y[i] + tk;
				U1[i*N + 0] = y[i] * y[i] + tk;
				Upr[i*N + 0] = y[i] * y[i] + tk;
				//Utoch[i*N + 0] = y[i] * y[i] + tk;

				U[i*N + (N - 1)] = 1 + y[i] * y[i] + tk;
				U1[i*N + (N - 1)] = 1 + y[i] * y[i] + tk;
				Upr[i*N + (N - 1)] = 1 + y[i] * y[i] + tk;
				//Utoch[i*N + (N - 1)] = 1 + y[i] * y[i] + tk;
			}
			
			//----sweep by x intermediate(k+0.5) layer----//
			clock_t t3 = clock();
			for (i = 1; i <= N - 2; i++) {
				//---from the boundary conditions we find alpha and beta---//
				alphax[1] = 0;
				betax[1] = U[i*N + 0];
				for (j = 1; j <= M - 2; j++) {
					//---coefficients of sweep---
					F[j] = 0.5*v*U[(i-1)*N + j] - v*U[i*N + j] + 0.5*v*U[(i+1)*N + j] + U[i*N + j] + 0.5*dt*f;
					
					//---find alpha and beta on x-- forward step
					alphax[j + 1] = (C[j]) / ((B[j])-(A[j]*alphax[j]));
					betax[j + 1] = (F[j] + A[j]*betax[j]) / ((B[j])-(A[j]*alphax[j]));
				}
				//---back step sweep---//
				for (j = M - 2; j >= 0; j--) {

					Upr[i*N + j] = alphax[j + 1] * Upr[i*N + (j + 1)] + betax[j + 1];

				}
			}
			clock_t t4 = clock();
			timer1 = timer1 + (t4 - t3);
			/*cout << "Upr = " << endl;
			for (i = 0; i < N; i++) {
				for (j = 0; j < M; j++) {
					printf("%.12f\t", Upr[i*N + j]);//cout << Upr[i][j] << " ";
				}
				cout << endl;
			}*/

			clock_t t5 = clock();
			//----sweep by x intermediate(k+1) layer----//
			for (j = 1; j <= M - 2; j++) {
				//---from the boundary conditions we find alpha and beta---//
				alphay[1] = 0;
				betay[1] = Upr[0 * N + j];//tk;
				for (i = 1; i <= N - 2; i++) {
					//---coefficients of progonka---
					F[i] = 0.5*v*Upr[i*N + (j-1)] - v*Upr[i*N + j] + 0.5*v*Upr[i*N + (j+1)] + Upr[i*N + j] + 0.5*dt*f;
					//---find alpha and beta on x----forward step
					alphay[i + 1] = (C[i]) / ((B[i])-(A[i]*alphay[i]));
					betay[i + 1] = (F[i] + A[i]*betay[i]) / ((B[i])-(A[i]*alphay[i]));
				}
				
				//---back step sweep---//
				for (i = N - 2; i >= 0; i--) {

					U1[i*N + j] = alphay[i + 1] * U1[(i+1)*N + j] + betay[i + 1];

				}
			}
			clock_t t6 = clock();
			timer2 = timer2 + (t6 - t5);

			for (i = 1; i < N - 1; i++) {
				for (j = 1; j < M - 1; j++) {
					U[i*N + j] = U1[i*N + j];
				}
			}

			cout << "\r" << iteration;
			iteration++;
			
			/*//cout << "Tochnoe rewenie = " << endl;
			for (i = 1; i < N - 1; i++) {
				for (j = 1; j < M - 1; j++) {
					Utoch[i*N + j] = x[i] * x[i] + y[j] * y[j] + tk;
				}
				cout << endl;
			}*/
			/*cout << "U1 = " << endl;
			for (i = 0; i < N; i++) {
			for (j = 0; j < M; j++) {
			printf("%.12f\t", U1[i*N + j]);
			}
			cout << endl;
			}*/
			/*cout << "U1 = " << endl;
			for (i = 0; i < N; i++) {
				for (j = 0; j < M; j++) {
					printf("%.12f\t", U1[i*N + j]);//cout << U1[i][j] << " ";
				}
				cout << endl;
			}

			//cout << "Tochnoe rewenie = " << endl;
			for (i = 1; i < N-1; i++) {
				for (j = 1; j < M-1; j++) {
					Utoch[i*N + j] = x[i] * x[i] + y[j] * y[j] + tk;
				}
				cout << endl;
			}

			cout << "Tochnoe rewenie = " << endl;
			for (i = 0; i < N; i++) {
				for (j = 0; j < M; j++) {
					printf("%.12f\t", Utoch[i*N + j]);
				}
				cout << endl;
			}*/

			/*if (t % 10 == 0) {
				FILE * outputfile;
				string fname1 = "heatSeq2DImp";
				string pre = "x";
				string suf = ".dat";
				string outname;
				outname = convertInt(N) + pre + convertInt(t) + fname1 + suf;
				outputfile = fopen(outname.c_str(), "w");
			
				fprintf(outputfile, "TITLE=\"USERData\"\nVARIABLES=i, j, U\n");
				fprintf(outputfile, "ZONE T=\"ZONE1\", i=%d j=%d f=Point\n", N, M);

				for (int i = 0; i < N; i++) {
					for (int j = 0; j < M; j++) {
						fprintf(outputfile, "%.8f\t%.8f\t%.8f\t\n",
							i*h, j*h, U1[i*N +j]);
					}
				}
				fclose(outputfile);

			}
			*/

		}

		/*cout << "U1 = " << endl;
		for (i = 0; i < N; i++) {
			for (j = 0; j < M; j++) {
				printf("%.12f\t", U1[i*N + j]);
			}
			cout << endl;
		}*/

		
		clock_t t2 = clock();
		printf("\n for NxM = %d", N);
		printf("\n Total_time : %f", (float)(t2 - t1) / CLOCKS_PER_SEC);
		printf("\n Step_1_Time : %f", (float)(timer1) / CLOCKS_PER_SEC);
		printf("\n Step_2_Time : %f", (float)(timer2) / CLOCKS_PER_SEC);

		delete[] alphax;
		delete[] alphay;
		delete[] betax;
		delete[] betay;
		delete[] U;
		delete[] Upr;
		delete[] U1;
		delete[] F;
		delete[] A;
		delete[] B;
		delete[] C;

	return 0;
	
}
