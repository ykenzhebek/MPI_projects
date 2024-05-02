#include<iostream>
#include<stdio.h>
#include<mpi.h>
#include<fstream>
#include<time.h>
#include<math.h>
#include"2Dcoefficients1Dec.h"
#include"2DparamUnkn1Dec.h"

#define time 100
#define N 512

using namespace std;

/*string convertInt(int number)
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
int main(int argc, char** argv) {
	double start, end;
	
	
	int i, j, n, k,t, ProcRank, ProcNum,iteration = 0, iterj = 0;
	n = N;

	MPI_Init(&argc, &argv);
	start = MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	k = (n / ProcNum);
	double h = 1.0 / (n - 1);
	double T = 1.0;
	double tk;
	double dt = T / time;
	double A, B, C, top, down;
	double *M = new double[n*n];
	double *U = new double[k*n];
	double *Upr = new double[k*n];
	double *U1 = new double[k*n];
	double *alphax = new double[n + 1];
	double *betax = new double[n + 1];
	double *F = new double[n + 1];
	double *x = new double[n + 2]; 
	double *y = new double[n + 2];
	double *x1 = new double[k + 2];
	double *y1 = new double[k + 2];
	double *v = new double[k + 1];
	double *z = new double[k + 1];
	double *w = new double[k + 1];
	double *Fx = new double[k + 1];
	double *vprevG = new double[ProcNum];
	double *zprevG = new double[ProcNum];
	double *wprevG = new double[ProcNum];
	double *v1G = new double[ProcNum];
	double *z1G = new double[ProcNum];
	double *w1G = new double[ProcNum];
	double *FpartG = new double[ProcNum];
	double *xp_part = new double[ProcNum + 1];
	double *xp = new double[ProcNum + 1];
	double *Ap = new double[ProcNum + 1];
	double *Bp = new double[ProcNum + 1];
	double *Cp = new double[ProcNum + 1];
	double *Fp = new double[ProcNum + 1];
	
	double f = -3;
	double V = dt / (h*h);

	A = 0.5*V;
	B = 1 + V;
	C = 0.5*V;
	
	if (ProcRank == 0) {
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				x[i] = i * h;
				y[j] = j * h;

			}
		}
	}

	if (ProcRank == ProcNum - 1) {
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				x[i] = i * h;
				y[j] = j * h;

			}
		}
	}
	
	//-------------M-------------------------------//
	if (ProcRank == 0) {
		for (i = 1; i < n-1; i++) {
			for (j = 1; j < n-1; j++) {
				M[i*n + j] = x[i] * x[i] + y[j] * y[j];
			}
		}
		//---------------boundary conditions-------------//	
		/*for (i = 0; i < n; i++) {
			M[i*n + 0] = 0;//y[i] * y[i] + tk;
			M[i*n + (n - 1)] = 0;// 1 + y[i] * y[i] + tk;
		}

		for (j = 1; j < n; j++) {
			M[j + n * 0] = 0;// x[j] * x[j] + tk;
			M[j + n * (n - 1)] = 0;// 1 + x[j] * x[j] + tk;
		}*/
		
		/*cout << "Nachalnie U:=" << endl;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				printf("%.12f\t", M[i*N + j]);
			}
			cout << endl;
		}*/

	}

	MPI_Scatter(M, k*n, MPI_DOUBLE, U, k*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(y, k, MPI_DOUBLE, y1, k, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/*if (ProcRank == 0) {
		cout << " Uscattered:=" << endl;
		for (i = 0; i < k; i++) {
			for (j = 0; j < n; j++) {
				printf("%f\t", U[i*n + j]);
			}
			cout << endl;
		}
	}*/

	double timer1 = 0, timer2 = 0, timer3 = 0, timer4 = 0;
	for (t = 1; t <= time; t++) {
		tk = t*dt;
		/*if (ProcRank == 0) {
			cout << "\r" << iteration<<":"<<1;
			//cout << "iterj = " << iterj;
		}*/

		//-- boundary conditions in each dinamic time layer --// 
		for (i = 0; i < k; i++) {
			U[i*n + 0] = y1[i] * y1[i] + tk;
			U[i*n + (n - 1)] = 1 + y1[i] * y1[i] + tk;
			U1[i*n + 0] = y1[i] * y1[i] + tk;
			U1[i*n + (n - 1)] = 1 + y1[i] * y1[i] + tk;
		}

		if (ProcRank == 0) {
			for (j = 0; j < n; j++) {
				U[0 * n + j] = x[j] * x[j] + tk; //po x
				U1[0 * n + j] = x[j] * x[j] + tk; //po x
			}
		}

		if (ProcRank == ProcNum - 1) {
			for (j = 1; j < n - 1; j++) {
				U[n * (k - 1) + j] = 1 + x[j] * x[j] + tk; //po x
				U1[n * (k - 1) + j] = 1 + x[j] * x[j] + tk; //po x
			}
		}


		/*if (ProcRank == 0) {
			cout << "Ubound=" << endl;
			for (i = 0; i < k; i++) {
				for (j = 0; j < n; j++) {
					printf("%.12f\t", U[i*n + j]);
				}
				cout << endl;
			}
		}*/
		double t1 = MPI_Wtime();
		for (i = 0; i < k; i++) {
			//---from the boundary conditions we find alpha and beta---//
			alphax[1] = 0;
			betax[1] = U[i*n + 0];
			for (j = 1; j <= n - 2; j++) {
				//---coefficients of progonka---(constanti)ne massiv//
				F[j] = 0.5*V * (U[i*n + (j - 1)] - 2 * U[i*n + j] + U[i*n + (j + 1)]) + U[i*n + j] + 0.5*dt*f;
				//---find alpha and beta po x-- pryamoy xod//
				alphax[j + 1] = (C) / ((B)-(A*alphax[j]));
				betax[j + 1] = (F[j] + A*betax[j]) / ((B)-(A*alphax[j]));
			}
			Upr[i*n + (n - 1)] = U[i*n + (n - 1)];
			//---obratniy xod progonka---//
			for (j = n - 2; j >= 0; j--) {

				Upr[i*n + j] = alphax[j + 1] * Upr[i*n + (j + 1)] + betax[j + 1];

			}
		}
		double t2 = MPI_Wtime();
		timer1 = timer1 + (t2 - t1);
		/*MPI_Barrier(MPI_COMM_WORLD);
		if (ProcRank == 0) {
			cout << "Upr=" << endl;
			for (i = 0; i < k; i++) {
				for (j = 0; j < n; j++) {
					printf("%.12f\t", Upr[i*n + j]);
				}
				cout << endl;
			}
		}*/

		
		/*if (ProcRank == 0) {
			cout << "\r" << iteration << ":" << 2;
			//cout << "iterj = " << iterj;
		}*/

		for (j = 1; j < n - 1; j++) {
			//--------------Thomas algo----------------//
			double t3 = MPI_Wtime(); //Time of Thomas algo

			

			v[0] = 1.0;
			v[k - 1] = 0.0;
			z[0] = 0.0;
			z[k - 1] = 1.0;
			w[0] = 0.0;
			w[k - 1] = 0.0;

			//-----find v,z,w coefficients-----//
			top = 1;
			down = 0;
			for (i = 0; i < k; i++) {
				Fx[i] = 0.0; 
			}
			find_xRun1(v, k, n, A, B, C, Fx, top, ProcRank);

			top = 0;
			down = 1;
			for (i = 0; i < k; i++) {
				Fx[i] = 0.0; 
			}
			find_xRun1(z, k, n, A, B, C, Fx, top, ProcRank);

			top = 0;
			down = 0;
			for (i = 1; i < k-1; i++) {
				Fx[i] = 0.5*V*(Upr[(i - 1)*n + j] - 2 * Upr[i*n + j] + Upr[(i + 1)*n + j])+ Upr[i*n + j] + 0.5*dt*f;
			}
			find_xRun1(w, k, n, A, B, C, Fx, top, ProcRank);

			double t4 = MPI_Wtime(); 
			timer2 = timer2 + (t4 - t3);

			double t5 = MPI_Wtime(); //time second step Thomas
			//--------gathers v,z,w boundary coeffiecients----------------------------//
			MPI_Gather(&v[1], 1, MPI_DOUBLE, v1G, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&v[k - 2], 1, MPI_DOUBLE, vprevG, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&z[1], 1, MPI_DOUBLE, z1G, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&z[k - 2], 1, MPI_DOUBLE, zprevG, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&w[1], 1, MPI_DOUBLE, w1G, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&w[k - 2], 1, MPI_DOUBLE, wprevG, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&Fx[1], 1, MPI_DOUBLE, FpartG, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			/*if (ProcRank == 0) {
				cout << "\r" << iteration << ":" << 3;
				//cout << "iterj = " << iterj;
			}*/

			//----------find parametric unknowns-----------//
			if (ProcRank == 0) {
				
				
				for (i = 1; i <= ProcNum - 1; i++) {
					Ap[i] = A*vprevG[i];
					Bp[i] = B - A*zprevG[i] - C*v1G[i];
					Cp[i] = C*z1G[i];
					Fp[i] = FpartG[i] + A*wprevG[i - 1] + C*w1G[i]; 
				}

				xp_part[0] = x[j] * x[j] + tk;
				xp_part[ProcNum] = 1 + x[j] * x[j] + tk;;
				top = xp_part[0];
				find_xp(xp_part, k, n, Ap, Bp, Cp, Fp, top, ProcNum);

				
			}
			/*if (ProcRank == 0) {
				cout << "\r" << iteration << ":" << 4;
				//cout << "iterj = " << iterj;
			}*/


			/*if (ProcRank == 0) {
				cout << "xp_part = " << endl;
				for (i = 0; i <= ProcNum; i++) {
					printf("%.12f\t", xp_part[i]);
				}
				cout << endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);*/

			//---------Bcast parametric xp for each process and computing other unknowns---------//
			MPI_Bcast(xp_part, ProcNum + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			double t6 = MPI_Wtime(); //second 2 step Yanenko
			timer3 = timer3 + (t6 - t5);

			double t7 = MPI_Wtime();
			for (i = 0; i < k; i++) {
				U1[i*n + j] = v[i] * xp_part[ProcRank] + z[i] * xp_part[ProcRank + 1] + w[i];
			}
			double t8 = MPI_Wtime();
			timer4 = timer4 + (t8 - t7);

			iterj++;
			
			
		}
		/*MPI_Barrier(MPI_COMM_WORLD);
		if (ProcRank == 0) {
			cout << "U1=" << endl;
			for (i = 0; i < k; i++) {
				for (j = 0; j < n; j++) {
					printf("%.12f\t", U1[i*n + j]);
				}
				cout << endl;
			}
		}*/

		/* MPI_Barrier(MPI_COMM_WORLD);
		if (ProcRank == 1) {
			cout << "U1:" << endl;
			for (j = 1; j < n - 1; j++) {
				for (i = 0; i < k; i++) {
					printf("%.12f\t", U1[i*n + j]);
				}
				cout << endl;
			}
			cout << endl;
		}*/
			//cout << "iterj = " << iterj;
		
		//------rewriting for next layer -----//
		for (i = 0; i < k; i++) {
			for (j = 0; j < n; j++) {
				U[i*n + j] = U1[i*n + j];
			}
		}

		
		/*MPI_Barrier(MPI_COMM_WORLD);
		if (ProcRank == 0) {
			cout << "M=" << endl;
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++)
					printf("%.12f\t", M[i*n + j]);
				cout << endl;
			}
		}*/
		
		iteration++;
		/*if (ProcRank == 0) {
			cout << "\r" << iteration;
			//cout << "iterj = " << iterj;
		}*/

		/*if (ProcRank == 0) {
			if (t % 10 == 0) {
				FILE * outputfile;
				string fname1 = "heatPar2DImp";
				string pre = "x";
				string suf = ".dat";
				string outname;
				outname = convertInt(N) + pre + convertInt(t) + fname1 + suf;
				outputfile = fopen(outname.c_str(), "w");
				//ofstream fout("poisson.txt");

				fprintf(outputfile, "TITLE=\"USERData\"\nVARIABLES=i, j, U\n");
				fprintf(outputfile, "ZONE T=\"ZONE1\", i=%d j=%d f=Point\n", n, n);

				for (int i = 0; i < n; i++) {
					for (int j = 0; j < n; j++) {
						fprintf(outputfile, "%.8f\t%.8f\t%.8f\t\n",
							i*h, j*h, M[i*n + j]);
					}
				}
				fclose(outputfile);

				
			}
		}*/

	}
	
	MPI_Gather(U1, k*n, MPI_DOUBLE, M, k*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	end = MPI_Wtime();

	/*MPI_Barrier(MPI_COMM_WORLD);
	if (ProcRank == 0) {
		cout << "M=" << endl;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++)
				printf("%.12f\t", M[i*n + j]);
			cout << endl;
		}
	}*/

	if (ProcRank == 0) {
		printf("\n for NxM : %d", N);
		printf("\n Alltime: %f", end - start);
		printf("\n Step_1_Time: %f", timer1);
		printf("\n Step_1_Thomas_Time: %f", timer2);
		printf("\n Step_2_Thomas_Time: %f", timer3);
		printf("\n Step_3_Thomas_Time: %f", timer4);
		printf("\n Total_time_of_Thomas: %f", (timer2+timer3+timer4));

		cout << endl << " iteration = " << iteration;
	}
	
	MPI_Finalize();

	delete[] M;
	delete[] alphax;
	delete[] betax;
	delete[] Upr;
	delete[] U;
	delete[] U1;
	delete[] F;
	delete[] x;
	delete[] y;
	delete[] x1;
	delete[] y1;
	delete[] v;
	delete[] z;
	delete[] w;
	delete[] Fx;
	delete[] vprevG;
	delete[] zprevG;
	delete[] wprevG;
	delete[] v1G;
	delete[] z1G;
	delete[] w1G;
	delete[] FpartG;
	delete[] xp_part;
	delete[] xp;
	delete[] Ap;
	delete[] Bp;
	delete[] Cp;
	delete[] Fp;
	

	return 0;
}
