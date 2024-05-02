#include<iostream>
#include<stdio.h>
#include<mpi.h>
#include<fstream>
#include<math.h>
#include"2Dcoefficients1Dec.h"
#include"2DparamUnkn1Dec.h"
#include<time.h>

#define time 100
#define N 1024
#define GsX 4
#define GsY 2

using namespace std;
int ProcRank = 0;
int ProcNum = 0;
int GridSizeX, GridSizeY;
MPI_Comm GridComm;
MPI_Comm ColComm;
MPI_Comm RowComm;
int GridCoords[2];

void CreateGridCommunicators() {
	int DimSize[2];
	int Periodic[2];
	int Subdims[2];
	DimSize[0] = GridSizeY;
	DimSize[1] = GridSizeX;
	Periodic[0] = 1;
	Periodic[1] = 1;
	MPI_Dims_create(ProcNum, 2, DimSize);
	MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 1, &GridComm);
	MPI_Cart_coords(GridComm, ProcRank, 2, GridCoords);
	Subdims[0] = 0;
	Subdims[1] = 1;
	MPI_Cart_sub(GridComm, Subdims, &RowComm);
	Subdims[0] = 1;
	Subdims[1] = 0;
	MPI_Cart_sub(GridComm, Subdims, &ColComm);
}


void CheckerboardMatrixScatter(double* pMatrix, double* pMatrixBlock, int Size, int BlockSizex, int BLockSizey) {
	double * pMatrixRow = new double[BLockSizey*Size];
	if (GridCoords[1] == 0) {
		MPI_Scatter(pMatrix, BLockSizey*Size, MPI_DOUBLE, pMatrixRow, BLockSizey*Size, MPI_DOUBLE, 0, ColComm);

	}
	for (int i = 0; i<BLockSizey; i++) {
		MPI_Scatter(&pMatrixRow[i*Size], BlockSizex, MPI_DOUBLE, &(pMatrixBlock[i*BlockSizex]), BlockSizex, MPI_DOUBLE, 0, RowComm);
	}
	delete[] pMatrixRow;
}

void ResultCollection(double* pCMatrix, double* pCblock, int Size, int BlockSizex, int BLockSizey) {
	double * pResultRow = new double[Size*BLockSizey];
	for (int i = 0; i < BLockSizey; i++) {
		MPI_Gather(&pCblock[i*BlockSizex], BlockSizex, MPI_DOUBLE, &pResultRow[i*Size], BlockSizex, MPI_DOUBLE, 0, RowComm);
	}
	if (GridCoords[1] == 0) {
		MPI_Gather(pResultRow, BLockSizey*Size, MPI_DOUBLE, pCMatrix, BLockSizey*Size, MPI_DOUBLE, 0, ColComm);
	}
	delete[] pResultRow;
}

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
}

int main(int argc, char** argv) {
	double start, end;
	int i, j, n, k;
	int iter = 0;
	int BlockSizeX, BlockSizeY;
	n = N;
	
	double *M = new double[n*n];

	MPI_Init(&argc, &argv);
	start = MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Status st;
	
	GridSizeX = GsX; 
	GridSizeY = GsY;
	CreateGridCommunicators();   
	BlockSizeX = n / GridSizeX;   
	BlockSizeY = n / GridSizeY; 

	double h = 1.0 / (n - 1);
	double T = 1.0;
	double dt = T / time;
	double V = dt / (h*h);
	double top, down, tk;
	double A, B, C;
	double *U, *U1, *U2;
	U = new double[BlockSizeX*BlockSizeY]; 
	U1 = new double[BlockSizeX*BlockSizeY]; 
	U2 = new double[BlockSizeX*BlockSizeY]; 
	double *x = new double[n + 2];
	double *y = new double[n + 2];
	double *x1 = new double[BlockSizeY + 2];
	double *y1 = new double[BlockSizeY + 2];
	double *v = new double[BlockSizeY + 1];
	double *z = new double[BlockSizeY + 1];
	double *w = new double[BlockSizeY + 1];
	double *Fx = new double[BlockSizeY + 1];
	double *vprevG = new double[ProcNum];
	double *zprevG = new double[ProcNum];
	double *wprevG = new double[ProcNum];
	double *v1G = new double[ProcNum];
	double *z1G = new double[ProcNum];
	double *w1G = new double[ProcNum];
	double *FpartG = new double[ProcNum];
	double *xp_part = new double[ProcNum + 1];
	double *xp = new double[ProcNum + 1];
	double *Ap = new double[GridSizeX + 1];
	double *Bp = new double[GridSizeX + 1];
	double *Cp = new double[GridSizeX + 1];
	double *Fp = new double[GridSizeX + 1];
	double *Apy = new double[GridSizeY + 1];
	double *Bpy = new double[GridSizeY + 1];
	double *Cpy = new double[GridSizeY + 1];
	double *Fpy = new double[GridSizeY + 1];
	double f = -3;
	
	for (i = 0; i <= n; i++)
	{
		A = 0.5*V;
		B = 1 + V;
		C = 0.5*V;
	}

	if (ProcRank == 0) {
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				x[i] = i * h;
				y[j] = j * h;

			}
		}
	}
	
	if (ProcRank == 0) {
		for (i = 1; i < n - 1; i++) {
			for (j = 1; j < n - 1; j++) {
				M[i*n + j] = x[i] * x[i] + y[j] * y[j];
			}
		}
		//---------------boundary conditions------------//	
		for (i = 0; i < n; i++) {
			M[i*n + 0] = 0.0;
			M[i*n + (n - 1)] = 0.0;
		}

		for (j = 1; j < n; j++) {
			M[j + n * 0] = 0.0;
			M[j + n * (n - 1)] = 0.0;
		}
		
		/*cout << " M = " << endl;
		for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
		printf("%.12f\t", M[i*n+ j]);
		}
		cout << endl;
		}*/
	}
	
	CheckerboardMatrixScatter(M, U, n, BlockSizeX,BlockSizeY);

												  /* if (ProcRank == 3) {
												   cout << endl << "Mblock = " << endl;
												   for (i = 0; i < BlockSizeY; i++) {
												   for (j = 0; j < BlockSizeX; j++) {
												   printf("%.12f\t", U[i*BlockSizeX + j]);
												   }
												   cout << endl;
												   }
												   }*/

	
	if (ProcRank == 0) {
		MPI_Send(y, n, MPI_DOUBLE, GridSizeX - 1, 77, MPI_COMM_WORLD);
	}

	if (ProcRank == GridSizeX - 1) {
		MPI_Recv(y, n, MPI_DOUBLE, 0, 77, MPI_COMM_WORLD, &st);
	}

	if (ProcRank == 0) {
		MPI_Send(x, n, MPI_DOUBLE, GridSizeX*(GridSizeY - 1), 87, MPI_COMM_WORLD);
	}
	if (ProcRank == GridSizeX*(GridSizeY - 1)) {
		MPI_Recv(x, n, MPI_DOUBLE, 0, 87, MPI_COMM_WORLD, &st);
	}
	
	/*if (ProcRank == GridSizeX-1){//*(GridSizeY-1)) {
	cout << "xr=" << endl;
	for (i = 0; i < n; i++) {
	printf("%.12f\t", y[i]);
	}
	}*/
	
	if ((GridCoords[1] == 0) || (GridCoords[1] == GridSizeX - 1)) {
		MPI_Scatter(y, BlockSizeY, MPI_DOUBLE, y1, BlockSizeY, MPI_DOUBLE, 0, ColComm);
	}
	
	if ((GridCoords[0] == 0) || (GridCoords[0] == GridSizeY - 1)) {
		MPI_Scatter(x, BlockSizeX, MPI_DOUBLE, x1, BlockSizeX, MPI_DOUBLE, 0, RowComm);
	}
	
	double timer1 = 0, timer2 = 0, timer3 = 0, timer4 = 0, timer5 = 0, timer6 = 0;
	//---computing time layer---//
	for (int t = 1; t <= time; t++) {
		tk = t*dt;
		if (GridCoords[1] == 0) {
			for (i = 0; i < BlockSizeY; i++) {
				U[i*BlockSizeX + 0] = y1[i] * y1[i] + tk;
				U1[i*BlockSizeX + 0] = y1[i] * y1[i] + tk;
			}
		}

		if (GridCoords[1] == GridSizeX - 1) {
			for (i = 0; i < BlockSizeY; i++) {
				U[i*BlockSizeX + (BlockSizeX - 1)] = 1 + y1[i] * y1[i] + tk;
				U1[i*BlockSizeX + (BlockSizeX - 1)] = 1 + y1[i] * y1[i] + tk;
			}
		}

		if (GridCoords[0] == 0) {
			for (j = 0; j < BlockSizeX; j++) {
				U[0 * BlockSizeX + j] = x1[j] * x1[j] + tk; //po x
				U1[0 * BlockSizeX + j] = x1[j] * x1[j] + tk; //po x
			}
		}

		if (GridCoords[0] == GridSizeY - 1) {
			for (j = 0; j < BlockSizeX; j++) {
				U[BlockSizeX * (BlockSizeY - 1) + j] = 1 + x1[j] * x1[j] + tk; //po x
				U1[BlockSizeX * (BlockSizeY - 1) + j] = 1 + x1[j] * x1[j] + tk; //po x
			}
		}

		/*if (ProcRank == 4) {
		for (i = 0; i < BlockSizeY; i++) {
		cout << "y1 = " << y1[i] << endl;
		}
		}*/

		/*for (i = 0; i < BlockSize; i++) {
		U[i*BlockSize + 0] = y1[i] * y1[i] + tk;
		U[i*BlockSize + (BlockSize - 1)] = 1 + y1[i] * y1[i] + tk;
		U1[i*BlockSize + 0] = y1[i] * y1[i] + tk;
		U1[i*BlockSize + (BlockSize - 1)] = 1 + y1[i] * y1[i] + tk;
		}

		for (j = 0; j < BlockSize; j++) {
		U[0 * BlockSize + j] = x1[j] * x1[j] + tk; //po x
		U1[0 * BlockSize + j] = x1[j] * x1[j] + tk; //po x
		}

		for (j = 1; j < BlockSize - 1; j++) {
		U[BlockSize * (BlockSize - 1) + j] = 1 + x1[j] * x1[j] + tk; //po x
		U1[BlockSize * (BlockSize - 1) + j] = 1 + x1[j] * x1[j] + tk; //po x
		}*/

		/*if (ProcRank == 7) {
		cout << endl << "U = " << endl;
		for (i = 0; i < BlockSizeY; i++) {
		for (j = 0; j < BlockSizeX; j++) {
		printf("%.12f\t", U[i*BlockSizeX + j]);
		}
		cout << endl;
		}
		}*/


		for (i = 0; i < BlockSizeY; i++) {
			//--------------method Thomas----------------//
			double t1 = MPI_Wtime(); //time of method Thomas i

			v[0] = 1.0;
			v[BlockSizeX - 1] = 0.0;
			z[0] = 0.0;
			z[BlockSizeX - 1] = 1.0;
			w[0] = 0.0;
			w[BlockSizeX - 1] = 0.0;

			//-----find v,z,w coefficients-----//
			top = 1;
			down = 0;
			for (j = 0; j < BlockSizeX; j++) {
				Fx[j] = 0.0;
			}
			find_xRun1(v, BlockSizeX, n, A, B, C, Fx, top, ProcRank);

			top = 0;
			down = 1;
			for (j = 0; j < BlockSizeX; j++) {
				Fx[j] = 0.0;
			}
			find_xRun1(z, BlockSizeX, n, A, B, C, Fx, top, ProcRank);

			top = 0;
			down = 0;
			for (j = 1; j < BlockSizeX - 1; j++) {
				Fx[j] = 0.5*V*(U[i*BlockSizeX + (j - 1)] - 2 * U[i*BlockSizeX + j] + U[i*BlockSizeX + (j + 1)]) + U[i*BlockSizeX + j] + 0.5*dt*f;
				//cout << "Fx[1] = " << Fx[1];
			}
			find_xRun1(w, BlockSizeX, n, A, B, C, Fx, top, ProcRank);
			double t2 = MPI_Wtime(); // time of method Thomas
			timer1 = timer1 + (t2 - t1);
			
			/*if (ProcRank == 3) {
				cout << endl << "vzw = " << endl;
				//for (i = 0; i < BlockSizeY; i++) {
					for (j = 0; j < BlockSizeX; j++) {
						printf("%.12f\t", w[j]);
					}
					//cout << endl;
				//}
			}*/




			double t3 = MPI_Wtime(); //time second step of Thomas algorithm on i
									 //--------gathers v,z,w boundary coeffiecients----------------------------//
			MPI_Gather(&v[1], 1, MPI_DOUBLE, v1G, 1, MPI_DOUBLE, 0, RowComm);
			MPI_Gather(&v[BlockSizeX - 2], 1, MPI_DOUBLE, vprevG, 1, MPI_DOUBLE, 0, RowComm);
			MPI_Gather(&z[1], 1, MPI_DOUBLE, z1G, 1, MPI_DOUBLE, 0, RowComm);
			MPI_Gather(&z[BlockSizeX - 2], 1, MPI_DOUBLE, zprevG, 1, MPI_DOUBLE, 0, RowComm);
			MPI_Gather(&w[1], 1, MPI_DOUBLE, w1G, 1, MPI_DOUBLE, 0, RowComm);
			MPI_Gather(&w[BlockSizeX - 2], 1, MPI_DOUBLE, wprevG, 1, MPI_DOUBLE, 0, RowComm);
			MPI_Gather(&Fx[1], 1, MPI_DOUBLE, FpartG, 1, MPI_DOUBLE, 0, RowComm);

			/*if (ProcRank == 0) {
			cout << "Gathers :" << endl;
			for (j = 0; j < GridSizeX; j++) {
			printf("%.12f\t", v1G[j]);
			}
			cout << endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);*/
			
			
			//----------find parametric unknowns-----------//
			if (GridCoords[1] == 0) {
				

				for (j = 1; j <= GridSizeX - 1; j++) {
					Ap[j] = A*vprevG[j];
					Bp[j] = B - A*zprevG[j] - C*v1G[j];
					Cp[j] = C*z1G[j];
					Fp[j] = FpartG[j] + A*wprevG[j - 1] + C*w1G[j];
				}

				xp_part[0] = y1[i] * y1[i] + tk;
				xp_part[GridSizeX] = 1 + y1[i] * y1[i] + tk;
				top = xp_part[0];
				find_xp(xp_part, BlockSizeX, n, Ap, Bp, Cp, Fp, top, GridSizeX);

				
			}

		
			//---------Scatterv parametric xp for each process and computing other unknowns---------//
			int *displs = new int[GridSizeX];
			int *sendcount = new int[GridSizeX];
			int offset = 0;
			for (j = 0; j < GridSizeX; j++) {
				displs[j] = offset;
				sendcount[j] = 2;
				offset++;
			}

			MPI_Scatterv(xp_part, sendcount, displs, MPI_DOUBLE, xp, 2, MPI_DOUBLE, 0, RowComm);
			double t4 = MPI_Wtime(); //second 2 step Thomas
			timer2 = timer2 + (t4 - t3);

			
			double t5 = MPI_Wtime(); //time of 3 step Thomas
			for (j = 0; j < BlockSizeX; j++) {
				U1[i*BlockSizeX + j] = v[j] * xp[0] + z[j] * xp[1] + w[j];
			}
			double t6 = MPI_Wtime();
			timer3 = timer3 + (t6 - t5);

		}
		/*MPI_Barrier(MPI_COMM_WORLD);
		if (ProcRank == 0) {
			cout << "M=" << endl;
			for (i = 0; i < BlockSizeY; i++) {
				for (j = 0; j < BlockSizeX; j++)
					printf("%.12f\t", U1[i*BlockSizeX + j]);
				cout << endl;
			}
		}*/
	
	//ResultCollection(M, U1, n, BlockSizeX, BlockSizeY);
	/*MPI_Barrier(MPI_COMM_WORLD);
	if (ProcRank == 0) {
		cout << "M1=" << endl;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++)
				printf("%.12f\t", M[i*n + j]);
			cout << endl;
		}
	}*/

	for (j = 0; j < BlockSizeX; j++) {
		//--------------Thomas algorithm----------------//
		double t7 = MPI_Wtime(); //time of Thomas

		v[0] = 1.0;
		v[BlockSizeY - 1] = 0.0;
		z[0] = 0.0;
		z[BlockSizeY - 1] = 1.0;
		w[0] = 0.0;
		w[BlockSizeY - 1] = 0.0;

		//-----find v,z,w coefficients-----//
		top = 1;
		down = 0;
		for (i = 0; i < BlockSizeY; i++) {
			Fx[i] = 0.0;
		}
		find_xRun1(v, BlockSizeY, n, A, B, C, Fx, top, ProcRank);

		top = 0;
		down = 1;
		for (i = 0; i < BlockSizeY; i++) {
			Fx[i] = 0.0;
		}
		find_xRun1(z, BlockSizeY, n, A, B, C, Fx, top, ProcRank);

		top = 0;
		down = 0;
		for (i = 1; i < BlockSizeY - 1; i++) {
			Fx[i] = 0.5*V*(U1[(i - 1)*BlockSizeX + j] - 2 * U1[i*BlockSizeX + j] + U1[(i + 1)*BlockSizeX + j]) + U1[i*BlockSizeX + j] + 0.5*dt*f;
		}
		find_xRun1(w, BlockSizeY, n, A, B, C, Fx, top, ProcRank);
		double t8 = MPI_Wtime(); // vremya metoda Yanenko
		timer4 = timer4 + (t8 - t7);

		/*if (ProcRank == 3) {
			cout << endl << "vzw = " << endl;
			//for (i = 0; i < BlockSizeY; i++) {
			for (i = 0; i < BlockSizeY; i++) {
				printf("%.12f\t", w[i]);
			}
			//cout << endl;
			//}
		}*/

		double t9 = MPI_Wtime(); //time second step Method Yanenko
								 //--------gathers v,z,w boundary coeffiecients----------------------------//
		MPI_Gather(&v[1], 1, MPI_DOUBLE, v1G, 1, MPI_DOUBLE, 0, ColComm);
		MPI_Gather(&v[BlockSizeY - 2], 1, MPI_DOUBLE, vprevG, 1, MPI_DOUBLE, 0, ColComm);
		MPI_Gather(&z[1], 1, MPI_DOUBLE, z1G, 1, MPI_DOUBLE, 0, ColComm);
		MPI_Gather(&z[BlockSizeY - 2], 1, MPI_DOUBLE, zprevG, 1, MPI_DOUBLE, 0, ColComm);
		MPI_Gather(&w[1], 1, MPI_DOUBLE, w1G, 1, MPI_DOUBLE, 0, ColComm);
		MPI_Gather(&w[BlockSizeY - 2], 1, MPI_DOUBLE, wprevG, 1, MPI_DOUBLE, 0, ColComm);
		MPI_Gather(&Fx[1], 1, MPI_DOUBLE, FpartG, 1, MPI_DOUBLE, 0, ColComm);

				//----------find parametric unknowns-----------//
		if (GridCoords[0] == 0) {
			

			for (i = 1; i <= GridSizeY - 1; i++) {
				Apy[i] = A*vprevG[i];
				Bpy[i] = B - A*zprevG[i] - C*v1G[i];
				Cpy[i] = C*z1G[i];
				Fpy[i] = FpartG[i] + A*wprevG[i - 1] + C*w1G[i];
			}

			xp_part[0] = x1[j] * x1[j] + tk;
			xp_part[GridSizeY] = 1 + x1[j] * x1[j] + tk;
			top = xp_part[0];
			find_xp(xp_part, BlockSizeY, n, Apy, Bpy, Cpy, Fpy, top, GridSizeY);

		
		}
	
			//---------Scatterv parametric xp for each process and computing other unknowns---------//
			int *displs = new int[GridSizeY];
			int *sendcount = new int[GridSizeY];
			int offset = 0;
			for (i = 0; i < GridSizeY; i++) {
				displs[i] = offset;
				sendcount[i] = 2;
				offset++;
			}
			MPI_Scatterv(xp_part, sendcount, displs, MPI_DOUBLE, xp, 2, MPI_DOUBLE, 0, ColComm);
			double t10 = MPI_Wtime(); 
			timer5 = timer5 + (t10 - t9);
		

			double t11 = MPI_Wtime();
			for (i = 0; i < BlockSizeY; i++) {
				U2[i*BlockSizeX + j] = v[i] * xp[0] + z[i] * xp[1] + w[i];
			}
			double t12 = MPI_Wtime();
			timer6 = timer6 + (t12 - t11);

		}

		for (i = 0; i < BlockSizeY; i++) {
			for (j = 0; j < BlockSizeX; j++) {
				U[i*BlockSizeX + j] = U2[i*BlockSizeX + j];
			}
		}

		iter++;
		if (ProcRank == 0) { cout << "\r" << iter; }
	}//time

	 //--------------------Result Collection and print in root processor--------------------//
	ResultCollection(M, U2, n, BlockSizeX, BlockSizeY);

	end = MPI_Wtime();

	/*MPI_Barrier(MPI_COMM_WORLD);
	if (ProcRank == 0) {
		cout << "M2=" << endl;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++)
				printf("%.12f\t", M[i*n + j]);
			cout << endl;
		}
	}*/
	
	if (ProcRank == 0) {
		printf("\n for NxM : %d", n);
		printf("\n Alltime: %f", end - start);
		printf("\n XStep_1_Thomas_Time: %f", timer1);
		printf("\n XStep_2_Thomas_Time: %f", timer2);
		printf("\n XStep_3_Thomas_Time: %f", timer3);
		printf("\n X_time_Thomas: %f", (timer1 + timer2 + timer3));
		printf("\n");
		printf("\n YStep1YanenkoTime: %f", timer4);
		printf("\n YStep2YanenkoTime: %f", timer5);
		printf("\n YStep3YanenkoTime: %f", timer6);
		printf("\n Y_time_Thomas: %f", (timer4 + timer5 + timer6));
		cout << endl << " iteration = " << iter;
	}
	
	MPI_Finalize();
	
	delete[] M;
	delete[] U;
	delete[] U1;
	delete[] U2;
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
	delete[] Apy;
	delete[] Bpy;
	delete[] Cpy;
	delete[] Fpy;

	return 0;
}