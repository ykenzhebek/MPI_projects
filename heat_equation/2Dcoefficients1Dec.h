#pragma once
#include<stdio.h>
#include<iostream>
void find_xRun1(double *x, int k, int n, double A, double B, double C, double *F, double topx, int rank) {

	double *alphax = new double[k + 1];
	double *betax = new double[k + 1];

	alphax[0] = 0;
	betax[0] = topx;
	
		for (int j = 1; j < k-1; j++) {
			alphax[j] = (C) / ((B) - (A * alphax[j-1]));
			betax[j] = (F[j] + (A * betax[j-1])) / ((B) - (A * alphax[j-1])); 
																						
		}
		
		for (int j = k-2; j >= 0; j--)
		{
			x[j] = alphax[j] * x[j+1] + betax[j];
		}

	delete[] alphax;
	delete[] betax;
}

