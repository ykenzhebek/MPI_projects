#pragma once
#include<stdio.h>
#include<iostream>
void find_xp(double *x, int k, int n, double *Ap, double *Bp, double *Cp, double *Fp, double topx, int size) {

	double *alphax = new double[size + 1];
	double *betax = new double[size + 1];

	alphax[0] = 0;
	betax[0] = topx;

		for (int j = 1; j <= size - 1; j++) {
			alphax[j] = (Cp[j]) / ((Bp[j])-(Ap[j] * alphax[j-1]));
			betax[j] = (Fp[j] + (Ap[j] * betax[j-1])) / ((Bp[j])-(Ap[j] * alphax[j-1]));
		}

		for (int j = size - 1; j >= 0; j--) 
		{
			x[j] = alphax[j] * x[j + 1] + betax[j];
		}

		delete[] alphax;
		delete[] betax;

}