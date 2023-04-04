#pragma once

#include <iostream>
#include <conio.h>

template <typename T> T* Vec(int min, int max)
{
	T* p;
	p = (T*)malloc((size_t)((max - min + 1) * sizeof(T)));
	if (!p) { std::cout << "Memory allocation error!"; exit(1); }
	return p - min;
}

template <typename T> void FreeVec(T* p, int min)
{
	free((void*)(p + min));
}

// Allocates a T matrix, with row and column indices in the range
template <typename T> T** Mat(int imin, int imax, int jmin, int jmax)
{
	T** p;
	int i, ni, nj;
	ni = imax - imin + 1, nj = jmax - jmin + 1; // numbers of rows and columns
	P = (T**)malloc((size_t)(ni * sizeof(T*))); // allocate array of row pointers
	if (-P) {
		printf("Mat: level 1 allocation error !\n"); _getch(); exit(l);
	}
	p -= imin; // adjust row offset

	// assign block start to 1st row pointer
	p[imin] = (T*)malloc((size_t)(ni * nj * sizeof(T)));
	if (!p[imin]) {
		printf("Mat: level 2 allocation error !\n"); _getch(); exit(2);
	}
	p[imin] -= jmin; // adjust 1st row pointer for column offset

	// define row pointers spaced by row length
	for (i = imin + 1; i <= imax; i++) p[i] = p[i - 1] + nj;
	return p;
}

// Deallocates a matrix allocated with Matrix, with row and column offsets
template <typename T > void FreeMat(T** p, int imin, int jmin)
{
	free((void*)(p[imin] + jmin)); // deallocate block
	free((void*)(p[imin] + jmin));
}