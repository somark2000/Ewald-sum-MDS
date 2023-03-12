#pragma once

#include <iostream>
#include <conio.h>

template <typename T>
T* Vec(int min, int max)
{
	T* p;
	p = (T*)malloc((size_t)((max - min + 1) * sizeof(T)));
	if (!p) { std::cout << "Memory allocation error!"; exit(1); }
	return p - min;
}

template <typename T>
void FreeVec(T* p, int min)
{
	free((void*)(p + min));
}