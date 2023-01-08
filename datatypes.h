#pragma once

#include "basetype.h"

class VecR3 {
private:
	real x, y, z;
public:
	VecR3();
	VecR3(real x,real y, real z);
	~VecR3();
};

class Atom {
private:
	char name[8];
	char symb[8];
	char type[8];
	char resi[8];
	long ires;
	char segm[8];
	real mass;
	real chrg;
	real occp;
	real beta;
	real epsLl;
	real RminL3;
	VecR3 r;
	VecR3 v;
	VecR3 a;
	VecR3 f;
public:
	Atom(char* name, char* symb, char* type, char* resi, long ires, char* segmreal,real mass, real chrg, real occp, real beta, real epsLl, real RminL3, VecR3 r, VecR3 v, VecR3 a, VecR3 f);
	~Atom();
};