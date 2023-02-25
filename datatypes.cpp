#include "datatypes.h"
#include <string.h>
#include <iostream>

VecR3::VecR3()
{
	x = 0;
	y = 0;
	z = 0;
}

VecR3::VecR3(real x1, real y1, real z1) {
	x = x1;
	y = y1;
	z = z1;
}

VecR3::~VecR3()
{
}

Atom::Atom()
{
	strcpy(name, "");
	strcpy(symb, "");
	strcpy(type, "");
	strcpy(resi, "");
	strcpy(segm, "");
	ires = 0;
	mass = 0;
	chrg = 0;
	occp = 0;
	beta = 0;
	epsLJ = 0;
	RminLJ = 0;
	r = VecR3();
	v = VecR3();
	a = VecR3();
	f = VecR3();
}

Atom::Atom(char* name1, char* symb1, char* type1, char* resi1, long ires1, char* segm1, real mass1, real chrg1, real occp1, real beta1, real epsLJ1, real RminLJ1, VecR3 r1, VecR3 v1, VecR3 a1, VecR3 f1)
{
	strcpy(name, name1);
	strcpy(symb, symb1);
	strcpy(type, type1);
	strcpy(resi, resi1);
	strcpy(segm, segm1);
	ires = ires1;
	mass = mass1;
	chrg = chrg1;
	occp = occp1;
	beta = beta1;
	epsLJ = epsLJ1;
	RminLJ = RminLJ1;
	r = r1;
	v = v1;
	a = a1;
	f = f1;
}

Atom::~Atom()
{
}
