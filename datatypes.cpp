#include "datatypes.h"
#include <string.h>
#include <iostream>

VecR3::VecR3()
{
	
}

VecR3::VecR3(real x1, real y1, real z1) {
	x = x1;
	y = y1;
	z = z1;
}

VecR3::~VecR3()
{
}

Atom::Atom(char* name1, char* symb1, char* type1, char* resi1, long ires1, char* segm1, real mass1, real chrg1, real occp1, real beta1, real epsLl1, real RminL31, VecR3 r1, VecR3 v1, VecR3 a1, VecR3 f1)
{
	strcpy(name, name1);
	strcpy(symb, symb1);
	strcpy(type, type1);
	strcpy(resi, resi1);
	strcpy(segm, segm1);
	ires = ires1;
	mass = mass1;



}

Atom::~Atom()
{
}
