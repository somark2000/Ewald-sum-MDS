#include "datatypes.h"
#include <string.h>
#include <iostream>

VecR3::VecR3()
{
	this->x = 0;
	this->y = 0;
	this->z = 0;
}

VecR3::VecR3(real x1, real y1, real z1) {
	this->x = x1;
	this->y = y1;
	this->z = z1;
}

VecR3::~VecR3()
{
}

Atom::Atom()
{
	strcpy_s(this->name, "");
	strcpy_s(this->symb, "");
	strcpy_s(this->type, "");
	strcpy_s(this->resi, "");
	strcpy_s(this->segm, "");
	this->ires = 0;
	this->mass = 0;
	this->chrg = 0;
	this->occp = 0;
	this->beta = 0;
	this->epsLJ = 0;
	this->RminLJ = 0;
	this->r = VecR3();
	this->v = VecR3();
	this->a = VecR3();
	this->f = VecR3();
}

Atom::Atom(char* name1, char* symb1, char* type1, char* resi1, long ires1, char* segm1, real mass1, real chrg1, real occp1, real beta1, real epsLJ1, real RminLJ1, VecR3 r1, VecR3 v1, VecR3 a1, VecR3 f1)
{
	strcpy_s(this->name, name1);
	strcpy_s(this->symb, symb1);
	strcpy_s(this->type, type1);
	strcpy_s(this->resi, resi1);
	strcpy_s(this->segm, segm1);
	this->ires = ires1;
	this->mass = mass1;
	this->chrg = chrg1;
	this->occp = occp1;
	this->beta = beta1;
	this->epsLJ = epsLJ1;
	this->RminLJ = RminLJ1;
	this->r = r1;
	this->v = v1;
	this->a = a1;
	this->f = f1;
}

Atom::~Atom()
{
}

Pair::Pair() {
	this->indi = 0;
	this->indj = 0;
	this->indx = 0;
}

Pair::Pair(int indi, int indj ,int indx) {
	this->indi = indi;
	this->indj = indj;
	this->indx = indx;
}

Pair::~Pair() {

}

Bond::Bond() {
	strcpy_s(this->atmi,"");
	strcpy_s(this->atmj,"");
	this->indi=0;
	this->indj = 0;
	this->ityp = 0;
	this->rigid = 0;
	this->b0 = 0;
	this->Kb = 0;
	this->len = 0;
	this->unit = VecR3();
}
Bond::Bond(char* atomi, char* atomj, int indi, int indj, int ityp, int rigid, real b0, real Kb, real len, VecR3 unit) {
	strcpy_s(this->atmi, atomi);
	strcpy_s(this->atmj, atomj);
	this->indi = indi;
	this->indj = indj;
	this->ityp = ityp;
	this->rigid = rigid;
	this->b0 = b0;
	this->Kb = Kb;
	this->len = len;
	this->unit = unit;
}
Bond::Bond() {

}