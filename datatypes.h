#pragma once

#include "basetype.h"
/*
class VecR3 {
public:
	real x, y, z;

	VecR3();
	VecR3(real x,real y, real z);
	~VecR3();
};

class Atom {
public:
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
	real epsLJ;
	real RminLJ;
	VecR3 r;
	VecR3 v;
	VecR3 a;
	VecR3 f;

	Atom();
	Atom(char* name, char* symb, char* type, char* resi, long ires, char* segmreal,real mass, real chrg, real occp, real beta, real epsLl, real RminL3, VecR3 r, VecR3 v, VecR3 a, VecR3 f);
	~Atom();
};


class Pair {
public:
	int indi, indj;
	int indx;

	Pair();
	Pair(int indi, int indj, int indx);
	~Pair();
};

class Bond {
public:
	char atmi[8], atmj[8];
	int indi, indj;
	int ityp;
	int rigid;
	real b0;
	real Kb;
	real len;
	VecR3 unit;

	Bond();
	Bond(char* atomi, char* atomj, int indi, int indj, int ityp, int rigid, real b0, real Kb, real len, VecR3 unit);
	~Bond();
};*/

typedef struct VecR3 {
	real x, y, z;                                         // Cartesian components
} VecR3;

typedef struct Atom {                            // strings: 8 characters + NULL
	char name[8];                                                    // atom name
	char symb[8];                                              // chemical symbol
	char type[8];                                                    // atom type
	char resi[8];                                                 // residue name
	long ires;                                            // residue sequence no.
	char segm[8];                                                 // segment name
	real mass;                                                 // atom mass (amu)
	real chrg;                                                 // atom charge (e)
	real occp;                                                       // occupancy
	real beta;                                              // temperature factor
	real epsLJ;                // sqare root Lennard-Jones epsilon sqrt(kcal/mol)
	real RminLJ;                                      // Lennard-Jones Rmin/2 (A)
	VecR3 r;                                                      // position (A)
	VecR3 v;                                                   // velocity (A/ps)       
	VecR3 a;                                             // acceleration (A/ps^2)
	VecR3 f;                                               // forces (kcal/mol/A)
} Atom;

typedef struct Pair {
	int indi, indj;                            // indexes of objects forming pair
	int indx;                                          // index of pair in a list
} Pair;

typedef struct Bond {
	char atmi[8], atmj[8];                                          // atom types
	int indi, indj;                                               // atom indexes
	int ityp;                                                  // bond type index
	int rigid;                                                  // 1 - rigid bond
	real b0;                              // equilibrium value of bond length (A)
	real Kb;                   // bond stretching force constant (kcal/mole/A**2)
	real len;                                                      // bond length
	VecR3 unit;                                                    // unit vector
} Bond;