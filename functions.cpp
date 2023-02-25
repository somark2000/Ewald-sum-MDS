#include "functions.h"

//force field interactions
void Forces(Atom atoms[], int natm, real Rcut, real& ELJ) {
	/*
	ELJ - Lennard-Jones energy
	RluL - cutoff distance for short-range interactions
	Return: the non-bonded forces acting on atoms and the corresponding energy
	*/
	
	real fpr, Rcut2, Rr2, Rr6, rx, ry, rz, r2;
	real epsLJ, RminLJ;
	Rcut2 = Rcut * Rcut;
	for (int i = 1; i <= natm; i++) {
		atoms[i].f = VecR3();
	}

	//non-bonded interactions
	ELJ = 0;
	for (int i = 1; i < natm; i++) {
		for (int j = i + 1; j <= natm; j++) {
			rx = atoms[i].r.x - atoms[j].r.x; // interatomic distance
			ry = atoms[i].r.y - atoms[j].r.y;
			rz = atoms[i].r.z - atoms[j].r.z;
			fpr = 0e0; // inițializa f/r
			r2 = rx * rx + ry * ry + rz * rz; // squared interatomic distance
			if (r2 < Rcut2) { // short-range Lennard-Jones interactions
				// Lorentz-Berthelot mixing rules: RminLJ = Rmin / 2, epsLJ = sqrt(eps)
				RminLJ = atoms[i].RminLJ + atoms[j].RminLJ;
				epsLJ = atoms[i].epsLJ * atoms[j].epsLJ;
				Rr2 = RminLJ * RminLJ / r2;
				Rr6 = Rr2 * Rr2 * Rr2;
				ELJ += epsLJ * Rr6 * (Rr6 - 2);
				//fpr = 12e0 * epsLJ * Rr6 * (Rr6 - l) / r2; H f / r
			}
			if (fpr) {
				atoms[i].f.x += fpr * rx; atoms[j].f.x -= fpr * rx;
				atoms[i].f.y += fpr * ry; atoms[j].f.y -= fpr * ry;
				atoms[i].f.z += fpr * rz; atoms[j].f.z -= fpr * rz;
			}
		}
	}
}

//time propagation
void Verlet(Atom atoms[], int natm, real Rcut, real dt, real& ELJ, real& Ekin) {

}
