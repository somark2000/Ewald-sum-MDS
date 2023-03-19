#include "functions.h"

//force field interactions
void Forces(Atom atoms[], int n, real Rcut, real& ELJ) {
	/*
	ELJ - Lennard-Jones energy
	RluL - cutoff distance for short-range interactions
	Return: the non-bonded forces acting on atoms and the corresponding energy
	*/
	
	real fpr, Rcut2, Rr2, Rr6, rx, ry, rz, r2;
	real epsLJ, RminLJ;
	Rcut2 = Rcut * Rcut;
	for (int i = 1; i <= n; ++i) {
		atoms[i].f = VecR3();
	}

	//non-bonded interactions
	ELJ = 0e0;
	for (int i = 1; i < n; ++i) {
		for (int j = i + 1; j <= n; j++) {
			// interatomic distance
			rx = atoms[i].r.x - atoms[j].r.x;
			ry = atoms[i].r.y - atoms[j].r.y;
			rz = atoms[i].r.z - atoms[j].r.z;
			// f/r
			fpr = 0e0;
			r2 = rx * rx + ry * ry + rz * rz; // squared interatomic distance
			if (r2 < Rcut2) { 
				// short-range Lennard-Jones interactions
				// Lorentz-Berthelot mixing rules: RminLJ = Rmin / 2, epsLJ = sqrt(eps)
				RminLJ = atoms[i].RminLJ + atoms[j].RminLJ;
				epsLJ = atoms[i].epsLJ * atoms[j].epsLJ;
				Rr2 = RminLJ * RminLJ / r2;
				Rr6 = Rr2 * Rr2 * Rr2;
				ELJ += epsLJ * Rr6 * (Rr6 - 2e0);
				fpr = 12e0 * epsLJ * Rr6 * (Rr6 - 1e0) / r2;
			}
			if (fpr) {
				atoms[i].f.x += fpr * rx; atoms[j].f.x -= fpr * rx;
				atoms[i].f.y += fpr * ry; atoms[j].f.y -= fpr * ry;
				atoms[i].f.z += fpr * rz; atoms[j].f.z -= fpr * rz;
			}
		}
	}
}

// Returns the non-bonded forces acting on atoms and the corresponding energy.
void Forces1(Atom atoms[], int natm, int PBC, VecR3 Box, real Rcut, real& ELJ,	real &virial)
	// Rcut - cutoff distance for short-range interactions
	// ELJ - Lennard-Jones energy
{
	real Box2x, Box2y, Box2z;
	real fpr, Rcut2, Rr2, Rr6, rx, ry, rz, r2;
	real epsLJ, RminLJ;
	int iatm, jatm;
	Rcut2 = Rcut * Rcut;
	Box2x = Box.x / 2e0; // simulation box half-size
	Box2y = Box.y / 2e0;
	Box2z = Box.z / 2e0;
	for (iatm = 1; iatm <= natm; iatm++) atoms[iatm].f = { 0e0, 0e0, 0e0 };
	
	//Non-Bounded interaction
	ELJ = 0e0;
	virial = 0e0;
	for (iatm = 1; iatm < natm; iatm++) {
		for (jatm = iatm + 1; jatm <= natm; jatm++) {
			rx = atoms[iatm].r.x - atoms[jatm].r.x; // interatomic distance
			ry = atoms[iatm].r.y - atoms[jatm].r.y;
			rz = atoms[iatm].r.z - atoms[jatm].r.z;
			if (PBC) { // minimum image criterion
				if (rx > Box2x) rx -= Box.x; else if (rx < -Box2x) rx += Box.x;
				if (ry > Box2y) ry -= Box.y; else if (ry < -Box2y) ry += Box.y;
				if (rz > Box2z) rz -= Box.z; else if (rz < -Box2z) rz += Box.z;
				fpr = 0e0; // initialize f/r
				r2 = rx * rx + ry * ry + rz * rz; // squared interatomic distance
				if (r2 < Rcut2) { // short-range Lennard-Jones interactions
					// Lorentz-Berthelot mixing rules: RminLJ = Rmin / 2, epsLJ = sqrt(eps)
					RminLJ = atoms[iatm].RminLJ + atoms[jatm].RminLJ;
					epsLJ = atoms[iatm].epsLJ * atoms[jatm].epsLJ;
					Rr2 = RminLJ * RminLJ / r2;
					Rr6 = Rr2 * Rr2 * Rr2;
					ELJ += epsLJ * Rr6 * (Rr6 - 2e0);
					fpr = 12e0 * epsLJ * Rr6 * (Rr6 - 1e0) / r2; // f/r
				}
				if (fpr) {
					atoms[iatm].f.x += fpr * rx; atoms[jatm].f.x -= fpr * rx;
					atoms[iatm].f.y += fpr * ry; atoms[jatm].f.y -= fpr * ry;
					atoms[iatm].f.z += fpr * rz; atoms[jatm].f.z -= fpr * rz;
					virial += fpr * r2; // add virial term
				}
			}
		}
	}
}

//time propagation
void Verlet(Atom atoms[], int natm, real Rcut, real dt, real& ELJ, real& Ec) {
	/*
	MD propagator based on the velocity Verlet method.
	Rcut - cutoff distance for short-range interactions
	dt - time step
	ELJ - Lennard-Jones energy
	Ec - kinetic energy
	*/

	real dt2 = dt / 2e0;
	// PREDICTOR STEP
	for (int i = 1; i <= natm; ++i) {
		atoms[i].v.x += dt2 * atoms[i].a.x; // v(t+dt/2)
		atoms[i].v.y += dt2 * atoms[i].a.y;
		atoms[i].v.z += dt2 * atoms[i].a.z;
		atoms[i].r.x += dt * atoms[i].v.x; // r(t+dt)
		atoms[i].r.y += dt * atoms[i].v.y;
		atoms[i].r.z += dt * atoms[i].v.z;

		Forces(atoms, natm, Rcut, ELJ); // update forces

		// CORRECTOR STEP
		Ec = 0e0;
		for (i = 1; i <= natm; ++i) {
			atoms[i].a.x = facc * atoms[i].f.x / atoms[i].mass; // a(t+dt)
			atoms[i].a.y = facc * atoms[i].f.y / atoms[i].mass;
			atoms[i].a.z = facc * atoms[i].f.z / atoms[i].mass;
			atoms[i].v.x += dt2 * atoms[i].a.x; // v(t+dt)
			atoms[i].v.y += dt2 * atoms[i].a.y;
			atoms[i].v.z += dt2 * atoms[i].a.z;
			Ec += atoms[i].mass * (atoms[i].v.x * atoms[i].v.x + atoms[i].v.y * atoms[i].v.y + atoms[i].v.z * atoms[i].v.z);
			Ec *= fkin / 2e0;
		}
	}
}

void Verlet1(Atom atoms[], int natm, int PBC, VecR3 Box, real Rcut, real dt, real& ELJ, real& Ekin, real& virial)
	// MD propagator based on the velocity Verlet method.
	// Rcut - cutoff distance for short-range interactions
	// dt - time step
	// ELJ - Lennard-Jones energy
	// Ekin - kinetic energy
{
	real Box2x, Box2y, Box2z, dt2, v2;
	int iatm;
	dt2 = dt / 2e0;
	Box2x = Box.x / 2e0; // simulation box half-size
	Box2y = Box.y / 2e0;
	Box2z = Box.z / 2e0;
	for (iatm = 1; iatm <= natm; iatm++) { // PREDICTOR STEP
		//v(t+dt/2)
		atoms[iatm].v.x += dt2 * atoms[iatm].a.x;
		atoms[iatm].v.y += dt2 * atoms[iatm].a.y;
		atoms[iatm].v.z += dt2 * atoms[iatm].a.z;
		//r(t+dt)
		atoms[iatm].r.x += dt * atoms[iatm].v.x;
		atoms[iatm].r.y += dt * atoms[iatm].v.y;
		atoms[iatm].r.z += dt * atoms[iatm].v.z;
	}
	if (PBC) //wrap coordinates back into the original box
	{
		for (iatm = 1; iatm < natm; iatm++)
		{
			if (atoms[iatm].r.x > Box2x) atoms[iatm].r.x -= Box.x;
			else if (atoms[iatm].r.x < -Box2x) atoms[iatm].r.x += Box.x;
			if (atoms[iatm].r.y > Box2y) atoms[iatm].r.y -= Box.y;
			else if (atoms[iatm].r.y < -Box2y) atoms[iatm].r.y += Box.y;
			if (atoms[iatm].r.z > Box2z) atoms[iatm].r.z -= Box.z;
			else if (atoms[iatm].r.z < -Box2z) atoms[iatm].r.z += Box.z;
		}
	}
	
	//Update forces
	Forces1(atoms, natm, PBC, Box, Rcut, ELJ, virial);
	
	//Correction step
	Ekin = 0e0;
	for (iatm = 1; iatm < natm; iatm++) {
		atoms[iatm].a.x = facc * atoms[iatm].f.x / atoms[iatm].mass;
		atoms[iatm].a.y = facc * atoms[iatm].f.y / atoms[iatm].mass;
		atoms[iatm].a.z = facc * atoms[iatm].f.z / atoms[iatm].mass;

		atoms[iatm].v.x += dt2 * atoms[iatm].a.x;
		atoms[iatm].v.y += dt2 * atoms[iatm].a.y;
		atoms[iatm].v.z += dt2 * atoms[iatm].a.z;

		v2 = atoms[iatm].v.x * atoms[iatm].v.x +
			atoms[iatm].v.y * atoms[iatm].v.y +
			atoms[iatm].v.z * atoms[iatm].v.z;

		atoms[iatm].beta = v2; // squared velocity for velocity distribution
		Ekin += atoms[iatm].mass * v2;
	}
	Ekin *= fkin / 2e0;
}

//Resets the velocities and accelerations of all the atoms  to 0
void ZeroVelAcc(Atom atoms[], int n) {
	for (int i = 0; i < n; ++i) {
		atoms[i].v = VecR3();
		atoms[i].a = VecR3();
		atoms[i].f = VecR3();
	}
}

//Rescales velocities to match total kinetic energy
void Vrescal(Atom atoms[], int n, real Ec) {
	real Ec0, f;
	Ec0 = 0e0; // inițial kinetic energy
	for (int i = 1; i <= n; ++i) {
		Ec0 += atoms[i].mass * (atoms[i].v.x * atoms[i].v.x + atoms[i].v.y * atoms[i].v.y + atoms[i].v.z * atoms[i].v.z);
		Ec0 *= fkin / 2e0;
		f = sqrt(Ec / Ec0);
		for (i = 1; i <= n; i++) { // rescale velocities
			atoms[i].v.x *= f;
			atoms[i].v.y *= f;
			atoms[i].v.z *= f;
		}
	}
}

//Generates random velocities to match total kinetic energy
void Heat(Atom atoms[], int n, real Ec) {
	VecR3 vCM=VecR3();
	real mass, v, theta, phi;
	for (int i = 1; i <= n; i++) { // generate random velocities
		v = random();
		theta = random() * pi;
		phi = random() * pi * 2e0;
		atoms[i].v.x = v * sin(theta) * cos(phi);
		atoms[i].v.y = v * sin(theta) * sin(phi);
		atoms[i].v.z = v * cos(theta);
	}

	mass = 0e0;
	for (int i = 1; i <= n; i++) {
		vCM.x += atoms[i].mass * atoms[i].v.x;
		vCM.y += atoms[i].mass * atoms[i].v.y;
		vCM.z += atoms[i].mass * atoms[i].v.z;
		mass += atoms[i].mass;
	}

	vCM.x /= mass; vCM.y /= mass; vCM.z /= mass;
	for (int i = 1; i <= n; i++) { // subtract CM velocity
		atoms[i].v.x -= vCM.x;
		atoms[i].v.y -= vCM.y;
		atoms[i].v.z -= vCM.z;
		Vrescal(atoms, n, c);
	}
}


// Rescales velocities and total kinetic energy using the Berendsen thermostat
void Thermostat(Atom atoms[], int natm, real Temp, real tauT, real dt, real& Ekin)
{
	real f, Tkin;
	int iatm;
	Tkin = Ekin / (1.5e0 * natm * kB); // kinetic temperature
	f = sqrt(1e0 + (dt / tauT) * (Temp / Tkin - 1e0)); // scaling factor
	for (iatm = 1; iatm <= natm; iatm++) {
		atoms[iatm].v.x *= f;
		atoms[iatm].v.y *= f;
		atoms[iatm].v.z *= f;
	}
	Ekin *= f * f;
}
