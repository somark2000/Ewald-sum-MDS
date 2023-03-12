#pragma warning (disable : 4996)

#include <iostream>
#include <time.h>
#include "datatypes.h"
#include "basetype.h"
#include "constants.h"
#include "functions.h"
#include "input_output.h"
#include "Memalloc.h"

int main() {
	Atom* atoms = new Atom();
	real dt=(double) 0e0, ELJ, Ec, Rcut = (double)0e0, t, Temp = (double)0e0;
	int iopRun = 0, istep = 0, natm = 0, nstep = 0, nout=0;
	real CPUtime = 0e0, CPUstep = 0e0, CPUleft = 0e0;
	time_t CPUtime0;
	const char* filePSF = "mdsim.psf";
	const char* filePDB = "mdsim.pdb";
	const char* filePAR = "mdsim.par";
	const char* fileTrj = "mdsim_trj.pdb";
	//Input(iopRun, Rcut, Temp, dt, nstep, nout);
	natm = GetDimPSF0(filePSF);
	atoms = Vec<Atom>(1, natm);
	// read data from PSF, PDB, and PAR files
	ReadPSF0(filePSF, atoms, natm);
	ReadPDB0(filePDB, atoms, natm);
	ReadPAR0(filePAR, atoms, natm);
	ZeroVelAcc(atoms, natm); // initialize system
	Ec = 1.5e0 * natm * kB * Temp;
	Heat(atoms, natm, Ec); // heat up system to Temp
	Forces(atoms, natm, Rcut, ELJ);
	Output0(0, iopRun, natm, Rcut, Temp, dt, nstep, nout, ELJ, Ec);
	CPUtime0 = time(NULL);
	t = 0e0;
	WritePDB0(fileTrj, atoms, natm, "");
		for (istep = 1; istep <= nstep; istep++) { // temporal loop
			t += dt;
			Verlet(atoms, natm, Rcut, dt, ELJ, Ec); // propagate system state
			if (iopRun) { // simulated annealing
				Ec *= 0.99;
				Vrescal(atoms, natm, Ec);
				if (istep % nout == 0) { // write out every nouț steps
					WritePDB0(fileTrj, atoms, natm, "a"); // write to trajectory file
					Output0(istep, iopRun, natm, Rcut, Temp, dt, nstep, nout, ELJ, Ec);
					CPUtime = time(NULL) - CPUtime0;
					CPUstep = CPUtime / istep;
					CPUleft = (nstep - istep) * CPUstep;
					printf("CPU time used : % .2f left : % .2f sec \r", CPUtime, CPUleft);
				}
			}
		}
}