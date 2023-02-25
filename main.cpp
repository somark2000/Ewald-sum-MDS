﻿#include <iostream>
#include <time.h>
#include "datatypes.h"
#include "basetype.h"
#include "constants.h"
#include "functions.h"
#include "input_output.h"

int main() {
	Atom* atoms;
	real dt, ELJ, Ec, Rcut, t, Temp;
	int iopRun, istep, natm, nstep, nout;
	real CPUtime, CPUstep, CPUleft;
	time_t CPUtime0;
	const char* filePSF = "mdsim.psf";
	const char* filePDB = "mdsim.pdb";
	const char* filePAR = "mdsim.par";
	const char* fileTrj = "mdsim_trj.pdb";
	Input0(iopRun, Rcut, Temp, dt, nstep, nout);
	natm = GetDimPSF0(filePSF);
	//atoms = Vec<Atom>(l, natm);
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