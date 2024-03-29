﻿#pragma warning (disable : 4996)

#include <iostream>
#include <time.h>
#include "datatypes.h"
#include "basetype.h"
#include "constants.h"
#include "functions.h"
#include "input_output.h"
#include "Memalloc.h"
#include "pairs.h"
#include "ewald.h"

int main() {
	Atom *atoms;
	Bond *bnds;
	Pair *pairs;
	VecR3 Box;
	real dt = (double)0e0, ELJ, Ec, Ebnd, Eele, Rcut = (double)0e0, t, Temp = (double)0e0, Tkin = (double)0e0, LBox = (double)0e0, Pres = (double)0e0, tauT = (double)0e0, Vbox = (double)0e0, virial = (double)0e0;
	int iopRun = 0, istep = 0, natm = 0, nstep = 0, nout=0, PBC=0, nbnd=0, npair=0;
	real CPUtime = 0e0, CPUstep = 0e0, CPUleft = 0e0;
	time_t CPUtime0;
	const char* filePSF = "C:/Users/smark/Documents/Uni/Licenta/Ewald-sum-MDS/mdsim.psf";
	//const char* filePSF = "C:/Users/smark/Documents/Uni/Licenta/Ewald-sum-MDS/CO_1000.psf";
	const char* filePDB = "C:/Users/smark/Documents/Uni/Licenta/Ewald-sum-MDS/mdsim.pdb";
	//const char* filePDB = "C:/Users/smark/Documents/Uni/Licenta/Ewald-sum-MDS/CO_1000.pdb";
	const char* filePAR = "C:/Users/smark/Documents/Uni/Licenta/Ewald-sum-MDS/mdsim.par";
	const char* fileTrj = "C:/Users/smark/Documents/Uni/Licenta/Ewald-sum-MDS/mdsim_trj.pdb";
	//Input(iopRun, Rcut, Temp, dt, nstep, nout);
	Input1(iopRun, PBC, LBox, Rcut, Temp, tauT, dt, nstep, nout);
	Box = { LBox, LBox, LBox };
	Vbox = Box.x * Box.y * Box.z;
	natm = GetDimPSF1(filePSF,nbnd);
	atoms = Vec<Atom>(1, natm);
	bnds = Vec<Bond>(1, nbnd);
	// read data from PSF, PDB, and PAR filesp
	ReadPSF1(filePSF, atoms, natm, bnds, nbnd);
	ReadPDB0(filePDB, atoms, natm);
	ReadPAR1(filePAR, atoms, natm, bnds, nbnd);
	npair = GetPairNum0(atoms, natm, bnds, nbnd);
	pairs = Vec<Pair>(1, npair);
	PairList0(atoms, natm, bnds, nbnd, pairs, npair);
	// initialize system
	ZeroVelAcc(atoms, natm);
	Ec = 1.5e0 * natm * kB * Temp;
	Heat(atoms, natm, Ec); // heat up system to Temp
	//Forces(atoms, natm, Rcut, ELJ);
	//Forces1(atoms, natm, PBC, Box, Rcut, ELJ, virial);
	//Forces2(atoms, natm, pairs, npair, bnds, nbnd, PBC, Box, Rcut, ELJ, Ebnd, virial);
	Forces3(atoms, natm, pairs, npair, bnds, nbnd, PBC, Box, Rcut, Ebnd, ELJ, Eele, virial);
	virial = 0e0;
	Pres = fpres * (2e0 * Ec + virial) / (3e0 * Vbox);
	//Output1(0, iopRun, natm, PBC,LBox, Rcut, Temp, tauT, dt, nstep, nout, ELJ, Ec, Temp, Pres,virial);
	//Output2(0, iopRun, natm, PBC, LBox, Rcut, Temp, tauT, dt, nstep, nout, Ebnd, ELJ, Ec, Temp, Pres, virial);
	Output3(0, iopRun, natm, PBC, LBox, Rcut, Temp, tauT, dt, nstep, nout, Ebnd, ELJ, Eele, Ec, Temp, Pres, virial);
	CPUtime0 = time(NULL);
	t = 0e0;
	WritePDB1(fileTrj, atoms, natm,Box, "w");
	for (istep = 1; istep <= nstep; istep++) { // temporal loop
		t += dt;
		//Verlet1(atoms, natm, PBC, Box, Rcut, dt, ELJ, Ec, virial); // propagate system state
		Verlet3(atoms, natm, pairs, npair, bnds, nbnd, PBC, Box, Rcut, dt, Ebnd, ELJ, Eele, Ec, virial);
		if (iopRun) { // simulated annealing
			Ec *= 0.99;
			Vrescal(atoms, natm, Ec);
		}

		Thermostat(atoms, natm, Temp, tauT, dt, Ec);
		Tkin = Ec/(1.5e0 * natm * kB);
		Pres = fpres * (2e0 * Ec + virial) / (3e0 * Vbox);

		if (istep % nout == 0) { // write out every nouț steps
			WritePDB1(fileTrj, atoms, natm, Box, "a"); // write to trajectory file
			//Output2(0, iopRun, natm, PBC, LBox, Rcut, Temp, tauT, dt, nstep, nout, Ebnd, ELJ, Ec, Temp, Pres, virial);
			Output3(istep, iopRun, natm, PBC, LBox, Rcut, Temp, tauT, dt, nstep, nout, Ebnd, ELJ, Eele, Ec, Tkin, Pres, virial);
			CPUtime = time(NULL) - CPUtime0;
			CPUstep = CPUtime / istep;
			CPUleft = (nstep - istep) * CPUstep;
			printf("CPU time used : % .2f left : % .2f sec \r", CPUtime, CPUleft);
		}
	}
}