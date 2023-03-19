#pragma once
#ifndef functions_h //assure that this file is defined just one time
#define functions_h

#include <iostream>
#include "basetype.h"
#include "datatypes.h"
#include "constants.h"

// Initializes the random number generator using the current system time
#define seed() srand((unsigned)time(NULL))
// Generates random floating point numbers in the range [0,1)
#define random() ((real)rand() / (RAND_MAX + 1))

//functions to be used
void Forces(Atom atoms[], int n, real Rcut, real& ELJ);
void Forces1(Atom atoms[], int natm, int PBC, VecR3 Box, real Rcut, real& ELJ, real& virial);
void Verlet(Atom atoms[], int n, real Rcut, real dt, real& ELJ, real& Ec);
void Verlet1(Atom atoms[], int natm, int PBC, VecR3 Box, real Rcut, real dt, real& ELJ, real& Ekin, real& virial);
void ZeroVelAcc(Atom atoms[], int n);
void Vrescal(Atom atoms[], int n, real Ec);
void Heat(Atom atoms[], int n, real Ec);
void Thermostat(Atom atoms[], int natm, real Temp, real tauT, real dt, real& Ekin);

#endif // !functions_h
