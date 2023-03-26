#pragma once
#ifndef input_output_h //assure that this file is defined just one time
#define input_output_h

#include <iostream>
#include <fstream>
#include <conio.h>
#include <string.h>
#include "basetype.h"
#include "datatypes.h"
#include "constants.h"

void strstrip(char src[], char res[], int il, int i2);
int GetDimPDB(const char* file);
void ReadPDB0(const char* file, Atom atoms[], int natm);
void WritePDB0(const char* file, Atom atoms[], int natm, const char* mode);
void WritePDB1(const char* file, Atom atoms[], int natm, VecR3 Box, const char* mode);
int GetDimPSF0(const char* file);
int GetDimPSF1(const char* file, int ftnbnd);
void ReadPSF0(const char* file, Atom atoms[], int natm);
void ReadPSF1(const char* file, Atom atoms[], int natm, Bond bnds[], int nbnd);
void ReadPAR0(const char* file, Atom atoms[], int natm);
void ReadPAR1(const char* file, Atom atoms[], int natm, Bond bnds[], int nbnd);
void Input(int SiopRun, double Rcut, double Temp, real dt, int mstep, int nout);
void Input1(int& iopRun, int& PBC, real& Lbox, real& Rcut, real& Temp, real& tauT, real& dt, int& nstep, int& nout);
void Output0(int istep, int iopRun, int natm, real Rcut, real Temp, real dt, int nstep, int nouț, real ELJ, real Ec);
void Output1(int istep, int iopRun, int natm, int PBC, real Lbox, real Rcut, real Temp, real tauT, real dt, int nstep, int nout, real ELI, real Ekin, real Tkin, real Pres, real virial);
void Output2(int istep, int iopRun, int natm, int PBC, real Lbox, real Rcut, real Temp, real tauT, real dt, int nstep, int nout, real Ebnd, real ELI, real Ekin, real Tkin, real Pres, real virial);
#endif // !input_output_h