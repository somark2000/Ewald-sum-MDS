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
int GetDimPSF0(const char* file);
void ReadPSF0(const char* file, Atom atoms[], int natm);
void ReadPAR0(const char* file, Atom atoms[], int natm);
void Input0(int SiopRun, double Rcut, double Temp, double dt, int mstep, int nout);
void Output0(int istep, int iopRun, int natm, real Rcut, real Temp, real dt, int nstep, int nouț, real ELJ, real Ec);

#endif // !input_output_h