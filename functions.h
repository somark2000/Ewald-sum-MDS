#pragma once
#ifndef functions_h //assure that this file is defined just one time
#define functions_h

#include <iostream>
#include "basetype.h"
#include "datatypes.h"
#include "constants.h"

void Forces(Atom atoms[], int natm, real Rcut, real& ELJ);
void Verlet(Atom atoms[], int natm, real Rcut, real dt, real& ELJ, real& Ekin);

#endif // !functions_h
