#pragma once
#ifndef input_output_h //assure that this file is defined just one time
#define input_output_h

#include <iostream>
#include "basetype.h"
#include "datatypes.h"
#include "constants.h"

void Forces0(Atom atoms[], int natm, real Rcut, real& ELJ);
void Verlet0(Atom atoms[], int natm, real Rcut, real dt, real& ELJ, real& Ekin);

#endif // !input_output_h